require 'peach'
require 'Tools'

class PredictionsController < ApplicationController

	MAX_SIZE = 52428800 # 50 MB
	BASE_PATH = Dir::tmpdir + "/cug/" # all data files are stored in /tmp
	REF_DATA = "alignment_gene_structure.json"
	SORT = "/usr/bin/sort"
	FORMATDB = "/usr/bin/formatdb"
	BLASTALL = "/usr/bin/blastall"
	FASTACMD = "/usr/bin/fastacmd"
	AUGUSTUS = "/usr/local/bin/augustus/src/augustus"
	PAIR_ALIGN =  Rails.root.join('lib', 'pair_align').to_s #"/usr/local/bin/pair_align"
	DIALIGN2 = "/usr/local/bin/dialign_package/src/dialign2-2"
	# CLUSTAW = "/usr/bin/clustalw"
	MAFFT = "/usr/local/bin/mafft"
	HYDORPHOBIC_AAS = ["V", "I", "L"]
	POLAR_AAS = ["S", "T"]

	# Render start page for prediction
	def search
		delete_old_data # delete old uploaded genome files
		prepare_new_session # a fresh session
		bg_dataprocessing # align reference data, background job
	end

	# Handel uploaded genome file: 
	# store it in /tmp/cug and add its id(new file name), original filename and path to the session
	# render partial showing the uploaded file
	# accessible params: @fatal_error [Array] Errors occured during file load
	def upload_file
		@fatal_error = []
		is_success = system("fromdos", params[:uploaded_file].path)
		@fatal_error << "Fatal error." << "Cannot load genome data. Please contact us!" if ! is_success
		# is the file too big? maybe jquery-check did not work, so better check again
		@fatal_error |= check_filesize(params[:uploaded_file].size)
		# check content, only if file is not too big
		@fatal_error |= check_fasta(params[:uploaded_file]) if @fatal_error.blank?
		if @fatal_error.blank? then
			# rename and save file
			file_id = rand(1000000000).to_s
			file_name = File.basename(params[:uploaded_file].original_filename)
			file_path = BASE_PATH + file_id + "_query.fasta"
			File.rename(params[:uploaded_file].path, file_path)
			session[:file] = { id: file_id, name: file_name, path: file_path }
		end
		render :upload_file_ajax, formats: [:js]
	end

	# Handel loading the example genome:
	# store it in /tmp/cug and add its id(new file name), original filename and path to the session
	# render partial showing the uploaded file
	# accessible params: @fatal_error [Array] Errors occured during file load
	def load_example
		# do not check content, but check file size, just to be sure
		example_file = "Candida_albicans_WO_1.fasta"
		example_path = Rails.root + "spec/fixtures/files/" + example_file
		@fatal_error = check_filesize(File.size(example_path))
		if @fatal_error.blank? then
			file_id = rand(1000000000).to_s
			file_path = BASE_PATH + file_id + "_query.fasta"
			File.copy_stream(example_path, file_path)
			session[:file] = { id: file_id, name: example_file, path: file_path }
		end
		render :upload_file_ajax, formats: [:js]
	end

# uncomment if alignment options should get an own submit-button instead of the big "Predict button"
	# def set_alignment_options
	# 	session[:align] = { algo: params[:algo], config: params[:config] }
	# 	render :set_options
	# end

	# Prepare alignment file for view show_alignment
	# store it in /tmp/cymobase_alignment_ ... for lucullus
	# render partial show_alignment
	# accessible params in view: @frame_id [Array] Identify the right div
	# accessible params in view: @file_id [Array] Part of file path, for lucullus
	def show_alignment
		prot = params[:prot].gsub(" ", "-").downcase
		hit = params[:hit].to_s
		@frame_id = prot + "_" + hit
		# copy file to /tmp/cymo_alignment_
		@file_id = "cug" + session[:file][:id] + rand(1000000000).to_s
		# set correct file extension for dialign/ seqan files
		file_scr = BASE_PATH + session[:file][:id] + "_" + prot + "-" + hit.to_s + "-aligned.fasta"
		file_dest = Dir::tmpdir + "/cymobase_alignment_" + @file_id + ".fasta"

		# write data again to file, leave line breaks after 80 chars out (lucullus will be much faster this way)
		fh_dest = File.new(file_dest, 'w')
		File.open(file_scr, 'r').each_line do |line|
			line.chomp!
			str = ""
			if line[0] == ">" then
				str = "\n" << line << "\n"
			else
				str = line
			end
			fh_dest.write(str)
		end
		fh_dest.close
		render :show_alignment
	end

	# Predict CUG usage for uploaded data
	# 1) extract reference proteins
	# 2) gene prediction foreach reference protein:
	# 2.1) BLAST
	# 2.2) AUGUSTUS
	# 3) Compare with reference data
	# 3.1) Compare with reference alignment
	# 3.2) Compare with reference genes
	# render partial predict_genes
	# accessible params in view: @predicted_prots [Hash] Prediction data for each reference protein
	# accessible params in view: @stats [Hash] Statistics over prediction data
	# accessible params in view: @fata_error [Array] Fatal errors leading to program abort
	# accessible params in view: @minor_error [Array] Errors during gene prediction of single proteins, no program abort
	def predict_genes
		# add options to session
		session[:align] = { algo: params[:algo], config: params[:config] }
		session[:augustus] = { species: params[:species] }

		ref_data, @fatal_error = load_ref_data

		@predicted_prots = {} # containing final results ...
		@minor_error = []
		if @fatal_error.empty? then
			# sucessfully loaded reference data file
			# setup blast database
			genome_db = session[:file][:path].gsub(".fasta", "_db")
			error = create_blast_db(genome_db)
			if ! error.blank? then
				# database setup failed; no cug-usage can be predicted
				@fatal_error << error
				return
			end
	
			# all other tasks need to be done for each ref protein
			# peach for paralell execution, maximum 10 threads
			ref_data.peach(10) do |prot, all_prot_data|
				@predicted_prots[prot] = {} # ... final results for each protein!
				prot_basename = prot.gsub(" ", "-").downcase
				file_basename = session[:file][:path].gsub("query.fasta", prot_basename)

				### 1) extract reference proteins if they do not already exist
				file_refseq = file_basename.gsub(/\d+_/, "") + "-refseq.fasta"
				# ref_prot_species, ref_prot_stat, ref_prot_seq, ref_prot_key = "", "", "", ""
				ref_prot_seq, ref_prot_key = "", ""

				if File.exists?(file_refseq) then
					ref_prot_key, ref_prot_seq = fasta2str(File.read(file_refseq))
					ref_prot_key = ref_prot_key[0]
					ref_prot_seq = ref_prot_seq[0]
				else
					# find best reference sequence and fill variables with these data
					sc_genes = all_prot_data["genes"].keys.select {|key| key =~ /Sc_/ && all_prot_data["genes"][key]["completeness"] =~ /[complete|partial]/}
					if sc_genes.any? then
						# Sc_b and complete gene structure
						ref_prot_key = sc_genes.find {|key| all_prot_data["genes"][key]["completeness"] == "complete" &&  all_prot_data["genes"][key] =~ /Sc_b/}
						# ... or ... other Sc and complete gene structure
						ref_prot_key = sc_genes.find {|key| all_prot_data["genes"][key]["completeness"] == "complete"} if ref_prot_key.blank?
						# ... or ... Sc_b and partial
						ref_prot_key = sc_genes.find {|key| key =~ /Sc_b/} if ref_prot_key.blank?
					end
					# ... or ... use first gene with complete gene structure
					ref_prot_key = all_prot_data["genes"].keys.find {|key| all_prot_data["genes"][key]["completeness"] =~ /complete/} if ref_prot_key.blank?
					# ... or ... or just any gene at all
					ref_prot_key = all_prot_data["genes"].keys.first if ref_prot_key.blank?


					# not necessary if file already exists ...
					ref_prot_seq = all_prot_data["genes"][ref_prot_key]["gene"].match(/prot_seq: (.*?)\n/)[1]
					# save reference protein in fasta format to file
					File.open(file_refseq, 'w'){|f| f.write(str2fasta(ref_prot_key, ref_prot_seq))}
				end
					
#				currently, no more information about the reference protein is needed
#				ref_prot_species = all_prot_data["genes"][ref_prot_key]["species"]
#				ref_prot_stat = all_prot_data["genes"][ref_prot_key]["completeness"]
#				ref_prot_geneseq = extract_gene_seq(all_prot_data["genes"][ref_prot_key]["gene"])

				### 2) gene prediction foreach reference protein:

				### 2.1) BLAST
				file_blast = file_basename + ".blast"
						# -p [PROGRAM] protein query against nt-db -d [DATABASE] -i [FILE] -m8 [OUTPUT FORMAT] -F [FILTERING] -s [SMITH-WATERMAN ALIGNMENT] T 
				stdin, stdout_err, wait_thr = Open3.popen2e(BLASTALL, "-p", "tblastn", "-d", genome_db, "-i", file_refseq, "-m8", "-F", "m S", "-s", "T")
				stdin.close
				output = stdout_err.read
				stdout_err.close
				# save blast hits in file for predict_more method
				File.open(file_blast, 'w') {|f| f.write(output)}

				if ! wait_thr.value.success? || output.include?("ERROR") || output.blank? then
					write_minor_error(prot, "BLASTALL")
					@predicted_prots[prot][:message] = "Sorry, an error occured"
					# delete file with blast hits, its contains only error messages
					File.delete(file_blast)
					next
				end
				@predicted_prots[prot][:n_hits] = output.lines.to_a.size # number of all blast hits
				@predicted_prots[prot][:hit_shown] = 1 # in this method, only blast best hit is analyzed

				### 2.2) AUGUSTUS
				is_success, pred_seq, pred_dnaseq, err = perform_gene_pred(output)
				if ! is_success || ! err.blank? then
					# an error occured
					write_minor_error(prot, err)
					@predicted_prots[prot][:message] = "Sorry, an error occured"
					next
				end
		
				# next comes alignment, so the alignment of all reference sequences is needed
				if session[:algo] == "mafft" then
					# check if reference data are already aligned
					file_refall = BASE_PATH + prot_basename + ".fasta"
					if ! FileTest.file?(file_refall) then
						# no, so calculate alignment
						is_success = prepare_ref_data(prot_data["alignment"], file_refall)
						if ! is_success then
							write_minor_error(prot, "MAFFT")
							@predicted_prots[prot][:message] = "Sorry, an error occured"
							next
						end
					end
				end
				is_success, aligned_fasta, err = align_pred_ref(file_basename, pred_seq)
				if ! is_success || ! err.blank? then
					# an error occured
					write_minor_error(prot, err)
					@predicted_prots[prot][:message] = "Sorry, an error occured"
					next
				end

				### 3) compare with reference data
				is_success, results, message = compare_pred_gene(ref_data, ref_prot_key, prot, aligned_fasta, pred_dnaseq)
				if ! is_success then
					@predicted_prots[prot][:message] = message
				else
					@predicted_prots[prot][:message] = ""
				end
				@predicted_prots[prot].merge!(results)

			end # ref_data.peach do |prot, des|
		end # if @fatal_error.empty?
		@stats = calc_stats(@predicted_prots)

		render :predict_genes, formats: [:js]
	end

	# redo gene prediction for next 10 blast hits (gene prediction for best blast hit is already done) for one specific protein
	# same workflow as predict_genes
	# renders predict_more
	# accessible params in view: @predicted_prots [Hash] Prediction data for reference protein @prot
	# accessible params in view: @stats [Hash] Up-to-date statistics over prediction data
	# accessible params in view: @prot [String] Protein for which gene prediction was restarted
	# accessible params in view: @fata_error [Array] Fatal errors leading to program abort
	# accessible params in view: @minor_error [Array] Errors during gene prediction of single proteins, no program abort
	def predict_more
		@prot = params[:prot]
		hit_start = params[:hit].to_i + 1 # use next hit for prediction

		# prepare stuff
		@predicted_prots = {}
		@minor_error = []
		ref_data, @fatal_error = load_ref_data

		if @fatal_error.empty? then
			prot_basename = @prot.gsub(" ", "-").downcase
			file_basename = session[:file][:path].gsub("query.fasta", prot_basename)
			file_refseq = file_basename.gsub(/\d+_/, "") + "-refseq.fasta"
			file_blast = file_basename + ".blast"
			headers, seqs = fasta2str(File.read(file_refseq))
			ref_prot_key = headers[0]
			blast_all_hits = File.read(file_blast)
			if blast_all_hits.blank? then
				write_minor_error(@prot, "No more BLAST hits.")
				@predicted_prots[key][:message] = "Sorry, an error occured"
				return
			end
			n_hits = blast_all_hits.lines.to_a.size
			hit_stop = hit_start == 2 ? [10, n_hits].min : [hit_start+9, n_hits].min # go in steps of 10, or to last hit

			(hit_start..hit_stop).peach(10) do |n_hit|
				key = @prot + "_" + n_hit.to_s
				@predicted_prots[key] = { n_hits: n_hits, hit_shown: n_hit}

				is_success, pred_seq, pred_dnaseq, err = perform_gene_pred(blast_all_hits, n_hit)
				if ! is_success || ! err.blank? then
					write_minor_error(@prot, err)
					@predicted_prots[key][:message] = "Sorry, an error occured"
					next
				end

				# no need to check again if alignment of all reference data exist
				is_success, aligned_fasta, err = align_pred_ref(file_basename, n_hit, pred_seq)
				if ! is_success || ! err.blank? then
					# an error occured
					write_minor_error(@prot, err)
					@predicted_prots[key][:message] = "Sorry, an error occured"
					next
				end

				is_success, results, message = compare_pred_gene(ref_data, ref_prot_key, @prot, aligned_fasta, pred_dnaseq)
				if ! is_success then
					@predicted_prots[key][:message] = message
				else
					@predicted_prots[key][:message] = ""
				end
				@predicted_prots[key].merge!(results)
			end
		end # if @fatal_error.empty?
		@stats = calc_stats(@predicted_prots, prev_stats=true) # stats needs to be combined with old data!
		render :predict_more, formats: [:js]
	end

	# prepare a new session 
	def prepare_new_session
		reset_session
	    session[:file] = {} 		# uploaded genome described by keys: :name, :id, :path
	    session[:align] = {}		# alignment options for seqan: pair_align or DIALIGN
	    session[:augustus] = {}		# options for augustus gene prediction
	end

	# delete old data in /tmp/cug/ and /tmp/cymobase_alignment_cug*
	# data from old predictions and show_alignment requests
	def delete_old_data(days = 1)
		# delete all files older than one day except file "alignment_gene_structure.json"
		# system("find", BASE_PATH, "-type", "f", "-not", "-name", REF_DATA, "-mtime", "+#{days}", "-delete")
		begin
			Dir.glob(BASE_PATH + "*").each do |file|
				if file == BASE_PATH + REF_DATA then
					next
				elsif File.mtime(file) <= days.day.ago
					FileUtils.rm(file) # FileUtils.rm(file, :noop => true, :verbose => true)
				end
			end
			# delete alignment file for show-alignment function
			FileUtils.rm Dir.glob(Dir::tmpdir + '/cymobase_alignment_cug*')
		rescue
		end
	end

	# load reference data from file alignment_gene_structure.json
	# @return [FalseClass] If an error occured
	# @return [TrueClass] If no error occured -> Array will be empty
	# @return [Array] Errors occured during file load
	def load_ref_data
		path = BASE_PATH + REF_DATA
		if ! FileTest.file?(path) then
			errors = ["Fatal error.", "Cannot load reference data. Please contact us!"]
			return false, errors
		end
		data = JSON.load(File.read(path))

		# separate actin, myosin and kinesin by class
		myo_data = separate_by_class(data, "Myosin heavy chain")
		actin_data = separate_by_class(data, "Actin related protein")
		kinesin_data = separate_by_class(data, "Kinesin")
		# delete unseparated data 
		data.delete("Myosin heavy chain")
		data.delete("Actin related protein")
		data.delete("Kinesin")
		# add separated (do this after deletion of unseparated !!! as old key is used as "default class")
		data.merge!(myo_data)
		data.merge!(actin_data)
		data.merge!(kinesin_data)

		return data, []
	end

	# separate reference data for action, myosin and kinesin by class
	# @param data [Hash] reference data
	# @param prot [String] key in data - hash, those data to separate
	def separate_by_class(data, prot)
		new_data = {}
		genes = data[prot]["genes"].keys
		old_alignment = data[prot]["alignment"]
		headers, seqs = fasta2str(old_alignment)

		genes.each do |name|
			if prot.include?("Myosin") then
				name =~ /[a-zA-Z_]+(Myo)?([0-9]+|Mhc)/
				klass = $2 if $2
			elsif prot.include?("Actin") 
				name =~ /[a-zA-Z_]+(Arp)([0-9]+)/
				klass = $2 if $2
			else
				name =~ /[a-zA-Z_]+(Kinesin)([0-9]+)/
				klass = $2 if $2
			end
			
			if klass then
				# could extract class name $2
				key = prot + " Class " + $2.to_s # convert to key name
			else
				# could not extract a class name, simply use protein name as class name
				key = prot
			end

			if ! new_data.has_key?(key) then
				# initialize hash for this class
				new_data[key] = {}
				new_data[key]["genes"] = {}
				new_data[key]["alignment"] = ""
			end

			# add gene structure and aligned sequence to new data hash
			new_data[key]["genes"][name] = data[prot]["genes"][name]

			if headers.include?(name) then
				ind = headers.index(name)
				new_data[key]["alignment"] += str2fasta(name, seqs[ind], true) << "\n" # true: no line break after 80 chars
				# no else needed: even if no sequence exists for this key, the other functions can handle missing data
			end
		end
		return new_data
	end

	# Prepare alignment of reference data using MAFFT
	# calls prepare_ref_data in a subprocess and detaches them
	def bg_dataprocessing
		ref_data, error = load_ref_data
		if error.empty? then
			ref_data.each do |prot, prot_data|
				data_file = BASE_PATH + prot.gsub(" ", "-").downcase + ".fasta"
				if ! File.exists?(data_file) then
					job = fork { prepare_ref_data(prot_data["alignment"], data_file) }
					Process.detach(job)
				end
			end
		end
	end

	# prepare alignment of reference data
	# necessary if ref_data[prot]["alignment"] is not aligned, e.g. sequences have not same length
	# @param data [Array] ref_data[prot]["alignment"]
	# @param file [String] path for MAFFT - alignment file
	# @return [Boolean] indicates if an error occured
	def prepare_ref_data(data, file)
		file_in = file.gsub(".", "_in.")
		File.open(file_in, 'w'){|f| f.write(data)}
		stdin, stdout_err, wait_thr = Open3.popen2e(MAFFT, "--auto", "--amino", "--anysymbol", "--quiet", file_in)
		if ! wait_thr.value.success? then
			return false
		else
			stdin.close
			output = stdout_err.read
			stdout_err.close
			File.open(file, 'w'){|f| f.write(output)}
			FileUtils.rm(file_in)
		end	
		return true
	end

	def extract_ref_data_alignment(data, nogaps=true)
		alignment = ""
		if nogaps then
			data["genes"].each do |k, v|
				seq = v["gene"].match(/prot_seq: (.*?)\n/)[1]
				if seq.nil? then 
					next
				end
				alignment << str2fasta(k, seq, true) # true: no linebreak after 80 chars
				alignment << "\n"
			end
			return alignment
		else
			return data["alignment"]
		end
	end

	# parse predicted data to get statistics, store them in a file
	# and combine stored ones with the ones from additionally predicted if neccessary
	# @param data [Hash] Prediction data 
	# @param combine [Boolean] Optional; Indicates if existing stats needs to be combined with additional predicted ones
	# @return [Hash] Up-to-date statistics
	def calc_stats(data, combine=false)
		stats = Hash.new(0)
		file = BASE_PATH + session[:file][:id] + ".stat"

		if combine then
			# for this dataset, statistics already exist
			old_stats = File.read(file)
			stats = Hash[old_stats.scan(/([\w_]+):(\d+)/)] # parse data and put them right into the hash
			stats.each{ |key,val| stats[key] = val.to_i } # convert all numbers (string representation) to integers
		end

		data.each do |prot, val|
			if ! val[:ref_chem].nil? && val[:ref_chem].any? then
				# number of positions with preferential chemical properties
				stats["sig_pos"] += val[:ref_chem].keys.size
				# kind of preferential chemical properity leads to predited codon usage
				all = val[:ref_chem].collect{|_,v| v[:transl]}
				stats["ser"] += all.count("S")
				stats["leu"] += all.count("L")
			end
			if ! val[:ref_ctg].nil? && val[:ref_ctg].any?
				# number of proteins with CTGs at same position
				stats["ref_ctg"] += 1 if val[:ref_ctg].values.max >= 0.05
			end
		end
		stats["n_prots"] += data.keys.length

		# save updated data statistics to file 
		str = stats.map {|obj| obj.join(":") }.join("\n")
		File.open(file, 'w'){|f| f.write(str)}
		return stats
	end

	# compare uploadedfile-size with maximal size (50MB)
	# @param size [Fixnum] File size 
	# @return [Array] Error message if size exceeds MAX_SIZE
	def check_filesize(size)
		errors = []
		if size > MAX_SIZE then
			errors << "Fatal error."
			errors << "File must be less than 50 MB."
			errors << "Please contact us to upload larger files."
		end
		return errors
	end

	# Read fasta and test if it is fasta-formatted genome file and contains an CTG codon
	# @param file [String] File handle
	# @return [Array] Error messages if invalid file
	def check_fasta(file)
		errors = []
		is_fasta = true
		contains_ctg = false
		is_first_line = true
		last_nucleotides = "" # need to store 2 last nts of each line to check for ctg!
		# read file line by line, since it might be large
		IO.foreach(file.path) do |line|
			line.chomp!
			if is_first_line then
				# expect very first line to be fasta header
				if ! line.starts_with?('>') then
					is_fasta = false
					# test failed, stop reading file!
					break
				end
				is_first_line = false
			else
				# it's not the very first line, might be header, comment or sequence
				if line.starts_with?('>') then
					# a new header: clear last_nucleotides-buffer
					last_nucleotides = ""
					next
				elsif line.starts_with?(';')
					# comment line, simply skip
					next
				elsif line =~ /[^ACGTURYMKSWBDHVN]/i
					# a sequence line, but containing charaters not part of iupac nucleotide definition
					# test failed, stop reading file!
					is_fasta = false
					break
				else
					# a sequence line, valid content
					# ctg might be splitted into 2 lines
					contains_ctg = true if (last_nucleotides + line).upcase.include?("CTG")
					if line.length >=2 then
						last_nucleotides = line[-2,2]
					else
						last_nucleotides = line
					end
				end
			end
		end

		# file read, check what happend
		if ! is_fasta then
			errors << "Fatal error."
			errors << "Invalid input."
			errors << "Expected fasta formatted genome file."
		end
		if  is_fasta && (! contains_ctg) then
			errors << "Fatal error."
			errors << "Genome sequence contains no \'CTG\'."
			errors << "Cannot predict codon usage."
		end
		return errors
	end

	# creating blast db from uploaded genome file with formatdb
	# for the blast search performed later
	# @param genome_db [String] Path for genome database
	# @return [String] Error if formatdb failed
	def create_blast_db(genome_db)
		errors = []
			# -i [INPUT: GENOME FILE] -p [PROTEIN] FALSE -n [DB NAME] -o [CREATE INDEX OVER SEQID] TRUE
		is_success = system(FORMATDB, "-i", session[:file][:path], "-p", "F", "-n", genome_db, "-o", "T")
		# no output needed, so system is sufficient
		if ! is_success then
			errors << "Fatal error."
			errors << "Cannot create BLAST database from genome data."
		end 
		return errors
	end

	# build up @minor_error [Array] Error messages
	# non-fatal errors, affecting only a single protein
	# on fab8, verbose error messages will be printed, outside non-verbose ones
	# @param prot [String] Protein for which gene prediction failed
	# @param message [String] More detailed information about the error
	def write_minor_error (prot, message)
		if @minor_error.empty? then
			@minor_error << "Cannot prediction genes in"
		end
		str = Verbose_error ? "#{prot}: #{message}" : prot
		@minor_error << "- #{str}"
	end

	# parse blast output 
	# @param hits [String] Complete blast output (Hit list)
	# @param nr [Fixnum] Number of blast hit of interest (human counting -> starts with 1)
	# @return [TrueClass] If the requested hit exists
	# @return [String] Contig containing hit
	# @return [Fixnum] Sequence start
	# @return [Fixnum] Sequence stop
	# @return [Fixnum] Strand (1 for plusstrand, 2 for minusstrand)
	# @return [FalseClass] AND NO OTHER VALUES If the requested hit does not exist
	def parse_blast_hits(hits, nr)
		hit = hits.lines.to_a[nr-1] # hit number count starts with 1 (human counting), ruby array count starts with 0
		if hit.nil? then
			# requested hit execceds number of hits
			return false 
		end
		# hit exits
		fields = hit.split("\t") 
		seq_id = fields[1] # sequence contig
		start = fields[8].to_i # sequence start
		stop = fields[9].to_i # sequence stop
		strand = 1 # set strand to plus for fastacmd
		if stop < start then
			# change start/ stop if prediction on minus strand 
			strand = 2 # set strand to minus for fastacmd
			tmp = start
			start = stop
			stop = tmp
		end
		return true, seq_id, start, stop, strand
	end

	# perform gene prediction with augustus: run augustus and parse output
	# input: 
	# @param seq_file [String] Path to file containing coding region 
	# @return [Boolean] True if no error occured, false otherwise
	# @return [String] Predicted protein sequence
	# @return [String] Predicted coding sequence (DNA)
	def run_augustus(file)
			# --species [REFERENCE SPEC] QUERY --genemodel=exactlyone [predict exactly one complete gene] --codingseq=on [output also coding sequence]
			# redirection: add stderr to stdout (screen)
		stdin, stdout_err, wait_thr = Open3.popen2e(AUGUSTUS, "--AUGUSTUS_CONFIG_PATH=/usr/local/bin/augustus/config/", 
			"--species=#{session[:augustus][:species]}", "--genemodel=exactlyone","--codingseq=on", file)
		stdin.close
		output = stdout_err.read
		stdout_err.close
		if (! wait_thr.value.success?) || output.include?("ERROR") || (! output.include?("coding sequence")) then
			return false
		end

		# no error, parse augustus
		pred_seq, pred_dnaseq = parse_augustus(output)
		if pred_seq.nil? || pred_dnaseq.nil? then
			return false
		end
		return true, pred_seq, pred_dnaseq, ""
	end


	# performing gene prediction for given protein (step 2.2)
	# calling augustus based on best blast hit 
	# @param blast_hits [String] Complete blast output (Hit list)
	# @param nr [Fixnum] Number of blast hit of interest (human counting -> starts with 1)
	# @return [Boolean] Indicating if an error occured
	# @return [String] Predicted protein sequence
	# @return [String] Predicted coding sequence (DNA)
	# @return [String] Error message, only set if an error occured
	def perform_gene_pred(blast_hits, blast_hit_nr=1)
		# parse blast output
		is_success, seq_id, start, stop, strand = parse_blast_hits(blast_hits, blast_hit_nr) 
		if ! is_success then
			# number exceeds blast hits!
			return false, "", "", "Number of BLAST hits exceeded"
		end 

		# 2.2) gene prediction
		# 2.2.1) get matching search sequence 
		file = BASE_PATH + session[:file][:id] + "_" + rand(1000000).to_s + ".fasta"
		genome_db = session[:file][:path].gsub(".fasta", "_db")
		# enlarge range of sequence to extract: add 1000 nucleodides to start/stop
		if start-1000 < 0 then
			# check is start is valid
			start = 0
		else
			start = start-1000
		end
		stop = stop+1000
				# -d [DB] -p [PROTEIN] false -s [SEARCH STRING] -L [START,STOP]
				# add 1000 nucleotides to start/stop 
				# redirection: stderr to stdout (at this moment: screen), stdout to file file
		stdin, stdout_err, wait_thr = Open3.popen2e(FASTACMD, "-d", genome_db, "-p", "F", "-s", seq_id, 
			"-L", "#{start},#{stop}", "-S", strand.to_s, "-o", file)
		stdin.close
		output = stdout_err.read
		stdout_err.close
		if ! wait_thr.value.success? then
			return false, "", "", "BLAST-FASTACMD"
		end

		# 2.2.2) AUGUSTUS
		is_success, pred_seq, pred_dnaseq = run_augustus(file)
		if ! is_success then
			# gene prediction failed
			return false, "", "", "AUGUSTUS"
		end

		return true, pred_seq, pred_dnaseq, ""
	end

	# align reference protein with the predicted protein
	# @param file_base [String] basename of files (containing protein name)
#	# @param hit [Fixnum] Number of blast hit of interest (human counting -> starts with 1)
	# @param pred_seq [String] Predicted protein sequence
	# @return [Boolean] True if no error occured
	# @return [Array] Aligned sequences
	# @return [String] Error indication which method failed only set if an error occured
	def align_pred_ref(file_base, hit=1, pred_seq)
		file_unaligned = file_base + "-" + rand(1000000).to_s + ".fasta" # input file: unaligned seqs
		file_aligned = file_base + "-" + hit.to_s + "-aligned.fasta" # output file: aligned seqs
		fasta = str2fasta("Prediction", pred_seq, true) # true: do not split after 80 chars

		if session[:align][:algo] == "mafft" then
			# file_unaligned contains only predicted sequence
			File.open(file_unaligned, 'w') {|f| f.write(fasta)}
			file_refall = file_base.gsub(/\d+\_/, "") + ".fasta"
			if File.zero?(file_refall) then
				# something went wrong during ref_data generation (cymo-api)
				return false, "", "reference data"
			end
		else
			# file_unaligned contains both reference and predicted sequence
			file_refseq = file_base.gsub(/\d+\_/, "") + "-refseq.fasta"
			FileUtils.cp(file_refseq, file_unaligned)
			File.open(file_unaligned, 'a') {|f| f.write("\n" + fasta)}
		end

		if session[:align][:algo] == "dialign" then
			# use system, as dialign has no easy-parsable error messages (like containing "ERROR")
			# check exit status for success 
			is_success = system(DIALIGN2, "-fa", "-fn", file_aligned, file_unaligned)
			if ! is_success then
				return false, "", "Dialign2"
			end
			begin
				File.rename(file_aligned + ".fa", file_aligned) # dialign automatically adds ".fa" to basic output filename
			rescue
				return false, "", "Dialign2"
			end
		# elsif session[:align][:algo] == "clustalw"
		# 	# use system as no easy parsealbe output
		# 	is_success = system(CLUSTAW, "-infile=", file_unaligned, "-quiet", "-quicktree", "-outfile=", file_aligned, "-output=fasta")
		# 	if ! is_success then
		# 		return false, "", "ClustalW"
		# 	end			
		elsif session[:align][:algo] == "mafft"
				# --addfragments: add non-full length sequence to existing MSA
				# --amino: input is protein
				# --anysymbol: replace unusal symbols by "X"
				# --quiet: output only alignment
			stdin, stdout_err, wait_thr = Open3.popen2e(MAFFT, "--addfragments", file_unaligned, "--amino", "--anysymbol", "--quiet", file_refall)
			stdin.close
			output = stdout_err.read
			stdout_err.close
			if ! wait_thr.value.success? then
				return false, "", "MAFFT"
			end
		else
				# -c tttt: initialize first row & column with zeros, search last row & column for maximum
				# 	=> end gap free alignment
				# --matrx matrix file
				# implicit options: gap open penalty: -11, gap extension penalty: -1
				# 					protein sequences
				# don't use system, as pair_align is quite verbose and has an easy parsable success-output
			stdin, stdout_err, wait_thr = Open3.popen2e(PAIR_ALIGN, "--seq", file_unaligned, "--matrix", Rails.root.join('lib', 'blosum62').to_s, 
				"--config", session[:align][:config], "--method", session[:align][:algo], "--outfile", file_aligned)
			stdin.close
			output = stdout_err.read
			stdout_err.close
			if ! wait_thr.value.success? || ! output.include?("Alignment score") then
				return false, "", "Seqan::pair_align"
			end
		end

		# parse multiple sequence alignment to get aligned reference and aligned predicted sequence
		# header, seqs (last seq is predicted seq)
		if session[:align][:algo] == "mafft" then
			File.open(file_aligned, 'w') {|f| f.write(output)}
			fasta = output
		else
			fasta = File.read(file_aligned)
		end
		# ensure that the methods returns an aligned fasta-formatted string!
		# or ends with "false" otherwise
		if fasta.blank? then
			return false, "", "Internal error"
		end
		return true, fasta
	end

	# preforming cug-usage prediction 
	# step 3.1) Compare with reference alignment
	# step 3.1.1) map predicted gene onto reference protein
	# step 3.2) Compare with reference genes
	# @param ref_data [Hash] Cymobase reference data
	# @param ref_prot_key [String] Reference protein identifier; key in ref_data
	# @param prot [String] Protein name; key in ref_data
	# @param aligned_fasta [String] Fasta-formatted aligned sequences (all reference seqs in case of MAFFT), only 2 otherwise
	# @param pred_dnaseq [String] Predicted coding sequence (DNA)
	# @return [TrueClass] If no error occured
	# @return [Hash] Prediction data for protein prot
	# @return [FalseClass] If an error occured
	# @return [Hash] Prediction data for protein prot ONLY IN CASE OF AN ERROR
	# @return [String] Error message only if an error occured
	def compare_pred_gene(ref_data, ref_prot_key, prot, aligned_fasta, pred_dnaseq)
		# results = {has_ctg: "", pred_prot: "", ctg_pos: [], ctg_transl: [], aa_comp: [], ctg_ref: []}
		results = {ctg_pos: [], ref_chem: {}, ref_ctg: {}}
		headers, seqs = fasta2str(aligned_fasta)
		pred_seq_aligned = seqs[-1]
		codons = pred_dnaseq.scan(/.{1,3}/)
		ctg_pos = codons.each_with_index.map{ |ele, ind| (ele == "CTG") ? ind : nil }.compact
		results[:ctg_pos] = ctg_pos

		# has predicted protein CTG codon(s)?
		if ctg_pos.empty? then
			# nope, no CTG => nothing to predict
			# results[:has_ctg] = false
			return false, results, "No CTG in predicted gene"
		else
			# yes, compare CTG positions with reference data
			# results[:has_ctg] = true
			results[:pred_prot] = pred_seq_aligned.gsub("-", "")
		end

		# Compare with reference alignment
		# => map CTG positions in predicted protein onto reference alignment columns
		# => compare with chemical properties of reference alignment columns
		# => Compare with reference genes: count usage of CTG for each mapped L and S 

		if session[:align][:algo] == "mafft" then
			# aligned_fasta is a MSA containing all reference data and the predicted sequence
			ref_alignment = Hash[headers.zip(seqs)]
			ref_alignment.delete("Prediction") # todo really delete it? but only without predictiion its an real ref_alignment
			is_mafft = true
		else
			is_mafft = false
			begin
			ref_alignment = Hash[*ref_data[prot]["alignment"].split("\n")]
			ref_alignment.keys.each {|k| ref_alignment[ k.sub(">", "") ] = ref_alignment.delete(k)}
			rescue => e
				puts "---"
				puts "ERROR (#{prot}):"
				puts e
				puts "---"
				return false, results, "Internal error"
			end
		end

		# CTG positions in ref_alignment and amino acids at respective positions
		ctg_pos_mapped, ref_alignment_cols = map_ctg_pos(ref_alignment, ctg_pos, aligned_fasta, is_mafft)

		if ctg_pos_mapped.any? then
			# mapping was possible, compare!

			# 1) preference for hydrophob/ polar residues?
			ref_alignment_cols.each_with_index do |col, ind|
				seq_num = ref_alignment.keys.size # total number of sequences

				if col.count("-") < seq_num / 2 then
					# less than half of all seqs have an gap at the CTG position 
					# => check the preference makes sence
					aa_freq, aa_num = word_frequency(col)
# todo delete all freqs < 5 % ???
					is_significant, prob_transl = predict_translation(aa_freq)
				else
					# half or more seqs have an gap => its not significant
					is_significant = false
				end

				if is_significant then
					results[:ref_chem][ctg_pos[ind]] = {transl: prob_transl, aa_comp: aa_freq, aa_num: aa_num, seq_num: seq_num}
				end
			end

			# 2) CTGs in reference data at CTG positions in predicted sequence?
			pct_ser, pct_leu = ref_ctg_usage(ref_alignment, ctg_pos_mapped, ref_data[prot]["genes"])
			if pct_ser > 0 || pct_leu > 0 then
				results[:ref_ctg] = {ser: pct_ser, leu: pct_leu}
			end

		else 
			# mapping was not possible
			return false, results, "Cannot match CTG position"
		end
		return true, results
	end

# alter kram
	def tmp_oldcompare_pred_gene
		ctg_pos.each do |pos|
			# map CTG positions, only neccessary if not all seqs were aligned
			if session[:align][:algo] == "mafft" then
				algnmnt_col_pos = pos 
				algnmnt_col_aa = seqs.collect{|seq| seq[pos]} 
				cymo_algnmnt = Hash[headers.zip(seqs)]
			else
				ref_seq_aligned = seqs[0]
				begin
					cymo_algnmnt = Hash[*ref_data[prot]["alignment"].split("\n")]
				rescue => e
					puts "---"
					puts "ERROR (#{prot}):"
					puts e
					puts "---"
					return false, results, "Internal error"
				end
				algnmnt_col_aa, algnmnt_col_pos = parse_alignment_by_ctgpos(cymo_algnmnt, pos, pred_seq_aligned, ref_seq_aligned, ref_prot_key)
			end
			if !(algnmnt_col_aa.nil? && algnmnt_col_pos.nil?) then
				# mapping was possible, compare with chemical properties
				aa_freq, aa_num = word_frequency(algnmnt_col_aa)
				seq_num = cymo_algnmnt.keys.length
				if aa_num <= seq_num/2 then
					# half or more of all seqs have an gap at position pos -> its not significant
					is_significant = false
				end
				is_significant, prob_transl = predict_translation(aa_freq)
				# compare with reference genes
				pct_ctg_ser, pct_ctg_leu = ref_ctg_usage(cymo_algnmnt, algnmnt_col_pos, ref_data[prot]["genes"])
				if is_significant then
					results[:ctg_pos] << pos
					results[:ctg_transl] << prob_transl
					results[:aa_comp] << [aa_freq, aa_num, seq_num]
					results[:ctg_ref] << {ser: pct_ctg_ser, leu: pct_ctg_leu}
				else
					return false, results, "No significant position"
				end # if is_significant
			else
				return false, results, "Cannot match CTG position"
			end # if !(algnmnt_col_aa.nil? && algnmnt_col_pos.nil?)
		end # ctg_pos.each do |pos|
		return true, results
	end


	# convert header and sequence into fasta-format
	# @param header [String] fasta header
	# @param seq [String] sequence
	# @param no_split [Boolean] include line break each 80 chars?
	# @return [String] fasta formatted header and sequence
	def str2fasta(header, seq, no_split=false)
		fasta = header.include?(">") ? header << "\n" : ">" << header << "\n"
		if no_split then
			fasta += seq
		else
			fasta += seq.scan(/.{1,80}/).join("\n")
		end
		return fasta
	end

	# extract headers and sequences from a fasta-formatted string
	# can deal with multiple sequence fasta
	# @param fasta [String] fasta-formatted sequence (MSA also possible)
	# @return [Array] All fasta header found
	# @return [Array] All fasta sequences found
	def fasta2str(fasta)
		headers = []
		seqs = []
		# get every record: everything between two ">"
		fasta.scan(/^>([^>]*)/m).flatten.each do |rec|
			rec.chomp!
			nl = rec.index("\n") # first match: separating header from seq
			header = rec[0..nl-1]
			seq = rec[nl+1..-1]
			seq.gsub!(/\n/,'') # removing all remaining \n from seq 
			headers.push(header)
			seqs.push(seq)
		end
		return headers, seqs
	end

	# extracting coding sequence out of gene entry 
	# @param gene [String] Gene entry from reference data
	# @return [String] Coding sequence of this gene
	def extract_gene_seq(gene)
		a = gene.scan(/\sseq:\s?(\w+)\n/).flatten
		# quick & dirty version: every second entry belongs to an exon
		a.values_at(*a.each_index.select(&:even?)).join("").upcase
	end

	# extraction protein and coding sequence from augustus output
	# @param output [String] Output from augustus
	# @return [String] Predicted protein sequence
	# @return [String] Predicted coding sequence (DNA)
	def parse_augustus(output)
		pred_prot = output.match(/protein sequence = \[(.*)\]/m)[1]
		pred_dna = output.match(/coding sequence = \[(.*?)\]/m)[1]
		# prepare output
		pred_prot.gsub!("\n# ", "").upcase!
		pred_dna.gsub!("\n# ", "").upcase!
		return pred_prot, pred_dna
	end

	# calculate position of residue in unaligned sequence based on the position in aligned sequence
	# @param apos [Fixnum] Position in aligned sequence
	# @param aseq [String] Aligned sequence
	# @return [Fixnum] Position in unaligned sequence 
	def alignment_pos2sequence_pos(apos, aseq)
		aseq[0..apos].gsub("-", "").length - 1
	end

	# calculate position of residue in aligned sequence based on the position in unaligned sequence
	# @param spos [Fixnum] Position in unaligned sequence
	# @param aseq [String] Aligned sequence
	# @return [Fixnum] Position in aligned sequence 
	def sequence_pos2alignment_pos(spos, aseq)
		pats = []
		aseq.gsub("-", "")[0..spos].split("").each {|chr| pats << ("-*" + chr)}
		pat = Regexp.new(pats.join)
		pat.match(aseq)[0].length - 1
	end


	# map CTG positions in predicted sequence (unaligned) onto the cymobase - alignment
	# mapping is done via the alignment between reference and predicted sequence
	# handels alignment method (in case of mafft only cymobase alignment at CTG pos needs to be collected!)
	# @param ref_alignment [Hash] cymobase alignment
	# @param ctg_pos [Array] CTG-Positions in unaligned predicted sequence
	# @param aligned_fasta [String] Alignment between reference and predicted seq
	# @param is_mafft [Boolean] Indicates if alignment method was mafft - influences aligned_fasta etc
	# @return [Array] CTG-Positions in the cymobase - alignment
	# @return [Array of Arrays] Cymobase - Alignment columns at CTG positions
	def map_ctg_pos(ref_alignment, ctg_pos, aligned_fasta, is_mafft)
		ref_pos = []
		ref_cols = []

		if is_mafft then
			# reference alignment covers all seqs from cymobase!
			ref_pos = ctg_pos # same position as in reference alignment
			ctg_pos.each do |pos|
				ref_cols << ref_alignment.collect {|cymo_header, cymo_prot| cymo_prot[pos]}
			end
		else
			header, seqs = fasta2str(aligned_fasta)
			ref_key = header[0]
			ref_seq = seqs[0] # aligned with pred_seq, but not with ref_alignment
			pred_seq = seqs[1]
			ctg_pos.each do |pos|
				# map CTG position in unaligned onto aligned predicted seq
				pred_apos = sequence_pos2alignment_pos(pos, pred_seq)

				# matching only possible if
					# 1) no gap in reference sequence at this position 
				if (ref_seq[pred_apos] == "-" || 
					# 2) this position is not aligned (only relevant in case of DIALIGN)
					ref_seq[pred_apos].match(/\p{Lower}/) || pred_seq[pred_apos].match(/\p{Lower}/)) then
					next
				end

				# continue matching: map aligend predicted seq onto reference alignment
				ref_spos = alignment_pos2sequence_pos(pred_apos, ref_seq)
				if ref_spos.nil? || ref_alignment[ref_key].nil? then
					next
				end
				ref_cymopos = sequence_pos2alignment_pos(ref_spos, ref_alignment[ref_key])

				# collect cymo column at CTG position
				col = ref_alignment.collect {|cymo_header, cymo_prot| cymo_prot[ref_cymopos]}
				
				# store results
				ref_pos << ref_cymopos
				ref_cols << col
			end
		end
		return ref_pos, ref_cols
	end

	# counting word frequencies; #doesnot count "-"
	# @param arr [Arr] Containing words (e.g. amino acids)
	# @return [Hash] Keys are the words, values their relative frequency
	# @return [Fixnum] Sum of all occurrences (e.g. total number of amino acids)
	def word_frequency(arr)
		res = Hash.new(0)
		arr.each { |a| res[a] += 1 }
		res.delete(nil) # delete "nil" key (some reference sequence were shorter than ctg positon)
		# res.delete("-") # delete count for gaps in alignment!
		# normalize word frequency
		sum = res.inject(0) {|s, (_,val)| s + val}.to_f
		sum = sum - res["-"] # substract gaps, to count only amino acids
		norm_res = res.each {|k,v| res[k] = v/sum} # normalize
		return norm_res, sum.to_i
	end

	# predict the most probable translation based on amino acids statistics
	# @param aa_freq [Hash] Amino acid frequencies
	# @return [Boolean] Are stats discriminative?
	# @return [String] Most probable CUG-translation ("S" or "L")
	def predict_translation(aa_freq)
		aa_freq.delete("-")
		if aa_freq.empty? then
			# only gaps, and the "-" key deleted
			return false
		end
		num_aas = aa_freq.inject(0) {|sum, (_,val)| sum + val} # total number of amino acids
		pol_aas = POLAR_AAS.collect{|aa| aa_freq[aa]}.sum # polar amino acids
		hyd_aas = HYDORPHOBIC_AAS.collect{|aa| aa_freq[aa]}.sum # hydrophobic amino acids
		pct_pol = pol_aas/num_aas.to_f
		pct_hyd = hyd_aas/num_aas.to_f
		# requirements for a discriminative position:
			# 1) occurence in more than half of sequences
			# 2) the other usage should occure in less than quater or sequences
		is_discrim = ([pct_hyd, pct_pol].max > 0.5 && [pct_hyd, pct_pol].min < 0.25) ? true : false 

		# default translation = ""
			# discriminative AND more hydrophobic aas: "L"
			# discriminative AND more polar aas: "S"
		transl = ""
		if is_discrim then
			transl = (pol_aas > hyd_aas) ? "S" : "L"
		end
		return is_discrim, transl
	end



	# check the CTG usage in reference sequences
	# @param alignment [Hash] Cymobase alignment
	# @param pos_list [Array] CTG positions in cymo alignment
	# @param genes [Hash] Gene structrues for cymo alignment
	# @return [Fixnum] Percentage of Serine (at CTG position in predicted seq) encoded by CTG (rounded)
	# @return [Fixnum] Percentage of Leucine (at CTG position in predicted seq) encoded by CTG (rounded)
	def ref_ctg_usage(alignment, pos_list, genes)
		# counts for Serine, Leucine encoded by CTG and in general (regardless codon)
		ser, leu, total = 0, 0, 0
		genes.each do |key, gene|
			dna_seq = extract_gene_seq(gene["gene"])
			ref_codons = dna_seq.scan(/.{1,3}/)
			pos_list.each do |pos|
				ref_seq = alignment[key]
				ref_aa = ref_seq[pos]
				spos = alignment_pos2sequence_pos(pos, ref_seq) # position in unaligned sequence
				if ref_aa == "S" && ref_codons[spos] == "CTG" then
					ser += 1
					total += 1
				elsif ref_aa == "L" && ref_codons[spos] == "CTG"
					leu += 1
					total += 1
				elsif ref_aa == "S" || ref_aa == "L"
					# any other codon 
					total += 1
				end # ser/leu count
			end # pos_list.each
		end # genes.each

		# convert to relative frequencies
		pct_ser = ser.to_f / total
		pct_ser = 0.to_f if pct_ser.nan?
		pct_leu = leu.to_f / total
		pct_leu = 0.to_f if pct_leu.nan?
		return pct_ser.round(2), pct_leu.round(2)
	end


	# # check for a given position in the cymo alignment, how many "S" and how many "L" are encoded by an CTG
	# # @param algnmnt [Hash] Cymobase alignment
	# # @param apos [Fixnum] CTG position in cymo alignment
	# # @param genes [Hash] Genes from reference data
	# # @return [Fixnum] Percentage of "Ser" encoded by CTGs (rounded to 2 decimal places)
	# # @return [Fixnum] Percentage of "Leu" encoded by CTGs (rounded to 2 decimal places)
	# def ref_ctg_usage_alt(algnmnt, apos, genes)
	# 	pct_ctg_S, pct_ctg_L = "", "" 
	# 	# Serine
	# 	ref_codons = []
	# 	algnmnt.collect{|k,v| k if v[apos]=="S"}.compact.each do |ref_key|
	# 		dna_seq = extract_gene_seq(genes[ref_key.gsub(">", "")]["gene"]) # access via "gene-name"
	# 		codons_ref = dna_seq.scan(/.{1,3}/)
	# 		spos = alignment_pos2sequence_pos(apos, algnmnt[ref_key]) # access via ">gene-name"
	# 		ref_codons << codons_ref[spos]
	# 	end
	# 	pct_ctg_S = ref_codons.count("CTG").to_f/ref_codons.length 
	# 	pct_ctg_S = 0.0 if pct_ctg_S.nan?

	# 	# same for Leucine
	# 	ref_codons = []
	# 	algnmnt.collect{|k,v| k if v[apos]=="L"}.compact.each do |ref_key|
	# 		dna_seq = extract_gene_seq(genes[ref_key.gsub(">", "")]["gene"]) # access via "gene-name"
	# 		codons_ref = dna_seq.scan(/.{1,3}/)
	# 		spos = alignment_pos2sequence_pos(apos, algnmnt[ref_key]) # access via ">gene-name"
	# 		ref_codons << codons_ref[spos]
	# 	end
	# 	pct_ctg_L = ref_codons.count("CTG").to_f/ref_codons.length
	# 	pct_ctg_L = 0.0 if pct_ctg_L.nan?
	# 	return pct_ctg_S.round(2), pct_ctg_L.round(2)
	# end

	# this method is not needed any more; 
	# # match CTG position to reference alignment
	# # @param algnmnt [Hash] Cymobase alignment
	# # @param pred_spos [Fixnum] Position in unaligned predicted sequence
	# # @param pred_aseq [String] Aligned predicted sequence
	# # @param ref_aseq [String] Aligned reference sequence
	# # @param ref_key [String] Reference protein; key in algnmnt
	# # @return [Array] Residues at matched alignment column
	# # @return [Fixnum] CTG position in reference sequence aligned to cymobase alignment
	# # @return [NilClass] If matching not possible (gap in cymobase alignment) AND
	# # @return [NilClass] If matching not possible (gap in cymobase alignment)
	# def parse_alignment_by_ctgpos(algnmnt, pred_spos, pred_aseq, ref_aseq, ref_key)
	# 	# 1) ctg pos in unaligned seq -> pos in aligned predicted sequence (aligned with reference seq)
	# 	pred_apos = sequence_pos2alignment_pos(pred_spos, pred_aseq) # = pos_ref_seq_aligned

	# 	# matching possible at all?
	# 	# NO: gap in reference sequence at ctg pos in predicted sequence
	# 	# or, in case of dialign: important positions are not aligned at all (= in lowercase letters)
	# 	if (ref_aseq[pred_apos] == "-" || 
	# 		ref_aseq[pred_apos].match(/\p{Lower}/) || pred_aseq[pred_spos].match(/\p{Lower}/)) then
	# 		return nil, nil
	# 	end

	# 	# YES
	# 	# 2) pos aligned -> cymo alignment
	# 	begin
	# 		ref_spos = alignment_pos2sequence_pos(pred_apos, ref_aseq)
	# 		ref_cymopos = sequence_pos2alignment_pos(ref_spos, algnmnt[">"<<ref_key])
	# 		# 3) return column
	# 		col = algnmnt.collect {|cymo_prot, cymo_seq| cymo_seq[ref_cymopos]}
	# 	rescue => e
	# 		puts "---"
	# 		puts "ERROR (#{ref_key}):"
	# 		puts e
	# 		puts "---"
	# 		return nil, nil
	# 	end
	# 	return col, ref_cymopos
	# end
end
