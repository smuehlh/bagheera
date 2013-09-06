require 'peach'
require 'Tools'
require 'open3'

class PredictionsController < ApplicationController

	MAX_SIZE = 15728640 # 15 MB
	HYDORPHOBIC_AAS = ["V", "I", "L", "M", "F"]
	POLAR_AAS = ["S", "T", "A", "C"]

	# Render start page for prediction
	def search
		delete_old_data # delete old uploaded genome files
		prepare_new_session # a fresh session
		# File.open(Progress_file, 'w') { |f| f.write(write_status(1, 1, false)) }# setup status iframe
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
		filesize = params[:uploaded_file].size
		@fatal_error |= check_filesize(filesize)
		# check content, only if file is not too big
		check_cug = (filesize < 5120 ) ? true : false  # 5 KB
		@fatal_error |= check_fasta(params[:uploaded_file], check_cug) if @fatal_error.blank?
		if @fatal_error.blank? then
			# rename and save file
			file_id = rand(1000000000).to_s
			file_name = File.basename(params[:uploaded_file].original_filename)
			file_path = BASE_PATH + file_id + "/query.fasta"
			begin
				FileUtils.mkdir(BASE_PATH + file_id, :mode => 0775)
				FileUtils.mv(params[:uploaded_file].path, file_path)
			rescue
				@fatal_error << "Cannot upload file."
			end
			session[:file] = { id: file_id, name: file_name, path: file_path }
		end
		render :upload_file_ajax, formats: [:js]
	end

	# Handel loading the example genome:
	# store it in /tmp/cug and add its id(new file name), original filename and path to the session
	# render partial showing the uploaded file
	# accessible params: @fatal_error [Array] Errors occured during file load
	def load_example
		example_file = "Candida_albicans_WO_1.fasta"
		example_path = Rails.root + "spec/fixtures/files/" + example_file
		@fatal_error = []
		if File.exist?(example_path) then
			# do not check content, but check file size, just to be sure
			@fatal_error |= check_filesize(File.size(example_path))
			if @fatal_error.blank? then
				file_id = rand(1000000000).to_s
				file_path = BASE_PATH + file_id + "/query.fasta"
				begin
					FileUtils.mkdir(BASE_PATH + file_id, :mode => 0775)
					FileUtils.cp(example_path, file_path)
				rescue
					@fatal_error << "Cannot upload file."
				end
				session[:file] = { id: file_id, name: example_file, path: file_path }
			end
		else
			@fatal_error << "Cannot upload file."
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
		file_scr = BASE_PATH + session[:file][:id] + "/" + prot + "-" + hit.to_s + "-aligned.fasta"
		file_dest = Dir::tmpdir + "/cymobase_alignment_" + @file_id + ".fasta"

		# write data again to file, leave line breaks after 80 chars out (lucullus will be much faster this way)
		headers, seqs = fasta2str(File.read(file_scr))
		headers.map!{ |e| ">" << e } # add ">" again to make it a valid header
		pred_seq_ind = headers.index(">Prediction")

		if ! pred_seq_ind.nil? then
			# found an element called "Prediction", remove it from its position and add it to the front
			if pred_seq_ind != 0 then
				# ... of course, unless it is already at the first pos
				pred_header = headers.delete_at(pred_seq_ind)
				pred_seq = seqs.delete_at(pred_seq_ind)
				# add them again at first position
				headers.unshift(pred_header)
				seqs.unshift(pred_seq)
			end
		end

		fh_dest = File.new(file_dest, 'w')
		fh_dest.puts( headers.zip(seqs).flatten.join("\n") )
	# QUICK-FIX to Lucullus-Bug
	# # write predicted seq at the end of alignment, to make sure no seq is missing (Prediction is already at the top)
	# fh_dest.puts( str2fasta("Prediction", seqs[headers.index(">Prediction")], true) )
		fh_dest.close
		render :show_alignment
	end

	# def read_status
	# 	render :text => File.open(Progress_file, "r"){ |f| f.read }, :layout => false
	# end

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
		use_low_comp_filter = params.has_key?(:blast) ? true : false
		session[:align] = { algo: params[:algo], config: params[:config] }
		session[:augustus] = { species: params[:species] }
		session[:blast] = { use_low_comp_filter: use_low_comp_filter}

		if session[:file][:name] == "Candida_albicans_WO_1.fasta" then
			ref_data, @fatal_error = load_ref_data(true) # true: use ref_data without candida albicans wo 1
			ref_files_path = BASE_PATH + PATH_REF_WO_EXAMPLE
		else
			ref_data, @fatal_error = load_ref_data
			ref_files_path = BASE_PATH
		end 

		@predicted_prots = {} # containing final results ...
		@minor_error = [] # containing external program errors
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

			# # setup up progress file for (nearly) synchronous status messaging
			# n_prot = ref_data.keys.size
			# this_prot = 1

			# File.open(Progress_file, 'w') { |f| f.write(write_status(this_prot, n_prot, false)) }
			

			# all other tasks need to be done for each ref protein
			# peach for paralell execution, maximum 10 threads
			ref_data.peach(10) do |prot, all_prot_data|

				@predicted_prots[prot] = {} # ... final results for each protein!
				prot_basename = prot.gsub(" ", "-").downcase
				file_basename = session[:file][:path].gsub("query.fasta", prot_basename)

				### 1) get "best" reference protein according to selected augustus model
				file_refseq = file_basename + "-refseq.fasta"
				ref_prot_seq, ref_prot_key = "", ""

				matched_sp_abbr = ""
				case session[:augustus][:species]
				when 'candida_albicans'
					if session[:file][:name] == "Candida_albicans_WO_1.fasta" then
						# example genome Ca_b was loaded, and deleted from reference data
						matched_sp_abbr = "Ca_a"
					else
						matched_sp_abbr = "Ca_b"
					end
				when 'candida_guilliermondii'
					matched_sp_abbr = "Mrg"
				when 'candida_tropicalis'
					matched_sp_abbr = "Ct_a"
				when 'debaryomyces_hansenii'
					matched_sp_abbr = "Deh"
				when 'eremothecium_gossypii'
					matched_sp_abbr = "Erg"
				when 'kluyveromyces_lactis'
					matched_sp_abbr = "Kl"
				when 'lodderomyces_elongisporus'
					matched_sp_abbr = "Loe"
				when 'pichia_stipitis'
					matched_sp_abbr = "Shs"
				when 'saccharomyces_cerevisiae_S288C'
					matched_sp_abbr = "Sc_c"
				when 'saccharomyces_cerevisiae_rm11-1a_1'
					matched_sp_abbr = "Sc_b"
				when 'yarrowia_lipolytica'
					matched_sp_abbr = "Yl"
				end

				genes = all_prot_data["genes"].keys.select {|key| key =~ /#{matched_sp_abbr}/ && all_prot_data["genes"][key]["completeness"] =~ /[complete|partial]/}

				if genes.size == 1 then
					# no choice at all
					ref_prot_key = genes[0]
				elsif genes.any?
					# correct species and complete gene structure
					ref_prot_key = genes.find {|key| all_prot_data["genes"][key]["completeness"] == "complete" &&  all_prot_data["genes"][key] =~ /#{matched_sp_abbr}/}
					# ... or ... partial
					ref_prot_key = genes.find {|key| key =~ /#{matched_sp_abbr}/} if ref_prot_key.blank?
				end
				# ... or ... use first gene with complete gene structure
				ref_prot_key = all_prot_data["genes"].keys.find {|key| all_prot_data["genes"][key]["completeness"] =~ /complete/} if ref_prot_key.blank?
				# ... or ... or just any gene at all
				ref_prot_key = all_prot_data["genes"].keys.first if ref_prot_key.blank?
				ref_prot_seq = extract_prot_seq(all_prot_data["genes"][ref_prot_key]["gene"])

				# save reference protein in fasta format to file
				header = all_prot_data["genes"][ref_prot_key]["species"] + "@" + ref_prot_key
				File.open(file_refseq, 'w'){|f| f.write(str2fasta(header, ref_prot_seq))}

#				currently, some more information about the reference protein is needed
				@predicted_prots[prot][:ref_species] = all_prot_data["genes"][ref_prot_key]["species"]
				@predicted_prots[prot][:ref_prot] = ref_prot_seq
#				ref_prot_stat = all_prot_data["genes"][ref_prot_key]["completeness"]
#				ref_prot_geneseq = extract_gene_seq(all_prot_data["genes"][ref_prot_key]["gene"])

				@predicted_prots[prot][:message] = []

				### 2) gene prediction foreach reference protein:

				### 2.1) BLAST
				file_blast = file_basename + ".blast"
				filter_option = session[:blast][:use_low_comp_filter] ? "m S" : "F"

						# -p [PROGRAM] protein query against nt-db -d [DATABASE] -i [FILE] -m8 [OUTPUT FORMAT] -F [FILTERING] -s [SMITH-WATERMAN ALIGNMENT] T 
						# -F "m S" mask protein low complexity filter for lookup table only
				stdin, stdout_err, wait_thr = Open3.popen2e(BLASTALL, "-p", "tblastn", "-d", genome_db, "-i", file_refseq, "-m8", "-F", filter_option, "-s", "T")
				stdin.close
				output = stdout_err.read
				stdout_err.close
				# save blast hits in file for predict_more method
				File.open(file_blast, 'w') {|f| f.write(output)}

				if ! wait_thr.value.success? || output.include?("ERROR") || output.blank? then
					write_minor_error(prot, "BLAST failed.")
					@predicted_prots[prot][:message] |= ["BLAST failed."] #["Sorry, an error occured"]
					# delete file with blast hits, its contains only error messages
					File.delete(file_blast)
					next
				end

				@predicted_prots[prot][:n_hits] = output.lines.to_a.size # number of all blast hits
				@predicted_prots[prot][:hit_shown] = 1 # in this method, only blast best hit is analyzed
				

				### 2.2) AUGUSTUS
				is_success, pred_seq, pred_dnaseq, err, no_prot_profile = perform_gene_pred(output, prot_basename)

				if ! is_success then
					# some big error occured
					write_minor_error(prot, err)
					@predicted_prots[prot][:message] |= [err] #["Sorry, an error occured"]
					next
				end
				if ! no_prot_profile.blank? then
					# special case: augustus could not use profiles, but still could predict something
					# keep error message but continue with CTG usage prediction

					@predicted_prots[prot][:message] |= [no_prot_profile]
				end

				# next comes alignment, so the alignment of all reference sequences is needed
				file_refall = ref_files_path + prot_basename + ".fasta"

				if session[:algo] == "mafft" then
					# check if reference data are already aligned
					if ! FileTest.file?(file_refall) then
						# no, so calculate alignment
						puts "Need to calculate reference alignment with MAFFT first!"
						is_success = prepare_mafft_ref_alignment(prot_data["alignment"], file_refall)
						if ! is_success then
							write_minor_error(prot, "Aligning sequences with MAFFT failed.")
							@predicted_prots[prot][:message] |= ["Aligning sequences failed."] #["Sorry, an error occured"]
							next
						end
					end
				end

				is_success, aligned_fasta, err = align_pred_ref(file_basename, file_refall, pred_seq)
				if ! is_success || ! err.blank? then
					# an error occured
					write_minor_error(prot, err)
					@predicted_prots[prot][:message] |= [err] #["Sorry, an error occured"]
					next
				end

				### 3) compare with reference data
				is_success, results, message = compare_pred_gene(ref_data, ref_prot_key, prot, aligned_fasta, pred_dnaseq)
				if ! is_success then
					@predicted_prots[prot][:message] |= [message]
				end
				@predicted_prots[prot].merge!(results)

				# # write nearly synchronos status file
				# File.open(Progress_file, 'w') { |f| f.sync = true; f.write(write_status(this_prot,n_prot,false)) }
				# this_prot += 1 # protein counter

			end # ref_data.peach do |prot, des|

		end # if @fatal_error.empty?
		# write nearly synchronos status file
		# File.open(Progress_file, "w") { |f| f.sync = true; f.write(write_status(n_prot,n_prot,true)) }
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
		@minor_error = [] # external program errors
		if session[:file][:name] == "Candida_albicans_WO_1.fasta" then
			ref_data, @fatal_error = load_ref_data(true)
			ref_files_path = BASE_PATH + PATH_REF_WO_EXAMPLE
		else
			ref_data, @fatal_error = load_ref_data
			ref_files_path = BASE_PATH
		end 

		if @fatal_error.empty? then
			prot_basename = @prot.gsub(" ", "-").downcase
			file_basename = session[:file][:path].gsub("query.fasta", prot_basename)

			file_refseq = file_basename + "-refseq.fasta"
			file_blast = file_basename + ".blast"
			file_refall = ref_files_path + prot_basename + ".fasta"

			headers, seqs = fasta2str(File.read(file_refseq))
			ref_species, ref_prot_key = headers[0].split("@")

			blast_all_hits = File.read(file_blast)
			if blast_all_hits.blank? then
				write_minor_error(@prot, "No more BLAST hits.")
				return
			end
			n_hits = blast_all_hits.lines.to_a.size
			hit_stop = hit_start == 2 ? [10, n_hits].min : [hit_start+9, n_hits].min # go in steps of 10, or to last hit

			(hit_start..hit_stop).peach(10) do |n_hit|
				key = @prot + "_" + n_hit.to_s
				@predicted_prots[key] = { n_hits: n_hits, hit_shown: n_hit}
				@predicted_prots[key][:ref_species] = ref_species
				@predicted_prots[key][:ref_prot] = seqs[0]
				@predicted_prots[key][:message] = []

				is_success, pred_seq, pred_dnaseq, err, no_prot_profile = perform_gene_pred(blast_all_hits, n_hit, prot_basename)
				if ! is_success && ! err.blank? then
					# some big error occured
					write_minor_error(key, err)
					@predicted_prots[key][:message] |= ["Sorry, an error occured"]
					next
				end
				if ! no_prot_profile then
					# special case: augustus could not use profiles, but still could predict something
					# keep error message but continue with CTG usage prediction
					@predicted_prots[key][:message] |= [no_prot_profile]
				end

				# no need to check again if alignment of all reference data exist
				is_success, aligned_fasta, err = align_pred_ref(file_basename, file_refall, n_hit, pred_seq)
				if ! is_success || ! err.blank? then
					# an error occured
					write_minor_error(key, err)
					@predicted_prots[key][:message] |= ["Sorry, an error occured"]
					next
				end

				is_success, results, message = compare_pred_gene(ref_data, ref_prot_key, @prot, aligned_fasta, pred_dnaseq)
				if ! is_success then
					@predicted_prots[key][:message] |= [message]
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
		# delete all files older than one day except file "alignment_gene_structure.json" and logfile
		# system("find", BASE_PATH, "-type", "f", "-not", "-name", REF_DATA, "-mtime", "+#{days}", "-delete")
		begin
			Dir.glob(BASE_PATH + "**/") do |dir|
				if dir == BASE_PATH then
					next
				elsif dir == BASE_PATH + PATH_REF_WO_EXAMPLE 
					next
				elsif File.mtime(dir) <= days.day.ago
					FileUtils.rm_rf(dir, :secure=>true)
				end
			end
			# delete alignment file for show-alignment function
			FileUtils.rm Dir.glob(Dir::tmpdir + '/cymobase_alignment_cug*')
		rescue
		end
	end

	def write_status(done, total, final=false)
		html = ["<html>", "<head>"]
		html += ['<meta http-equiv="refresh" content="5;url=./status.html">'] if ! final 
		html += ["</head>", "<body>", "Processed #{done}/#{total}", "</body>", "</html>"]
		return html.join(" ")
	end

	# load reference data from file alignment_gene_structure.json
	# @param use_ref_data_without_ca_b [Boolean] 
		# set true if the reference data without the species loaded via 'load example' should be used
		# default: false
	# @return [Hash] reference data
	# @return [Array] Errors occured during file load
	def load_ref_data(use_ref_data_without_ca_b=false)
		if use_ref_data_without_ca_b then
			path = BASE_PATH + PATH_REF_WO_EXAMPLE + REF_DATA
		else
			path = BASE_PATH + REF_DATA
		end

		if ! FileTest.file?(path) then
			errors = ["Fatal error.", "Cannot load reference data. Please contact us!"]
			return false, errors
		end
		data = JSON.load(File.read(path))

		return data, []
	end

	def prepare_mafft_ref_alignment(data, file)
		File.open(file, 'w') {|f| f.write(data)}
		stdin, stdout_err, wait_thr = Open3.popen2e("/usr/local/bin/mafft", "--auto", "--amino", "--anysymbol", "--quiet", file)
		if ! wait_thr.value.success? then
			return false
		end
		stdin.close
		output = stdout_err.read
		stdout_err.close
		File.open(file, 'w') {|f| f.write(output)}
		return true
	end

	# parse predicted data to get statistics, store them in a file
	# and combine stored ones with the ones from additionally predicted if neccessary
	# @param data [Hash] Prediction data 
	# @param combine [Boolean] Optional; Indicates if existing stats needs to be combined with additional predicted ones
	# @return [Hash] Up-to-date statistics
	def calc_stats(data, combine=false)
		stats = Hash.new(0)
		file = BASE_PATH + session[:file][:id] + "/stat"

		if combine then
			# for this dataset, statistics already exist
			old_stats = File.read(file)
			stats = Hash[old_stats.scan(/([\w_]+):(\d+)/)] # parse data and put them right into the hash
			stats.each{ |key,val| stats[key] = val.to_i } # convert all numbers (string representation) to integers
		end

		data.each do |prot, val|
			if ! val[:ref_chem].nil? && val[:ref_chem].any? then
				# number of CTG positions mapped to a CTG position in reference data -> is already conserved!
				stats["conserved_pos"] += val[:ref_ctg].keys.size

				# number of proteins which contain conserved CTG positions
				stats["prot_conserved_pos"] += 1 if val[:ref_ctg].keys.size > 0

				# number of CTGs in predicted prot
				stats["pred_ctg"] += val[:ctg_pos].size

				# kind of preferential chemical properity leads to predited codon usage
				arr = view_context.suggested_transl(val[:ref_chem], val[:ref_ctg], true) # true: return array containing alt_usage and std_usage
				stats["ser"] += arr[0] # alt_usage: predicted transl is serine
				stats["leu"] += arr[1] # std_usage: predicted transl is leucine
			end
			# if ! val[:ref_ctg].nil? && val[:ref_ctg].any?
			# 	stats["ref_ctg"] += 1 if val[:ref_ctg].values.collect{|k, _| k[:ctg_num].to_f/val[:ref_seq_num] }.max >= 0.05
			# end
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
	# @param check_cug [Boolean] Determine to check for presence of CUG codons too (only for really small files)
	# @return [Array] Error messages if invalid file
	def check_fasta(file, check_cug)
		errors = []
		is_fasta = true
		contains_ctg = false
		is_first_line = true
		last_nucleotides = "" # need to store 2 last nts of each line to check for ctg!
		counter = 0 # if a simple check for fasta needs, read only first 100 lines
		# read file line by line, since it might be large
		IO.foreach(file.path) do |line|
			line.chomp!
			counter += 1
			break if counter > 100 && ! check_cug # simple check is done
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
					# a new header
					# the quick check for fasta format is done
					break unless check_cug
					# continue only if check for cug - codon also needed: clear last_nucleotides-buffer
					last_nucleotides = ""
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
		@minor_error << prot + ": " + message
	end

	# parse blast output 
	# combines hit of interest with the following (!) blast hits, if they overlap
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

		# test if another hit on same contig exits, which overlaps

		# all indices of hits about same sequence id
		indices = hits.lines.to_a.each_with_index.map{ |ele, ind| (ele.include?(seq_id)) ? ind : nil }.compact
		# iterate over hits listed after (!) the hit of interest; 
		# this saves us from merging hits a second hit (only important for predict_more)!
		indices -= ( 0..(nr-1) ).to_a

		indices.each do |ind|
			this_hit = hits.lines.to_a[ind]
			this_parts = this_hit.split("\t")
			this_start = this_parts[8].to_i
			this_stop = this_parts[9].to_i
			this_strand = 1
			if this_stop < this_start then
				# change start/ stop if prediction on minus strand 
				this_strand = 2 # set strand to minus for fastacmd
				tmp = this_start
				this_start = this_stop
				this_stop = tmp
			end

			# by extending the hit of interest, it will overlap with this hit, so merge them
			if this_strand == strand && stop+500 >= this_start && stop <= this_start then
				# overlap with stop
				stop = this_stop
			end
			if this_strand == strand && start-500 <= this_stop && start >= this_stop then
				# overlap with start
				start = this_start
			end

		end
		return true, seq_id, start, stop, strand
	end


	# perform gene prediction with augustus: run augustus and parse output
	# input: 
	# @param seq_file [String] Path to file containing coding region 
	# @return [Boolean] True if no error occured, false otherwise
	# @return [String] Predicted protein sequence
	# @return [String] Predicted coding sequence (DNA)
	# @return [String] Error message (covers the availability of the protein profile)
	def run_augustus(file, prot_basename)


			# --species [REFERENCE SPEC] QUERY --genemodel=exactlyone [predict exactly one complete gene] --codingseq=on [output also coding sequence]
			# --strand=forward [use only forward strand for gene prediction; this works as fastacmd always translates into forward strand]
			# --proteinprofile profilefile [use protein profile about alignment as foundation for gene prediciton]
			# redirection: add stderr to stdout (screen)

		prfl_file = BASE_PATH + prot_basename + ".prfl"
		if File.exists?(prfl_file) && ! File.zero?(prfl_file) then
			# use protein profile for prediction
			stdin, stdout_err, wait_thr = Open3.popen2e(AUGUSTUS, "--AUGUSTUS_CONFIG_PATH=/usr/local/bin/augustus/config/", 
				"--species=#{session[:augustus][:species]}", "--genemodel=exactlyone", "--strand=forward", "--codingseq=on", 
				"--proteinprofile=#{prfl_file}", file)
		else
			# profile file does not exist, dont use it
			stdin, stdout_err, wait_thr = Open3.popen2e(AUGUSTUS, "--AUGUSTUS_CONFIG_PATH=/usr/local/bin/augustus/config/", 
				"--species=#{session[:augustus][:species]}", "--genemodel=exactlyone", "--strand=forward", "--codingseq=on", file)
			
			err = "No protein profile available for gene prediction. "
		end	

		stdin.close
		output = stdout_err.read
		stdout_err.close

		if wait_thr.value.termsig == 6 && err.nil? then
			# memory error while using protein profile, try again without protein profile
			stdin, stdout_err, wait_thr = Open3.popen2e(AUGUSTUS, "--AUGUSTUS_CONFIG_PATH=/usr/local/bin/augustus/config/", 
				"--species=#{session[:augustus][:species]}", "--genemodel=exactlyone", "--strand=forward", "--codingseq=on", file)
			
			err = "Cannot use protein profile for gene prediction."
			stdin.close
			output = stdout_err.read
			stdout_err.close
			if (! wait_thr.value.success?) || output.include?("ERROR") || ! output.include?("coding sequence") then
				return false, "", "", err
			end
		elsif (! wait_thr.value.success?) || output.include?("ERROR") || ! output.include?("coding sequence")
			return false, "", "", err
		end

		# no error, parse augustus
		pred_seq, pred_dnaseq = parse_augustus(output)
		if pred_seq.nil? || pred_dnaseq.nil? then
			return false, "", "", err
		end
		return true, pred_seq, pred_dnaseq, err||=""
	end


	# performing gene prediction for given protein (step 2.2)
	# calling augustus based on best blast hit 
	# @param blast_hits [String] Complete blast output (Hit list)
	# @param nr [Fixnum] Number of blast hit of interest (human counting -> starts with 1)
	# @return [Boolean] Indicating if an error occured
	# @return [String] Predicted protein sequence
	# @return [String] Predicted coding sequence (DNA)
	# @return [String] Error message, only set if an error occured
	# @return [String] Error message, only for protein profile (only set if no prot profile used)
	def perform_gene_pred(blast_hits, blast_hit_nr=1, prot_basename)
	
		# parse blast output
		is_success, seq_id, start, stop, strand = parse_blast_hits(blast_hits, blast_hit_nr) 
		if ! is_success then
			# number exceeds blast hits!
			return false, "", "", "Number of BLAST hits exceeded."
		end 

		# 2.2) gene prediction
		# 2.2.1) get matching search sequence 
		file = BASE_PATH + session[:file][:id] + "/" + rand(1000000).to_s + ".fasta"
		genome_db = session[:file][:path].gsub(".fasta", "_db")
		# enlarge range of sequence to extract: add 500 nucleodides to start/stop
		if start-500 < 0 then
			# check if start is valid
			start = 0
		else
			start = start-500
		end
		stop = stop+500 # if stop is out of range, fastacmd will cause an error
		# 1) try fastacmd with this stop value
				# -d [DB] -p [PROTEIN] false -s [SEARCH STRING] -L [START,STOP]
				# add 1000 nucleotides to start/stop 
				# redirection: stderr to stdout (at this moment: screen), stdout to file file
		stdin, stdout_err, wait_thr = Open3.popen2e(FASTACMD, "-d", genome_db, "-p", "F", "-s", seq_id, 
			"-L", "#{start},#{stop}", "-S", strand.to_s, "-o", file)
		stdin.close
		output = stdout_err.read
		stdout_err.close

		is_stop_outofrange = output.match(/From location cannot be greater than/)
		# 2) stop was out of range, use 0 instead -> this will be evaluated as maximal value by fastacmd!
		if is_stop_outofrange then
			stop = 0 # evaluated as maximal value (by fastacmd)
			stdin, stdout_err, wait_thr = Open3.popen2e(FASTACMD, "-d", genome_db, "-p", "F", "-s", seq_id, 
				"-L", "#{start},#{stop}", "-S", strand.to_s, "-o", file)
			stdin.close
			output = stdout_err.read
			stdout_err.close
		end

		if ! wait_thr.value.success? then
			return false, "", "", "Cannot start gene prediction."
		end

		# 2.2.2) AUGUSTUS
		is_success, pred_seq, pred_dnaseq, has_prot_profile = run_augustus(file, prot_basename)
		if ! is_success then
			# gene prediction failed
			return false, "", "", "Gene prediction failed."
		end

		return true, pred_seq, pred_dnaseq, "", has_prot_profile
	end

	# align reference protein with the predicted protein
	# @param file_base [String] basename of files (containing protein name)
	# @param ref_algn_path [String] path to reference alignment of this protein
#	# @param hit [Fixnum] Number of blast hit of interest (human counting -> starts with 1)
	# @param pred_seq [String] Predicted protein sequence
	# @return [Boolean] True if no error occured
	# @return [Array] Aligned sequences
	# @return [String] Error indication which method failed only set if an error occured
	def align_pred_ref(file_base, ref_algn_path, hit=1, pred_seq)

		file_unaligned = file_base + "-" + rand(1000000).to_s + ".fasta" # input file: unaligned seqs
		file_aligned = file_base + "-" + hit.to_s + "-aligned.fasta" # output file: aligned seqs
		fasta = str2fasta("Prediction", pred_seq, true) # true: do not split after 80 chars

		if session[:align][:algo] == "mafft" then
			# file_unaligned contains only predicted sequence
			File.open(file_unaligned, 'w') {|f| f.write(fasta)}
			if File.zero?(ref_algn_path) then
				# something went wrong during ref_data generation (cymo-api)
				return false, "", "Cannot align predicted protein with reference data."
			end
		else
			# file_unaligned contains both reference and predicted sequence
			file_refseq = file_base + "-refseq.fasta"
			FileUtils.cp(file_refseq, file_unaligned)
			File.open(file_unaligned, 'a') {|f| f.write("\n" + fasta)}
		end

		if session[:align][:algo] == "dialign" then
			# use system, as dialign has no easy-parsable error messages (like containing "ERROR")
			# check exit status for success 
			is_success = system(DIALIGN2, "-fa", "-fn", file_aligned, file_unaligned)
			if ! is_success then
				return false, "", "Cannot align predicted protein with reference data."
			end
			begin
				File.rename(file_aligned + ".fa", file_aligned) # dialign automatically adds ".fa" to basic output filename
			rescue
				return false, "", "Cannot align predicted protein with reference data."
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
			stdin, stdout_err, wait_thr = Open3.popen2e(MAFFT, "--addfragments", file_unaligned, "--amino", "--anysymbol", "--quiet", ref_algn_path)
			stdin.close
			output = stdout_err.read
			stdout_err.close
			if ! wait_thr.value.success? then
				return false, "", "Cannot align predicted protein with reference data."
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
				return false, "", "Cannot align predicted protein with reference data."
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
			return false, "", "Internal error."
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
		results[:pred_prot] = pred_seq_aligned.gsub("-", "")
		# results[:ref_prot] = seqs[headers.index(ref_prot_key)].gsub("-", "")

		# has predicted protein CTG codon(s)?
		if ctg_pos.empty? then
			# nope, no CTG => nothing to predict
			# results[:has_ctg] = false
			return false, results, "No CTG in predicted gene"
		else
			# yes, compare CTG positions with reference data
			# results[:has_ctg] = true
		end

		# Compare with reference alignment
		# => map CTG positions in predicted protein onto reference alignment columns
		# => compare with chemical properties of reference alignment columns
		# => Compare with reference genes: count usage of CTG for each mapped L and S 

		if session[:align][:algo] == "mafft" then
			# aligned_fasta is a MSA containing all reference data and the predicted sequence
			ref_alignment = Hash[headers.zip(seqs)]
			ref_alignment.delete("Prediction") # otherwise prediced seq and its CTGs counts twice! (and its not a real reference alignment)
		else
			begin
			ref_alignment = Hash[*ref_data[prot]["alignment"].split("\n")]
			ref_alignment.keys.each {|k| ref_alignment[ k.sub(">", "") ] = ref_alignment.delete(k)}
			rescue => e
				puts "---"
				puts "ERROR (#{prot}):"
				puts e
				puts "---"
				return false, results, "Internal error."
			end
		end

		results[:ref_seq_num] = ref_alignment.keys.size # total number of reference sequences

		# CTG positions in ref_alignment and amino acids at respective positions
		ind = headers.find_index{|str| str.include?("@")}
		if ind then
			# we used an alignment method other than mafft, so parse the header to extract real key
			seq = seqs[ind]
		else
			seq = seqs[headers.index(ref_prot_key)]
		end
		ref_aligned_fasta = str2fasta(ref_prot_key, seq, true) # true: no split after 80 chars
		pred_aligned_fasta = str2fasta("Prediction", pred_seq_aligned, true) # true: no split after 80 chars

		ctg_pos_mapped, ref_alignment_cols, ref_codons = map_ctg_pos(ref_alignment, ctg_pos, ref_aligned_fasta, pred_aligned_fasta, ref_data[prot]["genes"])

		if ctg_pos_mapped.compact.any? then
			# mapping of at least one position was possible, compare!

			# 1) preference for hydrophob/ polar residues?
			ref_alignment_cols.each_with_index do |col, ind|
				# handle unmatched CTG positions
				next if ctg_pos_mapped[ind].nil?

				# if col.count("-") < seq_num / 2 then
					# less than half of all seqs have an gap at the CTG position 
					# => check the preference makes sence
					aa_freq, aa_num = word_frequency(col)
					is_significant, prob_transl = predict_translation(aa_freq)
				# else
					# half or more seqs have an gap => its not significant
					# is_significant = false
				# end

				results[:ref_chem][ctg_pos[ind]] = {aa_comp: aa_freq, aa_num: aa_num}
				if is_significant then
					results[:ref_chem][ctg_pos[ind]][:is_significant] = true
					results[:ref_chem][ctg_pos[ind]][:transl] = prob_transl
				else
					results[:ref_chem][ctg_pos[ind]][:is_significant] = false
				end
			end

			# 2) CTGs in reference data at CTG positions in predicted sequence?

			# this method handels unmatched CTG positions on its own
			results[:ref_ctg] = ref_ctg_usage(ref_codons, ref_alignment_cols, ctg_pos)

		else 
			# mapping was not possible
			return false, results, "Cannot match CTG position"
		end
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
			seq.gsub!(/\r?\n/,'') || seq # removing all remaining \n from seq 
			headers.push(header)
			seqs.push(seq)
		end
		return headers, seqs
	end

	# extracting coding sequence from gene entry 
	# @param gene [String] Gene entry from reference data
	# @return [String] Coding sequence of this gene
	def extract_gene_seq(gene)
		# # this does not work in case of exons located at different contigs
		# a = gene.scan(/\sseq:\s?(\w+)\n/).flatten
		# # quick & dirty version: every second entry belongs to an exon
		# a.values_at(*a.each_index.select(&:even?)).join("").upcase

		if ! (defined? ScipioResult.new_from_string) then
			require 'scipio_result'
		end	
		sr = ScipioResult.new_from_string(gene)
		sr.cdna.upcase
	end

	# extracting protein sequence from gene entry
	# @param gene [String] Gene entry from reference data
	# @return [String] Protein sequence of this gene
	def extract_prot_seq(gene)
		# gene.match(/prot_seq: (.*?)\n/)[1].upcase

		if ! (defined? ScipioResult.new_from_string) then
			require 'scipio_result'
		end	
		sr = ScipioResult.new_from_string(gene)
		sr.translation.upcase
	end

	# extract a certain codon from gene entry
	# @param gene [String] Gene entry from reference data
	# @param pos [String] Position of the codon of interest
	# @return [String] codon of interest
	def get_codon(gene, pos)
		dna_seq = extract_gene_seq(gene) # dna sequence
		codons = dna_seq.scan(/.{1,3}/)
		return codons[pos]
	end


	# extraction protein and coding sequence from augustus output
	# @param output [String] Output from augustus
	# @return [String] Predicted protein sequence
	# @return [String] Predicted coding sequence (DNA)
	def parse_augustus(output)
		pred_prot = output.match(/protein sequence = \[(.*)\]/m)[1]
		pred_dna = output.match(/coding sequence = \[(.*?)\]/m)[1]
		# prepare output
		# using the non-destructive variant is important, as gsub! returns nil if no substituion made
		prot = pred_prot.gsub("\n# ", "").upcase
		dna = pred_dna.gsub("\n# ", "").upcase
		return prot, dna
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
	# also collects codons used in reference seq
	# placeholder are added to output for each unmapped CTG position
	# @param ref_alignment [Hash] cymobase alignment
	# @param ctg_pos [Array] CTG-Positions in unaligned predicted sequence
	# @param ref_fasta [String] fasta formatted reference sequence
	# @param pred_fasta [String] fasta formatted predicted sequence
	# @return [Array] CTG-Positions in the cymobase - alignment
	# @return [Array of Arrays] Cymobase - Alignment columns at CTG positions
	# @return [Array of Arrays] Cymobase - Codons at CTG positions
	def map_ctg_pos(ref_alignment, ctg_pos, ref_fasta, pred_fasta, genes)

		ref_pos = []
		ref_cols = []
		ref_codons = []

		header, seq = fasta2str(ref_fasta)
		ref_key = header[0]
		ref_seq = seq[0] # aligned with pred_seq, but not with ref_alignment
		header, seq = fasta2str(pred_fasta)
		pred_seq = seq[0]

		ctg_pos.each do |pos|
			# map CTG position in unaligned onto aligned predicted seq
			pred_apos = sequence_pos2alignment_pos(pos, pred_seq)

			# matching only possible if
				# 1) no gap in reference sequence at this position 
			if (ref_seq[pred_apos] == "-" || 
				# 2) this position is not aligned (only relevant in case of DIALIGN)
				ref_seq[pred_apos].match(/\p{Lower}/) || pred_seq[pred_apos].match(/\p{Lower}/)) then
				# add placeholder to results
				ref_pos << nil
				ref_cols << []
				ref_codons << []
				next
			end

			# continue matching: map aligend predicted seq onto reference alignment
			ref_spos = alignment_pos2sequence_pos(pred_apos, ref_seq)
			if ref_spos.nil? || ref_alignment[ref_key].nil? then
				# add placeholder to results
				ref_pos << nil
				ref_cols << []
				ref_codons << []
				next
			end
			ref_cymopos = sequence_pos2alignment_pos(ref_spos, ref_alignment[ref_key])

			col = []
			codons = []
			ref_alignment.each do |cymo_header, cymo_prot|

				# collect cymo column at CTG position
				col << cymo_prot[ref_cymopos]

				# collect codons at CTG position
				if cymo_prot[ref_cymopos] == "-" then
					# add nil as placeholder if its a gap
					codons << nil
				else
					# actually collect codon
					spos = alignment_pos2sequence_pos(ref_cymopos, cymo_prot) # position in sequece without gaps
					codons << get_codon(genes[cymo_header]["gene"], spos)
				end
			end

			# store results
			ref_pos << ref_cymopos
			ref_cols << col
			ref_codons << codons.compact # if only nil- placeholders, it will be an empty array
		end
		return ref_pos, ref_cols, ref_codons
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
			# 2) the other usage should occure in less than half or sequences
		is_discrim = ([pct_hyd, pct_pol].max >= 0.5 && [pct_hyd, pct_pol].min < 0.45) ? true : false 

		# default translation = ""
			# discriminative AND more hydrophobic aas: "L"
			# discriminative AND more polar aas: "S"
		transl = ""
		if is_discrim then
			transl = (pol_aas > hyd_aas) ? "S" : "L"
		end
		return is_discrim, transl
	end


	# parse reference codons 
	# @param codons [Array of Arrays] Codons for each ctg positions
	# @param aas [Array of Arrays] Alignment columns at ctg positions
	# @param ctg_pos [Array] List of CTG positions
	# @return [Hash] keys: ctg pos, value: {:ctg_usage => {"S" => xx, "L", xx}, ctg_num => xx, :seq_num => xx}
	def ref_ctg_usage(codons, aas, ctg_pos)
		counts = {}
		codons.each_with_index do |codons_thiscol, ind|
			# check codons at this alignment column (identical with iteration over CTG position)
			indices = codons_thiscol.find_each_index("CTG")
			ctg_num = indices.size
			if ctg_num > 0 then
				# ctg_num: all seqs, also those with a gap at this pos
				counts[ctg_pos[ind]] = {ctg_usage: Hash.new(0), ctg_num: ctg_num} 
				counts[ctg_pos[ind]][:ctg_usage]["S"] = 0
				counts[ctg_pos[ind]][:ctg_usage]["L"] = 0
				aa_thiscol = aas[ind].reject{|ele| ele == "-"} # now amino acids match codons
				indices.each do |i| 
					if aa_thiscol[i] == "S" then
						counts[ctg_pos[ind]][:ctg_usage]["S"] += 1
					elsif aa_thiscol[i] == "L"
						counts[ctg_pos[ind]][:ctg_usage]["L"] += 1
					end # ser/leu count
				end # indices.each

				# normalize counts and save them to results
				pct_ser = counts[ctg_pos[ind]][:ctg_usage]["S"] / counts[ctg_pos[ind]][:ctg_num].to_f
				pct_leu = counts[ctg_pos[ind]][:ctg_usage]["L"] / counts[ctg_pos[ind]][:ctg_num].to_f
				counts[ctg_pos[ind]][:ctg_usage]["S"] = pct_ser
				counts[ctg_pos[ind]][:ctg_usage]["L"] = pct_leu

				# check if something went wrong
				if (pct_ser + pct_leu) == 0.0 then
					# dont delete key from counts to make the mistake accessible
					counts[ctg_pos[ind]][:transl] = "???"
					next
				end

				# find most probable translation
				is_discrim, counts[ctg_pos[ind]][:transl] = predict_translation(counts[ctg_pos[ind]][:ctg_usage])

			end # if ctg_num
		end # codons.each_with_index
		return counts
	end


	# def ref_ctg_usage(alignment, pos_list, genes)
	# 	# counts for Serine, Leucine encoded by CTG and in general (regardless codon)
	# 	ser, leu, total = 0, 0, 0
	# 	genes.each do |key, gene|
	# 		dna_seq = extract_gene_seq(gene["gene"])
	# 		ref_codons = dna_seq.scan(/.{1,3}/)
	# 		pos_list.each do |pos|
	# 			ref_seq = alignment[key]
	# 			ref_aa = ref_seq[pos]
	# 			spos = alignment_pos2sequence_pos(pos, ref_seq) # position in unaligned sequence
	# 			if ref_aa == "S" && ref_codons[spos] == "CTG" then
	# 				ser += 1
	# 				total += 1
	# 			elsif ref_aa == "L" && ref_codons[spos] == "CTG"
	# 				leu += 1
	# 				total += 1
	# 			elsif ref_aa == "S" || ref_aa == "L"
	# 				# any other codon 
	# 				total += 1
	# 			end # ser/leu count
	# 		end # pos_list.each
	# 	end # genes.each

	# 	# convert to relative frequencies
	# 	pct_ser = ser.to_f / total
	# 	pct_ser = 0.to_f if pct_ser.nan?
	# 	pct_leu = leu.to_f / total
	# 	pct_leu = 0.to_f if pct_leu.nan?
	# 	return pct_ser.round(2), pct_leu.round(2)
	# end


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


#### Evaluation methods
	# evaluate reference data
	# get statistics about ALL CTG positions in ref data - takes some time!
	# for all proteins
		# for all genes
			# - CTG positions (per protein and per organism)
			# - amino acid usage of other genes at this position
			# - CTGs in other genes at this position
	# also writes: number of CTG codons per organism

	def sp2sp_abbr
		fh = File.new("/fab8/smuehlh/data/cugusage/sp_sp-abbr.csv", "w")
		ref_data, fatal_error = load_ref_data
		list = {}
		if fatal_error.empty? then
			ref_data.keys.each do |prot|
				ref_data[prot]["genes"].keys.each do |k|
					org = ref_data[prot]["genes"][k]["species"]
					if ! list.has_key?(org) then
						puts k
						list[org] = k.match(/([A-Z][_a-z]+)[A-Z]\w*/)[1]
					end # ! list.has_key?(org)
				end # ref_data[prot]["genes"].each
			end # ref_data.each
		end # if fatal_error.empty?

		list.each do |sp, sp_abbr|
			fh.puts sp + "," + sp_abbr
		end

		fh.close

		render :eval_ref_data, formats:[:js]
	end

def mean_protein_length
	ref_data, fatal_error = load_ref_data

	if fatal_error.empty? then
		ref_data.keys.sort_by { |key| key.to_s.naturalized }.each do |prot|
			puts prot
			ref_alignment = Hash[*ref_data[prot]["alignment"].split("\n")]
			ref_alignment.keys.each {|k| ref_alignment[ k.sub(">", "") ] = ref_alignment.delete(k)}
			lens = ref_alignment.values.collect {|v| v.gsub("-", "").length}
			n_seqs = lens.size
			puts (lens.sum.to_f/n_seqs).round
		end
	end
	render :eval_ref_data, formats:[:js]
end
	def eval_refdata_prot_pro_species
		# done: 20 august

		ref_data, fatal_error = load_ref_data
		fh = File.new("/fab8/smuehlh/data/cugusage/stats_2008/prots_per_org.csv", "w")
		orgs = {}
		if fatal_error.empty? then
			ref_data.keys.sort_by {|key| key.to_s.naturalized}.each do |prot|
				puts prot
				ref_data[prot]["genes"].each do |name, gene|
					if orgs.has_key?(gene["species"]) then
						if orgs[gene["species"]].has_key?(prot) then
							orgs[gene["species"]][prot] += 1
						else
							orgs[gene["species"]][prot] = 1
						end
					else
						orgs[gene["species"]] = {}
						orgs[gene["species"]][prot] = 1
					end
				end
			end
			fh.puts "Species,#{orgs.keys.collect{|sp| orgs[sp].keys}.flatten.uniq.sort.join(",")}"
			allprots = orgs.keys.collect{|sp| orgs[sp].keys}.flatten.uniq.sort
			orgs.keys.each do |sp|
				line_array = Array.new(allprots.size)
				# add every protein count to corresponding field in line_array
				orgs[sp].keys.sort_by {|key| key.to_s.naturalized}.each do |prot|
					ind = allprots.index(prot)
					line_array[ind] = orgs[sp][prot]
				end
				fh.puts "#{sp},#{line_array.join(",")}"
			end
		end
		fh.close
		render :eval_ref_data, formats:[:js]
	end

	def eval_ref_data2
# Achtung: geht nur haeppechenweise mit use_not, next und break !!!
		fh = File.new("/fab8/smuehlh/data/cugusage/stats_2008/statistics_about_ref_data.txt", "a")

		fh_org = File.new("/fab8/smuehlh/data/cugusage/stats_2008/ctg_org.csv", "w")
		fh_org.puts "Organism,# CTGs"

		ref_data, fatal_error = load_ref_data

		if fatal_error.empty? then
			results_per_org = {}
			use_not = [] #["Actin related protein",
			# 	"Actin related protein Class 1",
			# 	"Actin related protein Class 2",
			# 	"Actin related protein Class 3",
			# 	"Actin related protein Class 4",
			# 	"Actin related protein Class 5",
			# 	"Actin related protein Class 6",
			# 	"Actin related protein Class 7",
			# 	"Actin related protein Class 8",
			# 	"Actin related protein Class 9",
			# 	"Actin related protein Class 10",
			# 	"Calcineurin",
			# 	"Calmodulin",
			# 	"Capping Protein",
			# 	"Centrin",
			# 	"Coronin",
			# 	"Dynactin1 p150",
			# 	"Dynactin2 p50",
			# 	"Dynactin3 p24",
			# 	"Dynactin4 p62",
			# 	"Dynactin5 p25",
			# 	"Dynein Heavy Chain",
			# 	"Dynein Intermediate Chain",
			# 	"Dynein Light Chain 8",
			# 	"Dynein Light Intermediate Chain",
			# 	"Dynein TcTex",
			# 	"Formin",
			# 	"Frequenin",
			# 	"Kinesin Class 1",
			# 	"Kinesin Class 3",
			# 	"Kinesin Class 4",
			# 	"Kinesin Class 5",
			# 	"Kinesin Class 8",
			# 	"Kinesin Class 14",
			# 	"Kinesin Class 15",
			# 	"Kinesin Class 16",
			# 	"Myosin essential light chain"]
			# "Myosin heavy chain Class 1"

			regex = Regexp.union(*use_not)

			ref_data.keys.sort_by { |key| key.to_s.naturalized }.each do |prot|

				if regex.match(prot) then
					puts "Skipping #{prot}"
					next
				end
				# break if prot == "Myosin heavy chain Class 1"
				puts prot


				# statistics about all CTG-Positions in reference data
				results = {ref_seq_num: "", ctg_pos: [], ref_chem: {}, ref_ctg: {}, sc_c_aa: {}, ca_a_aa: {}}

				pos2key = {} # map from ctg position to a sequece containing this ctg

				# extract ref_alignment 
				begin
				ref_alignment = Hash[*ref_data[prot]["alignment"].split("\n")]
				ref_alignment.keys.each {|k| ref_alignment[ k.sub(">", "") ] = ref_alignment.delete(k)}
				rescue => e
					puts "---"
					puts "ERROR (#{prot}):"
					puts e
					puts "---"
					next
				end

				results[:ref_seq_num] = ref_alignment.keys.size

				# get all CTG positions in all reference sequences
				ref_data[prot]["genes"].keys.each do |k|

					dna_seq = extract_gene_seq(ref_data[prot]["genes"][k]["gene"])
					codons = dna_seq.scan(/.{1,3}/)
					ctg_pos = codons.each_with_index.map{ |ele, ind| (ele == "CTG") ? ind : nil }.compact

					# store ctg position
					# results[:ctg_pos] |= ctg_pos # set union
					# !!! 
					# pos should be position in reference alignments
					ctg_pos_mapped = []
					ctg_pos.each do |pos|
						ctg_pos_mapped << sequence_pos2alignment_pos(pos, ref_alignment[k])
					end
					ctg_pos = ctg_pos_mapped
					results[:ctg_pos] |= ctg_pos                         

					# store sequence key if it is a "new" ctg position
					ctg_pos.each {|pos| pos2key[pos] = k if ! pos2key.has_key?(pos)}

					# store ctg pos per organism
					org = ref_data[prot]["genes"][k]["species"]
					if ! results_per_org.has_key?(org) then
						results_per_org[org] = {}
					end
					results_per_org[org][prot] = ctg_pos

				end


				# get list of all amino acids and codons at all CTG -Positions
				ref_alignment_cols = []
				ref_codons = []

				results[:ctg_pos].each do |pos|
					# all amino acids and all codons for a given CTG position
					col = []
					codons = []
					# amino acids used by c.albicans and s.cerevisiae at this CTG position, default if prot is not encoded: ""
					ca_a_aa = ""
					sc_c_aa = ""
					
	 				# get reference key of a sequence containing this CTG position
	 				ref_key = pos2key[pos]
	 				# !!!
	 				# pos is now the position in reference alignments
					# ref_cymopos = sequence_pos2alignment_pos(pos, ref_alignment[ref_key])
					ref_cymopos = pos

					ref_alignment.each do |cymo_header, cymo_prot|

						# collect cymo column at CTG position
						col << cymo_prot[ref_cymopos]

						# collect codons at CTG position
						if cymo_prot[ref_cymopos] == "-" then
							# add nil as placeholder if its a gap
							codons << nil
						else
							# actually collect codon
							spos = alignment_pos2sequence_pos(ref_cymopos, cymo_prot) # position in sequece without gaps
							codons << get_codon(ref_data[prot]["genes"][cymo_header]["gene"], spos)
						end

						# collect amino acid used by reference species
						if cymo_header.match(/^(Sc_c)[A-Z]/) then
							sc_c_aa = cymo_prot[ref_cymopos]
						end
						if cymo_header.match(/^(Ca_a)[A-Z]/) then
							ca_a_aa = cymo_prot[ref_cymopos]
						end
					end

					# store results
					ref_alignment_cols << col
					ref_codons << codons.compact # if only nil- placeholders, it will be an empty array

					# amino acid used by s.cerevisiae and c.albicans
					results[:sc_c_aa][pos] = sc_c_aa
					results[:ca_a_aa][pos] = ca_a_aa
				end

				# amino acid usage 
				ref_alignment_cols.each_with_index do |col, ind|
					# check for each column (= each CTG position) the amino acid distribution in reference data
					aa_freq, aa_num = word_frequency(col)
					is_significant, prob_transl = predict_translation(aa_freq)

					# possibly, here the mapping on the CTG position is wrong (results[:ctg_pos][ind])
					results[:ref_chem][results[:ctg_pos][ind]] = {transl: prob_transl, aa_comp: aa_freq, aa_num: aa_num}

				end

				# CTG usage
				results[:ref_ctg] = ref_ctg_usage(ref_codons, ref_alignment_cols, results[:ctg_pos])

				fh.puts prot
				fh.puts "#{results[:ref_seq_num]} sequences."
				fh.puts "#{results[:ctg_pos].size} CTG positions."
				fh.puts ""
				fh.puts "CTG position\tDistribution of amino acids\tNumber of amino acids\tCTG usage\tNumber of CTG codons\tUsage Sc_c\t Usage Ca_c"
				results[:ctg_pos].sort.each do |pos|
					# check if for all ctg positsions info about aa distribution and ctg usage exists
					if results[:ref_chem].has_key?(pos) then
						# call method from PredictionsHelper (formerly with: view_context.format_aa_stats, without extend module)
						aa_stat = view_context.format_aa_stats(results[:ref_chem][pos][:aa_comp])
						n_aas = results[:ref_chem][pos][:aa_num]
					else 
						aa_stat = "Not discriminative."
						n_aas = "-"
					end
					if results[:ref_ctg].has_key?(pos) then
						ctg_stat = view_context.format_aa_stats(results[:ref_ctg][pos][:ctg_usage])
						n_ctgs = results[:ref_ctg][pos][:ctg_num]
					else
						ctg_stat = "Not discriminative."
						n_ctgs = "-"
					end
					if results[:sc_c_aa].has_key?(pos) then
						sc_c_aa = results[:sc_c_aa][pos]
					else
						sc_c_aa = ""
					end

					if results[:ca_a_aa].has_key?(pos) then
						ca_a_aa = results[:ca_a_aa][pos]
					else
						ca_a_aa = ""
					end
					# tab separated list of all statistics for this position
					# pos + 1 to convert between ruby and human counting! 
					fh.puts "#{view_context.ruby2human_counting(pos)}\t#{aa_stat}\t#{n_aas}\t#{ctg_stat}\t#{n_ctgs}\t#{sc_c_aa}\t#{ca_a_aa}"
				end
				fh.puts ""
				# "free" some memory
				results = nil
				ref_alignment_cols = nil
				ref_alignment = nil

			end # ref_data.each

			# store results per org to file
			results_per_org.each do |org, prots|
				sum = 0

				# uniq positions for actins, myosins, kinesins, tubulins and capping prots!
				# return total count with uniq positions
				sum += prots.select{|k,_| k.match("Actin")}.values.flatten.uniq.size
				sum += prots.select{|k,_| k.match("Capping")}.values.flatten.uniq.size
				sum += prots.select{|k,_| k.match("Kinesin")}.values.flatten.uniq.size
				sum += prots.select{|k,_| k.match("Myosin")}.values.flatten.uniq.size
				sum += prots.select{|k,_| k.match("Tubulin")}.values.flatten.uniq.size

				# delete these counts from hash and add the remaining to sum
				prots.delete_if{|k,_| k.match("Actin")}
				prots.delete_if{|k,_| k.match("Capping")}
				prots.delete_if{|k,_| k.match("Kinesin")}
				prots.delete_if{|k,_| k.match("Myosin")}
				prots.delete_if{|k,_| k.match("Tubulin")}

				sum += prots.values.flatten.size
				fh_org.puts org + "," + sum.to_s
			end

		end # if fatal_error.empty?

		# "free" some memory
		ref_alignment_cols = nil
		ref_data = nil
		ref_alignment = nil

		fh.close
		fh_org.close
		# stat_conserved_pos
		# stats_ctg_prot
		# parse_all_stats
		render :eval_ref_data, formats:[:js]
	end

	# number of CTG codons per protein

	# def stats_ctg_prot
	# 	fh = out_ctg_prot = "/fab8/smuehlh/data/cugusage/ctg_prot_neu.csv"
	# 	File.open(fh).each do |line|
	# 		next if line.include?("Protein")
	# 		line.chomp!
	# 		parts = line.split("\;")
	# 		puts "#{parts[0]},#{(parts[1].to_f/parts[2].to_f).round(2)}"
	# 	end
	# 	render :eval_ref_data, formats:[:js]
	# end
	def stats_ctg_prot 
		file_in = "/fab8/smuehlh/data/cugusage/stats_2008/statistics_about_ref_data.txt"

		# out-file 1: Verteilung CTG Positionen pro Protein
		out_ctg_prot = "/fab8/smuehlh/data/cugusage/stats_2008/ctg_prot_neu.csv"
		fh_ctg_prot = File.new(out_ctg_prot, "w")
		fh_ctg_prot.puts "Protein,# CTGs"

		prot = ""
		# n_ctg_a = 0
		# n_ctg_k = 0
		# n_ctg_m = 0
		File.open(file_in).each do |line|
			if ! line.include?(".") && ! line.include?("\t") then
				# its an protein
				# if line.include?("Actin") then
				# 	prot = "Actin"
				# elsif line.include?("Kinesin")
				# 	prot = "Kinesin"
				# elsif line.include?("Myosin heavy chain")
				# 	prot = "Myosin heavy chain"
				# else
					prot = line.chomp
				# end
			elsif line.include?("CTG positions.") 
				# its the number of CTG positions
				n_ctg = line.match(/\d+/)[0]
				# if prot == "Actin" then
				# 	n_ctg_a += n_ctg.to_i

				# elsif prot == "Kinesin"
				# 	n_ctg_k += n_ctg.to_i

				# elsif prot == "Myosin heavy chain"
				# 	n_ctg_m += n_ctg.to_i

				# else
					fh_ctg_prot.puts prot + "," + n_ctg
				# end

			end
		end
		# fh_ctg_prot.puts "Actin,"+n_ctg_a.to_s
		# fh_ctg_prot.puts "Kinesin,"+n_ctg_k.to_s
		# fh_ctg_prot.puts "Myosin heavy chain,"+n_ctg_m.to_s

		fh_ctg_prot.close

		render :eval_ref_data, formats:[:js]
	end


	# number of CTG codons per protein and number of conserved CTG codons (occuring in min. 2 and min. 5 genes) per protein
	# number of CTG codons broken down: number of genes vs. number of ctgs occuring in that many genes  
	def stat_conserved_pos

		file_in = "/fab8/smuehlh/data/cugusage/stats_2008/statistics_about_ref_data.txt"

		out_conserved = "/fab8/smuehlh/data/cugusage/stats_2008/prot_conserved_pos.csv"
		fh_out = File.new(out_conserved, "w")
		fh_out.puts "Protein,#CTGs,#CTGs in min.2,#CTGs in min.5,#CTG je 1000 aa,#aa,#CTG=Ser"

		# out_more_det = "/fab8/smuehlh/data/cugusage/stats_2008/prot_conserved_pos_detail.csv"
		# fh_out_detailed = File.new(out_more_det, "w")
		# fh_out_detailed.puts "Protein"
		# fh_out_detailed.puts "#genes,#CTGs"

		prot = ""
		n_ctg = 0
		n_conserved = 0
		n_conserved_in5 = 0
		n_used_as_ser = 0
		n_aas = 0
		n_ctg_per_1000aa = 0 # = n_ctg/n_aas * 1000
		# pos_n_genes = Hash.new(0) # use this hash to extract number of CTGs occuring in 1, 2, 3, ... genes (the detailed csv)

		File.open(file_in).each do |line|
			line.chomp!
			parts = line.split("\t")
			n_ctgs = parts[4]
			ctg_usage = parts[3]

			if ! line.include?(".") && ! line.include?("\t") && ! line.blank? then
				# its a protein

				# write data of last visited protein, if this is not the very first protein
				if ! prot.blank? then
					n_ctg_per_1000aa = (n_ctg.to_i / n_aas.to_f * 1000).round
					fh_out.puts prot + "," + n_ctg + "," + n_conserved.to_s + "," + n_conserved_in5.to_s + "," + n_ctg_per_1000aa.to_s + "," + n_aas.to_s + "," + n_used_as_ser.to_s

					# # write detailed output
					# fh_out_detailed.puts prot
					# pos_n_genes.keys.sort.each do |key|
					# 	fh_out_detailed.puts key.to_s + "," + pos_n_genes[key].to_s
					# end
					# fh_out_detailed.puts ""

				end

				# set up data for this protein
				prot = line
				puts prot
				n_ctg = 0
				n_conserved = 0
				n_conserved_in5 = 0
				# pos_n_genes = Hash.new(0)
				n_used_as_ser = 0
				n_aas = 0
				n_ctg_per_1000aa = 0

				# total number of amino acids in ref data: parse alignment file
				prot_mapped = prot.gsub(" ", "-").downcase
				Dir.glob(BASE_PATH + prot_mapped + ".fasta") do |file|
					headers, seqs = fasta2str(File.read(file))
					n_aas = seqs.inject(0) {|sum, seq| sum += seq.gsub("-", "").length}
				end

			elsif line.include?("CTG positions.") 
				# its the total number of CTG positions in this protein
				n_ctg = line.match(/\d+/)[0]
					puts line	
			elsif ! n_ctgs.nil? 
				# its a row with ctg usage data

				n_ctgs = n_ctgs.to_i # a string will be converted to zero, which is ok for now

				if n_ctgs > 0 then
					# number CTGs occuring in xx genes
					# pos_n_genes[n_ctgs] += 1
				end

				if n_ctgs > 1 then
					# its conserved, as it occurs in at least one prot
					n_conserved += 1
				end

				if n_ctgs > 4 
					# its also high conserved, as it occures in at least 5 prots
					n_conserved_in5 += 1
				end

				# get number of CTG codons translated only as Serine
				if ctg_usage =~ /S:\s?(\d+)/ then
					n_used_as_ser += 1
				end

			end
		end

		# write data of last protein
		n_ctg_per_1000aa = (n_ctg.to_i / n_aas.to_f * 1000).round
		fh_out.puts prot + "," + n_ctg + "," + n_conserved.to_s + "," + n_conserved_in5.to_s + "," + n_ctg_per_1000aa.to_s + "," + n_aas.to_s + "," + n_used_as_ser.to_s
	
		# write detailed output
		# fh_out_detailed.puts prot
		# pos_n_genes.keys.sort.each do |key|
		# 	fh_out_detailed.puts key.to_s + "," + pos_n_genes[key].to_s
		# end

		fh_out.close
		# fh_out_detailed.close

		render :eval_ref_data, formats:[:js]
	end

end
