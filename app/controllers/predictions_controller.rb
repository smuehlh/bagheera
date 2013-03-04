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
	HYDORPHOBIC_AAS = ["V", "I", "L"]
	POLAR_AAS = ["S", "T"]

	# Render start page for prediction
	def search
		delete_old_data # delete old uploaded genome files
		prepare_new_session # a fresh session
	end

	# Handel uploaded genome file: 
	# store it in /tmp/cug and add its id(new file name), original filename and path to the session
	# render partial showing the uploaded file
	# accessible params: @content_error [Array] Errors occured during file load
	def upload_file
		`fromdos #{params[:uploaded_file].path}` # remove \r, if uploaded from windows
		# is the file too big? maybe jquery-check did not work, so better check again
		@content_error = check_filesize(params[:uploaded_file].size)
		# check content, only if file is not too big
		@content_error |= check_fasta(params[:uploaded_file]) if @content_error.blank?
		if @content_error.blank? then
			# rename and save file
			file_id = rand(1000000000).to_s
			file_name = File.basename(params[:uploaded_file].original_filename)
			file_path = BASE_PATH + "uploaded_genomefile_" + file_id + ".fasta"
			File.rename(params[:uploaded_file].path, file_path)
			File.open(file_path, 'wb') {|file| file.write(params[:uploaded_file].read)}
			session[:file] = { id: file_id, name: file_name, path: file_path }
		end
		render :upload_file_ajax, formats: [:js]
	end

	# Handel loading the example genome:
	# store it in /tmp/cug and add its id(new file name), original filename and path to the session
	# render partial showing the uploaded file
	# accessible params: @content_error [Array] Errors occured during file load
	def load_example
		# do not check content, but check file size, just to be sure
		example_file = "Candida_albicans_WO_1.fasta"
		example_path = Rails.root + "spec/fixtures/files/" + example_file
		@content_error = check_filesize(File.size(example_path))
		if @content_error.blank? then
			file_id = rand(1000000000).to_s
			file_path = BASE_PATH + "uploaded_genomefile_" + file_id + ".fasta"
			File.copy_stream(example_path, file_path)
			@file = example_file
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
		ext = session[:align][:algo] == "dialign" ? "_aligned.fasta.fa" : "_aligned.fasta"
		file_scr = BASE_PATH + prot + "_" + hit + "_" + session[:file][:id] + ext
		file_dest = Dir::tmpdir + "/cymobase_alignment_" + @file_id + ".fasta"
		FileUtils::cp(file_scr, file_dest)
		# TODO use complete cymoalignment for prot instead of only one ref seq?
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
	# accessible params in view: @errors[Array] Errors occured during prediction
	def predict_genes
		# general workflow:
		# 1) extract reference proteins
		# 2) gene prediction foreach reference protein:
		# 2.1) BLAST
		# 2.2) AUGUSTUS
		# 3) Compare with reference data
		# 3.1) Compare with reference alignment
		# 3.2) Compare with reference genes

		# add options to session
		session[:align] = { algo: params[:algo], config: params[:config] }
		session[:augustus] = { species: params[:species] }

		ref_data, @errors = load_ref_data
		@predicted_prots = {} # containing final results ...
		if @errors.empty? then
			# sucessfully loaded reference data file
			# setup blast database
			genome_db = session[:file][:path].gsub(".fasta", "_db")
			error = create_blast_db(genome_db)
			if ! error.blank? then
				# database setup failed; no cug-usage can be predicted
				@errors << error
				return
			end
			
			# all other tasks need to be done for each ref protein
			# peach for paralell execution, maximum 10 threads
			ref_data.peach(10) do |prot, all_prot_data|
				@predicted_prots[prot] = {} # ... final results for each protein!

				### 1) extract reference proteins
				ref_prot_species, ref_prot_stat, ref_prot_seq, ref_prot_key = "", "", "", ""

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
				# fill other ref_data-types
				ref_prot_species = all_prot_data["genes"][ref_prot_key]["species"]
				ref_prot_stat = all_prot_data["genes"][ref_prot_key]["completeness"]
				ref_prot_seq = all_prot_data["genes"][ref_prot_key]["gene"].match(/prot_seq: (.*?)\n/)[1]
#				ref_prot_geneseq = extract_gene_seq(all_prot_data["genes"][ref_prot_key]["gene"])

				prot_basename = prot.gsub(" ", "-").downcase

				# save reference protein in fasta format in file (needed for blast alignment of reference and predicted protein)
				seq_file = BASE_PATH + prot_basename + "_" + session[:file][:id] + ".fasta"
				fasta = str2fasta(ref_prot_key, ref_prot_seq)
				File.open(seq_file, 'w'){|file| file.write(fasta)}

				### 2) gene prediction foreach reference protein:

				### 2.1) BLAST
				# 2.1.1) convert uploaded genome into blast-db
				# 2.1.2) blast-query
				blast_file = BASE_PATH + prot_basename + "_" + session[:file][:id] + ".blastout"  # store blast hits
						# -p [PROGRAM] protein query against nt-db -d [DATABASE] -i [FILE] -m8 [OUTPUT FORMAT] -F [FILTERING] -s [SMITH-WATERMAN ALIGNMENT] T 
				stdin, stdout_err, wait_thr = Open3.popen2e(BLASTALL, "-p", "tblastn", "-d", genome_db, "-i", seq_file, "-m8", "-F", "m S", "-s", "T")
				stdin.close
				output = stdout_err.read
				stdout_err.close
				# write output to also to file to store it for later gene prediction
				File.open(blast_file, 'w') {|f| f.write(output)}

				if ! wait_thr.value.success? || output.include?("ERROR") then
					@errors << "#{prot}: Gene prediction (BLAST) failed"
					@predicted_prots[prot][:message] = "Sorry, an error occured"
					# delete file with blast hits, its contains only error messages
					File.delete(blast_file)
					next
				end
				@predicted_prots[prot][:n_hits] = output.lines.to_a.size
				@predicted_prots[prot][:hit_shown] = 1

				is_success, pred_seq, pred_dnaseq = perform_gene_pred(output)
				if ! is_success then
					# an error occured
					# variable called "pred_seq" actually contains the program whith caused the error
					@errors << "#{prot}: Gene prediction (#{pred_seq}) failed"
					@predicted_prots[prot][:message] = "Sorry, an error occured"
					next
				end
		
				is_success, ref_seq_aligned, pred_seq_aligned = align_pred_ref(seq_file, pred_seq)
				if ! is_success then
					# an error occured
					# variable called "ref_seq_aligned" actually contains the program which caused the error 
					@errors << "#{prot}: Gene prediction (#{ref_seq_aligned}) failed"
					@predicted_prots[prot][:message] = "Sorry, an error occured"
					next
				end
			
				is_success, results, message = compare_pred_gene(ref_data, ref_prot_key, prot, ref_seq_aligned, pred_seq_aligned, pred_dnaseq)
				if ! is_success then
					@predicted_prots[prot][:message] = message
				else
					@predicted_prots[prot][:message] = ""
				end
				@predicted_prots[prot].merge!(results)

			end # ref_data.peach do |prot, des|
		end # if @errors.empty?
		@stats = calc_stats(@predicted_prots)

		render :predict_genes, formats: [:js]
	end

	# redo gene prediction for next 10 blast hits (gene prediction for best blast hit is already done) for one specific protein
	# same workflow as predict_genes
	# renders predict_more
	# accessible params in view: @predicted_prots [Hash] Prediction data for reference protein @prot
	# accessible params in view: @stats [Hash] Up-to-date statistics over prediction data
	# accessible params in view: @prot [String] Protein for which gene prediction was restarted
	# accessible params in view: @errors[Array] Errors occured during prediction
	def predict_more
		@prot = params[:prot]
		hit_start = params[:hit].to_i + 1 # use next hit for prediction

		# prepare stuff
		@predicted_prots = {}
		ref_data, @errors = load_ref_data
		basename = @prot.gsub(" ", "-").downcase
		seq_file = BASE_PATH + basename + "_" + session[:file][:id] + ".fasta"
		blast_file = BASE_PATH + basename + "_" + session[:file][:id] + ".blastout"
		headers, seqs = fasta2str(File.read(seq_file))
		ref_prot_key = headers[0]
		blast_all_hits = File.read(blast_file)
		n_hits = blast_all_hits.lines.to_a.size
		hit_stop = hit_start == 2 ? [10, n_hits].min : [hit_start+9, n_hits].min # go in steps of 10, or to last hit

		(hit_start..hit_stop).each do |n_hit|
			key = @prot + "_" + n_hit.to_s
			@predicted_prots[key] = { n_hits: n_hits, hit_shown: n_hit}

			is_success, pred_seq, pred_dnaseq = perform_gene_pred(blast_all_hits, n_hit)

			if ! is_success then
				# an error occured
				# variable called "pred_seq" actually contains the program whith caused the error
				if pred_seq.nil? then
					# simply run out of BLAST hits
					break
				else
					# an "real" error
					@errors << "#{@prot}: Gene prediction (#{pred_seq}) failed"
					@predicted_prots[key][:message] = "Sorry, an error occured"
					next
				end
			end
			
			is_success, ref_seq_aligned, pred_seq_aligned = align_pred_ref(seq_file, n_hit, pred_seq)
			if ! is_success then
				# an error occured
				# variable called "ref_seq_aligned" actually contains the program which caused the error 
				@errors << "#{@prot}: Gene prediction (#{ref_seq_aligned}) failed"
				@predicted_prots[key][:message] = "Sorry, an error occured"
				next
			end

			is_success, results, message = compare_pred_gene(ref_data, ref_prot_key, @prot, ref_seq_aligned, pred_seq_aligned, pred_dnaseq)
			if ! is_success then
				@predicted_prots[key][:message] = message
			else
				@predicted_prots[key][:message] = ""
			end
			@predicted_prots[key].merge!(results)
		end
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
		# delete all files older than one day but file "alignment_gene_structure.json"
		system("find", BASE_PATH, "-type", "f", "-not", "-name", REF_DATA, "-mtime", "+#{days}", "-delete")
		# delete alignment file for show-alignment function
		FileUtils.rm Dir.glob(Dir::tmpdir + '/cymobase_alignment_cug*')
	end

	# load reference data from file alignment_gene_structure.json
	# @return [FalseClass] If an error occured
	# @return [TrueClass] If no error occured -> Array will be empty
	# @return [Array] Errors occured during file load
	def load_ref_data
		path = BASE_PATH + REF_DATA
		if ! FileTest.file?(path) then
			errors = ["Sorry, cannot load reference data. Please contact us!"]
			return false, errors
		end
		return JSON.load(File.read(path)), []
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

		# updata statistics (or create new ones)
		stats["sig_pos"] += data.values.inject(0) {|res, val| val[:has_ctg]? res + val[:ctg_pos].length : res + 0}
		stats["n_prots"] += data.keys.length
		stats["ser"] += data.values.inject(0) {|res, val| val[:has_ctg]? res + val[:ctg_transl].count("S") : res + 0}
		stats["leu"] += data.values.inject(0) {|res, val| val[:has_ctg]? res + val[:ctg_transl].count("L") : res + 0}
		stats["ref_ctgpos"] += data.values.inject(0) {|res, val| 
			val[:has_ctg] ? res + val[:ctg_ref].collect {|hash| hash.values}.flatten.count{|x| x >= 0.5} : res + 0}

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
					last_nucleotides = line[-2,2]
				end
			end
		end

		# file read, check what happend
		if ! is_fasta then
			errors << "Invalid input."
			errors << "Expected fasta formatted genome file."
		end
		if  is_fasta && (! contains_ctg) then
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
		error = ""
			# -i [INPUT: GENOME FILE] -p [PROTEIN] FALSE -n [DB NAME] -o [CREATE INDEX OVER SEQID] TRUE
		is_success = system(FORMATDB, "-i", session[:file][:path], "-p", "F", "-n", genome_db, "-o", "T")
		# no output needed, so system is sufficient
		if ! is_success then
			error = "Gene prediction failed: Cannot create BLAST database"
		end 
		return error
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
	# @return [TrueClass] If no error occured
	# @return [String] Predicted protein sequence
	# @return [String] Predicted coding sequence (DNA)
	# @return [FalseClass] AND NO OTHER VALUES If an error occured
	def run_augustus(seq_file)
			# --species [REFERENCE SPEC] QUERY --genemodel=exactlyone [predict exactly one complete gene] --codingseq=on [output also coding sequence]
			# redirection: add stderr to stdout (screen)
		stdin, stdout_err, wait_thr = Open3.popen2e(AUGUSTUS, "--AUGUSTUS_CONFIG_PATH=/usr/local/bin/augustus/config/", 
			"--species=#{session[:augustus][:species]}", "--genemodel=exactlyone","--codingseq=on", seq_file)
		stdin.close
		output = stdout_err.read
		stdout_err.close
		if ! wait_thr.value.success? || output.include?("ERROR") || output.blank? then
			return false
		end

		# no error, parse augustus
		pred_seq, pred_dnaseq = parse_augustus(output)
		return true, pred_seq, pred_dnaseq
	end


	# performing gene prediction for given protein (step 2.2)
	# calling augustus based on best blast hit 
	# @param blast_hits [String] Complete blast output (Hit list)
	# @param nr [Fixnum] Number of blast hit of interest (human counting -> starts with 1)
	# @return [TrueClass] If no error occured
	# @return [String] Predicted protein sequence
	# @return [String] Predicted coding sequence (DNA)
	# @return [FalseClass] AND NO OTHER VALUES If an error occured
	def perform_gene_pred(blast_hits, blast_hit_nr=1)
		# parse blast output
		is_success, seq_id, start, stop, strand = parse_blast_hits(blast_hits, blast_hit_nr) 
		if ! is_success then
			# number exceeds blast hits!
			return false
		end 

		# 2.2) AUGUSTUS
		# 2.2.1) get matching search sequence 
		file = BASE_PATH + session[:file][:id] + "_" + rand(1000000000).to_s + ".fasta"
		genome_db = session[:file][:path].gsub(".fasta", "_db")
				# -d [DB] -p [PROTEIN] false -s [SEARCH STRING] -L [START,STOP]
				# add 1000 nucleotides to start/stop 
				# redirection: stderr to stdout (at this moment: sreen), stdout to file file
		stdin, stdout_err, wait_thr = Open3.popen2e(FASTACMD, "-d", genome_db, "-p", "F", "-s", seq_id, 
			"-L", "#{start-1000},#{stop+1000}", "-S", strand.to_s, "-o", file)
		stdin.close
		output = stdout_err.read
		stdout_err.close
		if ! wait_thr.value.success? || output.include?("ERROR") then
			return false, "BLAST-FASTACMD"
		end

		is_success, pred_seq, pred_dnaseq = run_augustus(file)
		if ! is_success then
			# gene prediction failed
			return false, "AUGUSTUS"
		end

		return true, pred_seq, pred_dnaseq
	end

	# align reference protein with the predicted protein
	# @param file [String] file path containing reference sequence
	# @param hit [Fixnum] Number of blast hit of interest (human counting -> starts with 1)
	# @param pred_seq [String] Predicted protein sequence
	# @return [TrueClass] If no error occured
	# @return [String] Aligned reference sequence
	# @return [String] Aligned predicted sequence
	# @return [FalseClass] If an error occured
	# @return [String] Method which failed ONLY IN CASE OF AN ERROR
	def align_pred_ref(file, hit=1, pred_seq)

		# # prepare input data for alignment program: needs to contain both ref and pred. seq
		# file_in = BASE_PATH + session[:file][:id] + "_" + rand(1000000000).to_s + ".fasta"
		# file_out = file_in.gsub(".fasta", "_aligned.fasta")
		file_in = file.gsub("_", "_#{hit}_")
		file_out = file_in.gsub(".fasta", "_aligned.fasta")

		FileUtils.cp(file, file_in)
		fasta = str2fasta("Prediction", pred_seq)
		File.open(file_in, 'a') {|file| file.write("\n" + fasta)}

		# use pairalign or dialign
		if session[:align][:algo] == "dialign" then
			# use system, as dialign has no easy-parsable error messages (like containing "ERROR")
			# check exit status for success 
			is_success = system(DIALIGN2, "-fa", "-fn", file_out, file_in)
			if ! is_success then
				return false, "Dialign2"
			end
			file_out = file_out << ".fa" # dialign automatically adds ".fa" to basic output filename
		else
				# -c tttt: initialize first row & column with zeros, search last row & column for maximum
				# 	=> end gap free alignment
				# --matrx matrix file
				# implicit options: gap open penalty: -11, gap extension penalty: -1
				# 					protein sequences
				# don't use system, as pair_align is quite verbose and has an easy parsable success-output
			stdin, stdout_err, wait_thr = Open3.popen2e(PAIR_ALIGN, "--seq", file_in, "--matrix", Rails.root.join('lib', 'blosum62').to_s, 
				"--config", session[:align][:config], "--method", session[:align][:algo], "--outfile", file_out)
			stdin.close
			output = stdout_err.read
			stdout_err.close
			if ! wait_thr.value.success? || ! output.include?("Alignment score") then
				return false, "Seqan::pair_align"
			end
		end

		# parse multiple sequence alignment to get aligned reference and aligned predicted sequence
		# first is reference seq
		# second is predicted seq (human counting)
		dummy, seqs = fasta2str(File.read(file_out))
		return true, seqs[0], seqs[1]
	end

	# preforming cug-usage prediction 
	# step 3.1) Compare with reference alignment
	# step 3.1.1) map predicted gene onto reference protein
	# step 3.2) Compare with reference genes
	# @param ref_data [Hash] Cymobase reference data
	# @param ref_prot_key [String] Reference protein identifier; key in ref_data
	# @param prot [String] Protein name; key in ref_data
	# @param ref_seq_aligned [String] Aligned reference sequence
	# @param pred_seq_aligned [String] Aligned predicted sequence
	# @param pred_dnaseq [String] Predicted coding sequence (DNA)
	# @return [TrueClass] If no error occured
	# @return [Hash] Prediction data for protein prot
	# @return [FalseClass] If an error occured
	# @return [Hash] Prediction data for protein prot ONLY IN CASE OF AN ERROR
	# @return [String] Error message ONLY IN CASE OF AN ERROR
	def compare_pred_gene(ref_data, ref_prot_key, prot, ref_seq_aligned, pred_seq_aligned, pred_dnaseq)
		results = {has_ctg: "", pred_prot: "", ctg_pos: [], ctg_transl: [], aa_comp: [], ctg_ref: []}

		# has predicted protein CTG codons?
		codons = pred_dnaseq.scan(/.{1,3}/)
		if !codons.include?("CTG")
			# nope, no CTG => nothing to predict
			results[:has_ctg] = false
			return false, results, "No CTG in prediced gene"
		else
			results[:pred_prot] = pred_seq_aligned.gsub("-", "")
			results[:has_ctg] = true
		end
		ctg_pos = codons.each_with_index.map{ |ele, ind| (ele == "CTG") ? ind : nil }.compact

		# Compare with reference alignment
		# => map CTG positions in predicted protein onto reference alignment columns
		# => compare with chemical properties of reference alignment columns
		# => Compare with reference genes: count usage of CTG for each mapped L and S 
		cymo_algnmnt = Hash[*ref_data[prot]["alignment"].split("\n")]
		ctg_pos.each do |pos|
			# map CTG positions
			algnmnt_col_aa, algnmnt_col_pos = parse_alignment_by_ctgpos(cymo_algnmnt, pos, pred_seq_aligned, ref_seq_aligned, ref_prot_key)
			if !(algnmnt_col_aa.nil? && algnmnt_col_pos.nil?) then
				# mapping was possible, compare with chemical properties
				aa_freq, aa_num = word_frequency(algnmnt_col_aa)
				seq_num = cymo_algnmnt.keys.length
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
	# @param str [String] sequence
	# @return [String] fasta formatted header and sequence
	def str2fasta(header, str)
		fasta = header.include?(">") ? header << "\n" : ">" << header << "\n"
		# breaking string at every 80th character, joining by newline
		fasta += str.scan(/.{1,80}/).join("\n")
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

	# match CTG position to reference alignment
	# @param algnmnt [Hash] Cymobase alignment
	# @param pred_spos [Fixnum] Position in unaligned predicted sequence
	# @param pred_aseq [String] Aligned predicted sequence
	# @param ref_aseq [String] Aligned reference sequence
	# @param ref_key [String] Reference protein; key in algnmnt
	# @return [Array] Residues at matched alignment column
	# @return [Fixnum] CTG position in reference sequence aligned to cymobase alignment
	# @return [NilClass] If matching not possible (gap in cymobase alignment) AND
	# @return [NilClass] If matching not possible (gap in cymobase alignment)
	def parse_alignment_by_ctgpos(algnmnt, pred_spos, pred_aseq, ref_aseq, ref_key)
		# 1) ctg pos in unaligned seq -> pos in aligned predicted sequence (aligned with reference seq)
		pred_apos = sequence_pos2alignment_pos(pred_spos, pred_aseq) # = pos_ref_seq_aligned

		# matching possible at all?
		# NO: gap in reference sequence at ctg pos in predicted sequence
		# or, in case of dialign: important positions are not aligned at all (= in lowercase letters)
		if (ref_aseq[pred_apos] == "-" || 
			ref_aseq[pred_apos].match(/\p{Lower}/) || pred_aseq[pred_spos].match(/\p{Lower}/)) then
			return nil, nil
		end

		# YES

		# 2) pos aligned -> cymo alignment
		ref_spos = alignment_pos2sequence_pos(pred_apos, ref_aseq)
		ref_cymopos = sequence_pos2alignment_pos(ref_spos, algnmnt[">"<<ref_key])
		# 3) return column
		col = algnmnt.collect {|cymo_prot, cymo_seq| cymo_seq[ref_cymopos]}
		return col, ref_cymopos
	end

	# counting word frequencies; doesnot count "-"
	# @param arr [Arr] Containing words (e.g. amino acids)
	# @return [Hash] Keys are the words, values their relative frequency
	# @return [Fixnum] Sum of all occurrences (e.g. total number of amino acids)
	def word_frequency(arr)
		res = Hash.new(0)
		arr.each { |a| res[a] += 1 }
		res.delete("-") # delete count for gaps in alignment!
		# normalize word frequency
		sum = res.inject(0) {|s, (_,val)| s + val}.to_f
		norm_res = res.each {|k,v| res[k] = v/sum} # normalize
		return norm_res, sum.to_i
	end

	# predict the most probable translation based on amino acids statistics
	# @param aa_freq [Hash] Amino acid frequencies
	# @return [Boolean] Are stats discriminative?
	# @return [String] Most probable CUG-translation ("S" or "L")
	def predict_translation(aa_freq)
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

	# check for a given position in the cymo alignment, how many "S" and how many "L" are encoded by an CTG
	# @param algnmnt [Hash] Cymobase alignment
	# @param apos [Fixnum] CTG position in cymo alignment
	# @param genes [Hash] Genes from reference data
	# @return [Fixnum] Percentage of "Ser" encoded by CTGs (rounded to 2 decimal places)
	# @return [Fixnum] Percentage of "Leu" encoded by CTGs (rounded to 2 decimal places)
	def ref_ctg_usage(algnmnt, apos, genes)
		pct_ctg_S, pct_ctg_L = "", "" 
		# Serine
		ref_codons = []
		algnmnt.collect{|k,v| k if v[apos]=="S"}.compact.each do |ref_key|
			dna_seq = extract_gene_seq(genes[ref_key.gsub(">", "")]["gene"]) # access via "gene-name"
			codons_ref = dna_seq.scan(/.{1,3}/)
			spos = alignment_pos2sequence_pos(apos, algnmnt[ref_key]) # access via ">gene-name"
			ref_codons << codons_ref[spos]
		end
		pct_ctg_S = ref_codons.count("CTG").to_f/ref_codons.length 
		pct_ctg_S = 0.0 if pct_ctg_S.nan?

		# same for Leucine
		ref_codons = []
		algnmnt.collect{|k,v| k if v[apos]=="L"}.compact.each do |ref_key|
			dna_seq = extract_gene_seq(genes[ref_key.gsub(">", "")]["gene"]) # access via "gene-name"
			codons_ref = dna_seq.scan(/.{1,3}/)
			spos = alignment_pos2sequence_pos(apos, algnmnt[ref_key]) # access via ">gene-name"
			ref_codons << codons_ref[spos]
		end
		pct_ctg_L = ref_codons.count("CTG").to_f/ref_codons.length
		pct_ctg_L = 0.0 if pct_ctg_L.nan?
		return pct_ctg_S.round(2), pct_ctg_L.round(2)
	end
end
