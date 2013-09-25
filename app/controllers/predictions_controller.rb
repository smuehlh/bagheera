# load modules and classes
require 'helper.rb'

class PredictionsController < ApplicationController

	MAX_SIZE = 26214400 # 25 MB

	# Render start page for prediction
	def search
		delete_old_data # delete old uploaded genome files
		prepare_new_session # a fresh session
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
			file_path = File.join(BASE_PATH, file_id, "query.fasta")
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
		example_path = Rails.root + "lib/" + example_file
		@fatal_error = []
		if File.exist?(example_path) then
			# do not check content, but check file size, just to be sure
			@fatal_error |= check_filesize(File.size(example_path))
			if @fatal_error.blank? then
				file_id = rand(1000000000).to_s
				file_path = File.join(BASE_PATH, file_id, "query.fasta")
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
		require 'progCall.rb'
		require 'proteinFamily.rb'
		require 'refProtein.rb'
		require 'prediction.rb'
		# require 'status.rb' # for some reason, this class needs not to be loaded separately

		require 'open3'

		# add options to this session
		ProgCall.blast_filtering = "m S" if params.has_key?(:blast) # use low complexity filter
		Prediction.class_variable_set(:@@align_method, params[:algo])
		Prediction.class_variable_set(:@@align_config, params[:config])

		if session[:file][:name] == "Candida_albicans_WO_1.fasta" then
			ProteinFamily.class_variable_set(:@@ref_data_path, File.join(BASE_PATH,PATH_REF_WO_EXAMPLE) )
		else
			ProteinFamily.class_variable_set(:@@ref_data_path, BASE_PATH)
		end 

		# initialize "results" variables
		@fatal_error = []
		@minor_error = []
		@predicted_prots = {} # Hash of hashes
		@stats = Hash.new(0)

		ref_data = Helper.load_ref_data
		file_basename = File.dirname(session[:file][:path])

		# sucessfully loaded reference data file
		# setup blast database
		ProgCall.genome_db = session[:file][:path].gsub(".fasta", "_db")
		ProgCall.create_blast_db(session[:file][:path])

		# predict genes for every protein

		sorted_ref_prots = ref_data.keys.sort_by {|key| key.to_s.naturalized}
		# Important: do NOT use in_processes here, as then @stats and @predicted_prots will be lost

		Parallel.each( sorted_ref_prots, :in_threads => 10 ) do |prot|

			# initialize 
			all_prot_data = ref_data[prot]

			begin
				prot_obj = RefProtein.new(prot, all_prot_data, file_basename)
			rescue RuntimeError => exp
				# only profile missing or fasta file? 
				# nothing to do if fasta is missing
				next if ! exp.message.include?(".prfl") 
			end

			pred_obj = Prediction.new(prot_obj, 1, file_basename) # use BLAST hit nr 1
			pred_obj.predict

			# save results
			pred_obj.save(@predicted_prots)

			# save error messages
			if ! pred_obj.err_msg.blank? then
				if pred_obj.err_msg.include?("failed") then
					# messages like program failures go into @minor_error
					write_minor_error(prot, pred_obj.err_msg)
				else
					# messages concerning prediction itself go into :message
					pred_obj.save_message(@predicted_prots)
				end
			end

			if ! pred_obj.used_prfl then
				pred_obj.save_message(@predicted_prots, "No protein profile available for gene prediction." )
			end

			Status.update(@predicted_prots[prot], @stats)

		end # Parallel

		# save status to use it again in predict_more
		f_stats = File.join(file_basename, "stat")
		Status.save(f_stats, @stats)

	rescue RuntimeError => exp
		@fatal_error = [exp.message]

	rescue NoMethodError, TypeError => exp
		@fatal_error = ["Sorry, an error occured. Please contact us."]

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

		require 'progCall.rb'
		require 'proteinFamily.rb'
		require 'refProtein.rb'
		require 'prediction.rb'
		# require 'status.rb' # for some reason, this class needs not to be loaded separately

		require 'open3'

		# add options from view 

		ProgCall.blast_filtering = "m S" if params.has_key?(:blast) # use low complexity filter
		Prediction.class_variable_set(:@@align_method, params[:algo]) if params.has_key?(:algo)
		Prediction.class_variable_set(:@@align_config, params[:config]) if params.has_key?(:config)
		hit_start = params[:hit].to_i + 1 # use next hit for prediction
		@prot = params[:prot]
		file_basename = File.dirname(session[:file][:path]) 

		if session[:file][:name] == "Candida_albicans_WO_1.fasta" then
			ProteinFamily.class_variable_set(:@@ref_data_path, File.join(BASE_PATH,PATH_REF_WO_EXAMPLE) )
		else
			ProteinFamily.class_variable_set(:@@ref_data_path, BASE_PATH)
		end 

		# initialize "results" variables
		@fatal_error = []
		@minor_error = []
		@predicted_prots = {} # Hash of hashes
		@stats = Hash.new(0)

		# initialize reference data
		ref_data = Helper.load_ref_data
		all_prot_data = ref_data[@prot]
		begin
			prot_obj = RefProtein.new(@prot, all_prot_data, file_basename)
		rescue RuntimeError => exp
			# only profile missing or fasta file? 
			# cannot do anything if fasta is missing
			path_to_fasta = File.join(
				ProteinFamily.class_variable_get(:@@ref_data_path), "#{prot_obj.prot_basename}.fasta"
				)
			Helper.file_exist_or_die(path_to_fasta)
		end

		n_hits = Prediction.new(prot_obj, hit_start, file_basename).sneak_n_blast_hits
		hit_stop = hit_start == 2 ? [10, n_hits].min : [hit_start+9, n_hits].min # go in steps of 10, or to last hit

		# predict genes for every blast hit

		# Important: do NOT use in_processes here, as then @stats and @predicted_prots will be lost
		Parallel.each( (hit_start..hit_stop), :in_threads => 10 ) do |n_hit|
			# initialize 
			key = @prot + "_" + n_hit.to_s

			pred_obj = Prediction.new(prot_obj, n_hit, file_basename) 

			pred_obj.predict

			# save results
			pred_obj.save(@predicted_prots, key)

			# save error messages
			if ! pred_obj.err_msg.blank? then
				if pred_obj.err_msg.include?("failed") then
					# messages like program failures go into @minor_error
					write_minor_error(key, pred_obj.err_msg)
				else
					# messages concerning prediction itself go into :message
					pred_obj.save_message(@predicted_prots, key)
				end
			end

			if ! pred_obj.used_prfl then
				pred_obj.save_message(@predicted_prots, key, "No protein profile available for gene prediction." )
			end

			Status.update(@predicted_prots[key], @stats)

		end # Parallel

		# save status to use it again in predict_more
		f_stats = File.join(file_basename, "stat")
		Status.save(f_stats, @stats)

	rescue RuntimeError => exp
		@fatal_error = [exp.message]

	rescue NoMethodError, TypeError => exp
		@fatal_error = ["Sorry, an error occured. Please contact us."]

		render :predict_more, formats: [:js]

	end

	# prepare a new session 
	def prepare_new_session
		reset_session
	    session[:file] = {} 		# uploaded genome described by keys: :name, :id, :path
	    # session[:align] = {}		# alignment options for seqan: pair_align or DIALIGN
	    # session[:augustus] = {}		# options for augustus gene prediction
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

	# def write_status(done, total, final=false)
	# 	html = ["<html>", "<head>"]
	# 	html += ['<meta http-equiv="refresh" content="5;url=./status.html">'] if ! final 
	# 	html += ["</head>", "<body>", "Processed #{done}/#{total}", "</body>", "</html>"]
	# 	return html.join(" ")
	# end


	# compare uploadedfile-size with maximal size (50MB)
	# @param size [Fixnum] File size 
	# @return [Array] Error message if size exceeds MAX_SIZE
	def check_filesize(size)
		errors = []
		if size > MAX_SIZE then
			errors << "Fatal error."
			errors << "File must be less than 25 MB."
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

	# build up @minor_error [Array] Error messages
	# non-fatal errors, affecting only a single protein
	# on fab8, verbose error messages will be printed, outside non-verbose ones
	# @param prot [String] Protein for which gene prediction failed
	# @param message [String] More detailed information about the error
	def write_minor_error (prot, message)
		@minor_error << prot + ": " + message
	end


end
