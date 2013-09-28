class TreeController < ApplicationController

	def calc_tree

		@selected_prots = {}
		@error = ""
		@id = "cugtree" + rand(1000000000).to_s # for lucullus

		path = File.join(BASE_PATH, session[:file][:id])

		# if mafft was not used to create the alignment, files will contain only predicted seq and one reference seq
		# => useless for tree calucation
		if params[:algo] == "mafft" then
			# will use first prediction (more if predict_more was called)
			tree_obj = Tree.new(params[:n_prot].to_i, path)

			# collect sequences, concatenate them, reduce with gblocks and finally calculate tree
			tree_obj.calc_tree

			if tree_obj.f_tree then
				session[:tree] = tree_obj.f_tree
			end
			if tree_obj.prots then
				@selected_prots = tree_obj.prots
			end

			if tree_obj.err_msg.blank? then

				# everything worked, move tree file to location where it is accessible to lucullus
				# TODO do not do this if lucullus does not work!
				f_lucullus = File.join(Dir::tmpdir, "cymobase_alignment_" + @id + ".fasta")
				Helper.move_or_copy_file(tree_obj.f_tree, f_lucullus, "0444", "copy")

			else
				@error = tree_obj.err_msg
			end

		else
			@error = "Cannot calculate tree. Please restart gene prediction with MAFFT selected."
		end

	rescue NoMethodError, TypeError, NameError, RuntimeError
		@error = "Cannot calculate tree."

	ensure
		render :show_tree, formats: [:js]
	end

	def download
		file = Tempfile.new(['tree', '.phb'], 'tmp') # locate file in tmp dir of this rails app
		# FileUtils.cp session[:tree], file
		Helper.move_or_copy_file(session[:tree], file,"0444","copy")
		send_file file, :x_sendfile=>true

	rescue RuntimeError => exc
		# sufficient to initialize @error for template :show_tree, because if this is not empty, only error message will be shown
		@error = "Cannot prepare file for download."
		render :show_tree, formats: [:js]
	end

	def gblocks(file)
		stdout = IO.popen([GBLOCKS, file, "-b5=h", "-p=n"])
		output = stdout.read
		stdout.close
		# return value is always false, so parse output to find out if it was successful
		if output.include?("selected block(s)") then
			return true
		else
			return false
		end
	end

	def fasttree(file_in, file_out)
		is_success = system(FastTree, "-out", file_out, file_in)
		return is_success
	end

	def concat_seqs(concat_seqs, max_length, file)
		# concatenate sequences from files for every species
		fh = File.new(file, "w")

		concat_seqs.each do |abbr, seqs|
			conatenated = ""

			# double check that for every protein a sequence/ placeholder was collected
			if seqs.size != params[:n_prot].to_i then
				@error = "Internal error. Concatenate sequences failed."
				return false
			end

			# actual concatenation
			seqs.each_with_index do |seq, ind|
				if seq.length < max_length[ind] then
					# ensure lenght, fill up with gaps if neccessary
					seq = seq + "-" * (max_length[ind] - seq.length)
				end
				# concatenate
				conatenated += seq
			end

			# save concatenated sequences
			fh.puts PredictionsController.new.str2fasta(abbr, conatenated, true) # true: no linebreaks after 80 chars
		end
		fh.close
		return true
	end

	def filename2protname(file)
		basename = file =~ /.+\/([-a-z0-9]+)-[0-9]+-aligned.fasta$/ ? $1 : file
		prot_name = basename.gsub("-", " ").titleize
	end

	def sp_abbr2sp_name
		ref_data, fatal_error = load_ref_data
		list = {}
		if fatal_error.empty? then
			ref_data.keys.each do |prot|
				ref_data[prot]["genes"].keys.each do |k|
					org = ref_data[prot]["genes"][k]["species"]
					if ! list.has_key?(org) then
						# puts k
						list[org] = k.match(/([A-Z][_a-z]+)[A-Z]\w*/)[1]
					end # ! list.has_key?(org)
				end # ref_data[prot]["genes"].each
			end # ref_data.each
		end # if fatal_error.empty?
		return list
	end


# 	def calc_tree_alt

# 		@selected_prots = {}
# 		@error = ""
# 		@minor_error = []
# 		@id = ""

# 		f_out = File.join(BASE_PATH, session[:file][:id], "concat.fasta")
# 		concat_seqs = {}
# 		max_length = []


# 		# if mafft was not used to create the alignment, files will contain only predicted seq and one reference seq
# 		# => useless for tree calucation
# 		if session[:align][:algo] == "mafft" then

# 			# if predict_more was used, there will be more than one alignment for same class/family
# 			# -> use always first prediction (-1-aligned.fasta)
# 			file_list = Dir.glob("/tmp/cug/#{session[:file][:id]}/*-1-aligned.fasta")
# 			# random collection of files
# 			selected_files = file_list.sample(params[:n_prot].to_i)

# 			selected_files.each_with_index do |file, i_file|

# 				headers, seqs = PredictionsController.new.fasta2str(File.read(file))

# 				sp_abbrs = headers.collect{ |he| he.match(/
# 					(Prediction|[A-Z][_a-z]+) # match Prediction OR species abbreviation -- everything until the 2. uppercase char
# 					([A-Z][\w]+|) # match protein class and variant -- everything after 2. uppercase char OR nothing -- if Prediction was matched before
# 					/x)[1] }

# 				# find new species abbreviations and prepare concat_seqs entry
# 				new_sp_abbrs = sp_abbrs.uniq - concat_seqs.keys
# 				new_sp_abbrs.each do |abbr|
# 					# puts "new sp_abbr: #{abbr}"
# 					concat_seqs[abbr] = []
# 					i_file.times { concat_seqs[abbr].push("") }
# 					# puts concat_seqs.values.first.length
# 					# puts concat_seqs[abbr].length
# 				end

# 				# collect sequences
# 				concat_seqs.each do |sp, seq_arr|
# 					ind = headers.index { |he| he =~ /#{sp}([A-Z]\w+|)$/ }
# 					if ind.nil? then 
# 						# protein not encoded in this species
# 						# add gap 
# 						concat_seqs[sp].push("")
# 						# puts "Not found: #{sp} in #{headers.sort.join(", ")}"
# 					else
# 						# protein is encoded
# 						# add sequence
# 						concat_seqs[sp].push(seqs[ind])
# 						# puts "#{sp} => #{headers[ind]}"
# 					end
# 				end
				
# 				# store maximal sequence length for every file
# 				max_length.push seqs.collect {|s| s.length}.max

# 				# add prot name + number of sequences to info - hash
# 				prot_name = filename2protname(file)
# 				@selected_prots[prot_name] = sp_abbrs.uniq.size
# 			end # selected_files.each_with_index

# 			# concatenate sequences
# 			is_success = concat_seqs(concat_seqs, max_length, f_out)
# 			if is_success then
# 				# reduce alignment with gblocks
# 				is_success = gblocks(f_out)
# 				if is_success then
# 					# calculate tree
# 					f_gb = f_out + "-gb"
# 					# TODO clean up if lucullus does not work! use /tmp/cug/session/concat.phb then instead
# 					@id = "cugtree" + rand(1000000000).to_s
# 					f_tree = Dir::tmpdir + "/cymobase_alignment_" + @id + ".fasta"
# 					session[:tree] = f_tree
# puts "---"
# puts f_tree
# 					is_success = fasttree(f_gb, f_tree)
# 					if ! is_success then
# 						@error = "FastTree failed."
# 					end # is_success fasttree
# 				else
# 					@error = "Gblocks failed."
# 				end # is_success gblocks
# 			else
# 				@error = "Sequence concatenation failed."
# 			end # is_success concat_seqs

# 		else
# 			@error = "Cannot calculate tree. Please restart gene prediction with MAFFT selected."
# 		end	
# 		render :show_tree, formats: [:js]
# 	end


end
