class TreeController < ApplicationController


# TODO

# 2) vllt proteine ausschliessen, von denen es zu wenig referenz-seqs gibt
# 3) abbr durch speziesnamen ersetzen in concat.fasta 
# -> nur falls gblocks, fasttree und lucullus mit den leerzeichen klarkommen
# 4) capture variants?
# 5) nimm # dateien falls #datein < 5 oder < 10

	def calc_tree

		@selected_prots = {}
		@error = ""
		@minor_error = []
		@id = ""

		f_out = "#{BASE_PATH}#{session[:file][:id]}/concat.fasta"
		concat_seqs = {}
		max_length = []


		# if mafft was not used to create the alignment, files will contain only predicted seq and one reference seq
		# => useless for tree calucation
		if session[:align][:algo] == "mafft" then

			# if predict_more was used, there will be more than one alignment for same class/family
			# -> use always first prediction (-1-aligned.fasta)
			file_list = Dir.glob("/tmp/cug/#{session[:file][:id]}/*-1-aligned.fasta")
			# random collection of files
			selected_files = file_list.sample(params[:n_prot].to_i)

			selected_files.each_with_index do |file, i_file|

				headers, seqs = PredictionsController.new.fasta2str(File.read(file))

				sp_abbrs = headers.collect{ |he| he.match(/
					(Prediction|[A-Z][_a-z]+) # match Prediction OR species abbreviation -- everything until the 2. uppercase char
					([A-Z][\w]+|) # match protein class and variant -- everything after 2. uppercase char OR nothing -- if Prediction was matched before
					/x)[1] }

				# find new species abbreviations and prepare concat_seqs entry
				new_sp_abbrs = sp_abbrs.uniq - concat_seqs.keys
				new_sp_abbrs.each do |abbr|
					# puts "new sp_abbr: #{abbr}"
					concat_seqs[abbr] = []
					i_file.times { concat_seqs[abbr].push("") }
					# puts concat_seqs.values.first.length
					# puts concat_seqs[abbr].length
				end

				# collect sequences
				concat_seqs.each do |sp, seq_arr|
					ind = headers.index { |he| he =~ /#{sp}([A-Z]\w+|)$/ }
					if ind.nil? then 
						# protein not encoded in this species
						# add gap 
						concat_seqs[sp].push("")
						# puts "Not found: #{sp} in #{headers.sort.join(", ")}"
					else
						# protein is encoded
						# add sequence
						concat_seqs[sp].push(seqs[ind])
						# puts "#{sp} => #{headers[ind]}"
					end
				end
				
				# store maximal sequence length for every file
				max_length.push seqs.collect {|s| s.length}.max

				# add prot name + number of sequences to info - hash
				prot_name = filename2protname(file)
				@selected_prots[prot_name] = sp_abbrs.uniq.size
			end # selected_files.each_with_index

			# concatenate sequences
			is_success = concat_seqs(concat_seqs, max_length, f_out)
			if is_success then
				# reduce alignment with gblocks
				is_success = gblocks(f_out)
				if is_success then
					# calculate tree
					f_gb = f_out + "-gb"
					# TODO clean up if lucullus does not work! use /tmp/cug/session/concat.phb then instead
					@id = "cugtree" + rand(1000000000).to_s
					f_tree = Dir::tmpdir + "/cymobase_alignment_" + @id + ".fasta"
					session[:tree] = f_tree
puts "---"
puts f_tree
					is_success = fasttree(f_gb, f_tree)
					if ! is_success then
						@error = "FastTree failed."
					end # is_success fasttree
				else
					@error = "Gblocks failed."
				end # is_success gblocks
			else
				@error = "Sequence concatenation failed."
			end # is_success concat_seqs

		else
			@error = "Cannot calculate tree. Please restart gene prediction with MAFFT selected."
		end	
		render :show_tree, formats: [:js]
	end

	def download
		file = Tempfile.new(['tree', '.phb'], 'tmp') # locate file in tmp dir of this rails app
		FileUtils.cp session[:tree], file
		send_file file, :x_sendfile=>true
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


# 		def calc_tree_alt
# 		@selected_prots = []
# 		@error = ""
# 		@minor_error = []

# 		f_out = "#{BASE_PATH}#{session[:file][:id]}/concat.fasta"
# TODO
# # 1) concat funktioniert nicht: meistens nicht fuer jedes prot eine sequenz!
# # -> doch die concat methode aus concat_neu verwenden?
# # 2) vllt proteine ausschliessen, von denen es zu wenig referenz-seqs gibt
# # -> auf jeden fall die anzahl ref sequenzen in view auflisten!
# # 3) abbr durch speziesnamen ersetzen in concat.fasta 
# # -> nur falls gblocks, fasttree und lucullus mit den leerzeichen klarkommen

# 		# if mafft was not used to create the alignment, files will contain only predicted seq and one reference seq
# 		# => useless for tree calucation
# 		if session[:align][:algo] == "mafft" then
# 			# if predict_more was used, there will be more than one alignment for same class/family
# 			# -> use always first prediction (-1-aligned.fasta)
# 			file_list = Dir.glob("/tmp/cug/#{session[:file][:id]}/*-1-aligned.fasta")
# 			@selected_prots = file_list.sample(params[:n_prot].to_i)

# 			# concatenate sequences
# 			# 1) get and order all sequences for all species
# 			all_sp_abbrs = []
# 			seqs_unstruct = {}
# 			max_length = [] # maximal sequence length for each protein

# 			@selected_prots.each_with_index do |file, i_file|
# puts file
# 				headers, seqs = PredictionsController.new.fasta2str(File.read(file))

# 				sp_abbrs = headers.collect{ |he| he.match(/
# 					(Prediction|[A-Z][_a-z]+) # match Prediction OR species abbreviation -- everything until the 2. uppercase char
# 					([A-Z][\w]+|) # match protein class and variant -- everything after 2. uppercase char OR nothing -- if Prediction was matched before
# 					/x)[1] }
# 				all_sp_abbrs |= sp_abbrs

# 				max_length << seqs.collect {|s| s.length}.max

# 				# iterate over all species found so far to ensure missing seqs in fasta can be handeled
# debugger
# 				all_sp_abbrs.each do |abbr|
# 					# match abbr to headers to find correct sequence

# 					i_seq = headers.index{ |he| he =~ /#{abbr}([A-Z]\w+|)$/ }

# 					if i_seq.nil? then
# puts "not found #{abbr} in #{headers.sort.join(", ")}"
# 						prot_name = filename2protname(file)
# 						@minor_error |= ["Error while processing #{prot_name}"]
# 						next
# 					end
# puts i_seq
# 					# add sequence to "results" hash
# 					if ! seqs_unstruct.has_key?(abbr) then
# puts "Case 1: add key #{abbr}"
# 						# (this abbr was not seen before, but it is seen this time)
# 						# first, prepare "results" and fill up with placeholder(s)
# 						seqs_unstruct[abbr] = []
# 						i_file.times { seqs_unstruct[abbr] << "" }
# 						seqs_unstruct[abbr] << seqs[i_seq]
# 					end
# puts "Case 2: abbr was seen before and this time"
# # puts "#{abbr} => #{headers[i_seq]}"

# puts seqs_unstruct[abbr].size
# 					seqs_unstruct[abbr] << seqs[i_seq]	

# 					# end # if new_keys.include?(abbr)
# 				end # all_sp_abbrs.each
# 			end # @selected_prots.each_with_index
# debugger
# 			# concatenate sequences from files for every species
# 			is_success = concat_seqs(seqs_unstruct, max_length, f_out)
# 			if is_success then
# 				# reduce alignment with gblocks
# 				is_success = gblocks(file)
# 			end # is_success concat_seqs

# 		else
# 			@error = "Cannot calculate tree. Please restart gene prediction with MAFFT selected."
# 		end	
# 		render :show_tree, formats: [:js]
# 	end

end
