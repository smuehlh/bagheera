class Tree

	def initialize(n_prots, file_basename)
		@file_basename = file_basename
		@files = collect_proteins(n_prots)
	end

	def collect_proteins(x)
		file_list = Dir.glob( @file_basename + "*-1-aligned.fasta")
		# random collection of files
		file_list.sample(x)
	end

	def calc_tree
		f_fasta = @file_basename + "concat.fasta"
		@f_tree = @file_basename + "tree.tre"

		# collect sequences
		print "Concatenate sequences"
		seqs, max_length = collect_seqs

		# concatenate sequences
		is_success = concat_seqs(seqs, max_length, f_fasta)

		if is_success then
			Helper.done

			# reduce alignment with gblocks
			print "Reduce alignment "
			is_success = ProgCall.gblocks(f_fasta)

			if is_success then
				Helper.done

				# calculate tree
				f_gb = f_fasta + "-gb"

				print "Calculate tree "
				is_success = ProgCall.fasttree(@f_tree, f_gb)
				Helper.done
				puts "Tree is stored in file #{@f_tree}"

				if ! is_success then
					puts "Failed."
				end # is_success fasttree
			else
				print "... failed.\n"
			end # is_success gblocks

		else
			print " ... failed.\n"
		end # is_success concat_seqs
	end


	def collect_seqs
		collected_seqs = {}
		max_length = []
		@prots = []
		@files.each_with_index do |file, i_file|

			headers, seqs = Helper::Sequence.fasta2str(File.read(file))

			sp_abbrs = headers.collect{ |he| he.match(/
				(Prediction|[A-Z][_a-z]+) # match Prediction OR species abbreviation -- everything until the 2. uppercase char
				([A-Z][\w]+|) # match protein class and variant -- everything after 2. uppercase char OR nothing -- if Prediction was matched before
				/x)[1] }

			# find new species abbreviations and prepare collected_seqs entry
			new_sp_abbrs = sp_abbrs.uniq - collected_seqs.keys
			new_sp_abbrs.each do |abbr|
				collected_seqs[abbr] = []
				i_file.times { collected_seqs[abbr].push("") }
			end

			# collect sequences
			collected_seqs.each do |sp, seq_arr|
				ind = headers.index { |he| he =~ /#{sp}([A-Z]\w+|)$/ }
				if ind.nil? then 
					# protein not encoded in this species
					# add gap 
					collected_seqs[sp].push("")
				else
					# protein is encoded
					# add sequence
					collected_seqs[sp].push(seqs[ind])
				end
			end
			
			# store maximal sequence length for every file
			max_length.push seqs.collect {|s| s.length}.max

			# add prot name + number of sequences to info - hash
			protname = filename2protname(file)
			@prots << protname
			puts "#{protname}: #{sp_abbrs.uniq.size} sequences"
		end # @files.each_with_index
		return collected_seqs, max_length
	end

	def concat_seqs(collected_seqs, max_length, file)
		# concatenate sequences from files for every species
		fh = File.new(file, "w")

		collected_seqs.each do |abbr, seqs|
			conatenated = ""

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
			fh.puts Helper::Sequence.str2fasta(abbr, conatenated, true) # true: no linebreaks after 80 chars
		end
		fh.close
		return true
	end

	def filename2protname(file)
		basename = file =~ /.+\/([-a-z0-9]+)-[0-9]+-aligned.fasta$/ ? $1 : file
		prot_name = basename.gsub("-", " ").capitalize
	end

end