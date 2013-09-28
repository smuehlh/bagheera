class Tree

	attr_reader :file_basename, :f_tree, :f_fasta,
		:files, :prots, :err_msg

	def initialize(n_prots, file_basename)
		@file_basename = file_basename
		@files = collect_proteins(n_prots)
		@err_msg = ""
	end

	def collect_proteins(x)
		file_list = Dir.glob( File.join( @file_basename, "*-1-aligned.fasta") )
		# use every file if number of requested files exceeds actual number of files
		if file_list.size < x then
			x = file_list.size
		end
		# random collection of files
		file_list.sample(x)
	end

	def calc_tree
		@f_fasta = File.join(@file_basename, "concat.fasta") 
		@f_tree = File.join(@file_basename, "tree.tre")

		# collect sequences

		@err_msg = catch(:problem) {

			seqs, max_length = collect_seqs

			# concatenate sequences
			concat_seqs(seqs, max_length)

			# reduce alignment with gblocks
			run_gblocks

			# calculate tree
			f_gb = @f_fasta + "-gb"

			run_fasttree(f_gb)

			"" # default if no error was thrown
		}
	end


	def collect_seqs
		collected_seqs = {}
		max_length = []
		@prots = {}

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
			@prots[protname] = sp_abbrs.uniq.size

		end # @files.each_with_index
		return collected_seqs, max_length
	end

	def concat_seqs(collected_seqs, max_length)
		# concatenate sequences from files for every species
		fh = File.new(@f_fasta, "w")

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
		rescue => exc
			Helper.worked_or_throw_error(false, "Sequence concatenation failed.")
		end

	def filename2protname(file)
		basename = file =~ /.+\/([-a-z0-9]+)-[0-9]+-aligned.fasta$/ ? $1 : file
		prot_name = basename.gsub("-", " ").capitalize
	end

	def run_gblocks
		is_success = ProgCall.gblocks(@f_fasta)
		Helper.worked_or_throw_error(is_success, "Gblocks failed.")
	end

	def run_fasttree(f_gb)
		is_success = true
		Helper.file_exist_or_die(f_gb)
		is_success = ProgCall.fasttree(@f_tree, f_gb)

	rescue RuntimeError
		# file f_gb does not exist
		is_success = false
	ensure
		Helper.worked_or_throw_error(is_success, "FastTree failed.")
	end

end