class ProteinFamily
	attr_reader :prot, :prot_filesave_name, 
		:ref_alignment, :ref_genes

	def initialize(prot, data)
		@prot = prot 
		@prot_filesave_name = get_filesave_name 
		# @ref_prfl_file = get_fname("prfl")
		@ref_alignment, @ref_genes = get_refdata(data)
	end

	def get_filesave_name
		@prot.gsub(" ", "-").downcase
	end

	def get_refdata(data)
		begin
			seqs = Hash[*data["alignment"].split("\n")]
		rescue => exc
			Helper.abort("Alignment of protein #{@prot} is corrupt.")
		end
		seqs.keys.each {|k| seqs[ k.sub(">", "") ] = seqs.delete(k)}
		genes = data["genes"]
		return seqs, genes
	end

	def ensure_same_length

		seq_lines = @ref_alignment.values
		len = seq_lines.collect {|line| line.size}

		if ! (len.min == len.max) then
			# add gaps to end of each line shorter than longest line
			max = len.max 
			@ref_alignment.each do |key, val|
				
				if val.size != max then
					@ref_alignment[key] = val + "-" * (max - val.size)
				end

			end

		end

	end

	# method is based on /fab8/db_scripts/alignment/edit.rb
	def remove_common_gaps
  
		seqs = @ref_alignment.values
		max_len = seqs.collect do |s| s.length end.max

		#look for common gaps
		gaps = []
		max_len.times do |pos|
			if (seqs.select do |s| (!s[pos] || s[pos].chr == "-") end.length == seqs.length) then
				gaps << pos
			end
		end
		gaps = gaps.reverse

		#delete gaps
		@ref_alignment.each_pair do |k, v|
			v = v.split("")
			gaps.each do |pos| v.delete_at(pos) end
			v = v.join("")
			@ref_alignment[k] = v
		end
	end

	def ensure_mafft_is_fine(file)
		# test if mafft recognizes alignment by adding first sequence again to alignment
		tmp_file = File.join("/tmp", "test.fasta")

		File.open(tmp_file, 'w') { |f| f.write( Helper::Sequence.str2fasta(">test",@ref_alignment.values[0]) ) }
		is_success = Helper.run_mafft(file, tmp_file)
		Helper.del_file_or_dir(tmp_file)

		# if this is not successful, mafft does not recognize file as aligned
		if ! is_success then
			puts "#{@prot}:: need to calc alignment again"
			# align reference sequences with mafft
			Helper.run_mafft(file,"",true) # true: save output
			# Helper.worked_or_throw_error(is_success, "Mafft failed: #{file}")
		end
	end
end