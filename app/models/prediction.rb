HYDORPHOBIC_AAS = ["V", "I", "L", "M", "F"] unless defined?(HYDORPHOBIC_AAS)
POLAR_AAS = ["S", "T", "A", "C"] unless defined?(POLAR_AAS)


class Prediction
	attr_reader :pred_cug, :pred_seq, :pred_dnaseq, :this_hit_nr, :n_blast_hits, :err_msg, :used_prfl
	
	@@align_method = "mafft"
	@@align_config = "tttt"
	
	def initialize(prot_obj, hit_nr=1, file_basename)
		@protein = prot_obj
		@this_hit_nr = hit_nr
		@file_basename = file_basename
		@err_msg = ""
		@used_prfl = true # over-enthusiastic ...
	end

	# get number of blast hits before predict() was called
	def sneak_n_blast_hits
		if @this_hit_nr > 1 then
			f_blasthits = File.join(@file_basename, "#{@protein.prot_basename}.blast")
			Helper.file_exist_or_die(f_blasthits)
			hits = File.read(f_blasthits)
			n_blast_hits = hits.lines.to_a.size
		else
			Helper.raise_runtime_error
		end
		return n_blast_hits
	end

	def predict

		@err_msg = catch(:problem) {
			if @this_hit_nr > 1 then
				# blast hits already calculated
				@f_blasthits = File.join(@file_basename, "#{@protein.prot_basename}.blast")
			else
				# 1) run BLAST
				run_blast
			end

			# 2) run augustus
			seq_id, start, stop, strand = parse_blast_hits

			start, stop = extend_blast_hit(start, stop)

			output = run_fastacmd(seq_id, start, stop, strand)
			# hit extended beyond its lenght?
			if output.match(/From location cannot be greater than/) then
				stop = 0 # this will cause fastacmd to return everything from start till end of contig
				run_fastacmd(seq_id, start, stop, strand)
			end

			run_augustus

			parse_augustus

			align

			# 3) map predicted CUG codons
			res = compare_pred_ref

			@pred_cug = res
		}
		if @pred_cug then
			@err_msg = ""
		end

	end


	def run_blast
		@f_blasthits = File.join(@file_basename, "#{@protein.prot_basename}.blast")
		is_success = ProgCall.blast(@f_blasthits, @protein.ref_file)
		Helper.worked_or_throw_error(is_success, "BLAST failed.")
	end

	def run_fastacmd(seq_id, start, stop, strand)
		@tmp_file = File.join(@file_basename, Helper.get_tmp_file)
		is_success, output = ProgCall.fastacmd(@tmp_file, seq_id, start, stop, strand)
		Helper.worked_or_throw_error(is_success, "Gene prediction (fastacmd) failed.")
		return output
	end

	def run_augustus
		@f_augustus = File.join(@file_basename, Helper.get_tmp_file("augustus"))
		is_success = true # being over-enthusiastic ...
		# test if profile file exists
		if Helper.does_file_exist(@protein.ref_prfl_file) then
			is_success, exit_status = ProgCall.augustus(@f_augustus, @tmp_file, @protein.ref_prfl_file)

			# test if there was a memory error while using protein profile
			if exit_status == 6 then
				# try again without protein profile
				@used_prfl = false
				is_success, exit_status = ProgCall.augustus(@f_augustus, @tmp_file)
			end

		else
			@used_prfl = false
			is_success, exit_status = ProgCall.augustus(@f_augustus, @tmp_file)
		end
		
		Helper.worked_or_throw_error(is_success, "Gene prediction failed.")
	end

	def run_mafft(f_in_pred)
		is_success = ProgCall.mafft(@f_align, f_in_pred, @protein.ref_alignment_file)
		Helper.worked_or_throw_error(is_success, "Aligning sequences (MAFFT) failed.")
	end

	def run_pairalign(f_in)
		is_success = ProgCall.pairalign(@f_align, f_in, @@align_config, @@align_method)
		Helper.worked_or_throw_error(is_success, "Aligning sequences (PairAlign) failed.")
	end

	# parse blast output 
	# combines hit of interest with the following (!) blast hits, if they overlap
	def parse_blast_hits
		hits = File.read(@f_blasthits)
		hit = hits.lines.to_a[@this_hit_nr-1]
		@n_blast_hits = hits.lines.to_a.size
		if hit.nil? then
			Helper.worked_or_throw_error(false, "Requested hit #{@this_hit_nr} execceds number of BLAST hits.")
		end

		fields = hit.split("\t") 
		seq_id = fields[1] # sequence contig
		start = fields[8].to_i # sequence start
		stop = fields[9].to_i # sequence stop
		strand = "plus" # set strand to plus for fastacmd
		if stop < start then
			# change start/ stop if prediction on minus strand 
			strand = "minus" # set strand to minus for fastacmd
			tmp = start
			start = stop
			stop = tmp
		end

		# test if another hit on same contig exits, which overlaps

		indices = hits.lines.to_a.each_with_index.map{ |ele, ind| (ele.include?(seq_id)) ? ind : nil }.compact
		indices -= ( 0..(@this_hit_nr-1) ).to_a

		indices.each do |ind|
			this_hit = hits.lines.to_a[ind]
			this_parts = this_hit.split("\t")
			this_start = this_parts[8].to_i
			this_stop = this_parts[9].to_i
			this_strand = "plus"
			if this_stop < this_start then
				# change start/ stop if prediction on minus strand 
				this_strand = "minus" # set strand to minus for fastacmd
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
		return seq_id, start, stop, strand
	end

	def extend_blast_hit(start, stop)
		if start-500 < 0 then
			# check if start is valid
			start = 0
		else
			start = start-500
		end
		return start, stop
	end

	def parse_augustus
		output = File.read(@f_augustus)
		@pred_seq = output.match(/protein sequence = \[(.*)\]/m)[1].gsub("\n# ", "").upcase
		@pred_dnaseq = output.match(/coding sequence = \[(.*?)\]/m)[1].gsub("\n# ", "").upcase

		if @pred_seq.blank? || @pred_dnaseq.blank? then
			puts "parse_augustus -> something blank?"
			Helper.worked_or_throw_error(false, "AUGUSTUS failed.")
		end
	end

	# align reference protein with the predicted protein
	def align
		# align sequences
		@f_align = File.join(@file_basename, "#{@protein.prot_basename}-#{@this_hit_nr}-aligned.fasta")
		f_in = File.join(@file_basename, Helper.get_tmp_file)
		File.open(f_in, 'w') {|f| f.write( Helper::Sequence.str2fasta("Prediction", @pred_seq, true) )}
		case @@align_method
		when "mafft"
			# f_in contains only predicted seq
			run_mafft(f_in)
		else
			# add ref seq to f_in: it must contain both seqs
			File.open(f_in, 'a') {|f| f.write("\n" + Helper::Sequence.str2fasta(@protein.ref_key, @protein.ref_seq, true) )}
			run_pairalign(f_in)
		end

		# collect aligned predicted and ref sequences
		headers, seqs = Helper::Sequence.fasta2str(File.read(@f_align))
		ind = headers.find_index{|str| str.include?(@protein.ref_key)}
		if ind then
			@refseq_aligned = seqs[ind]
		else
			Helper.worked_or_throw_error(false, "Aligning sequences failed.")
		end
		ind = headers.find_index{|str| str.include?("Prediction")}
		if ind then
			@predseq_aligned = seqs[ind]
		else
			Helper.worked_or_throw_error(false, "Aligning sequences failed.")
		end
	end

	def compare_pred_ref
		results = {ctg_pos: [], ref_chem: {}, ref_ctg: {}, pred_ctg_unaligned: [], no_ref_ctg: []}
		headers, seqs = Helper::Sequence.fasta2str(File.read(@f_align))
		pred_seq_aligned = seqs[-1]
		codons = pred_dnaseq.scan(/.{1,3}/)
		ctg_pos = codons.each_with_index.map{ |ele, ind| (ele == "CTG") ? ind : nil }.compact
		results[:ctg_pos] = ctg_pos
		results[:pred_prot] = pred_seq_aligned.gsub("-", "")
		# has predicted protein CTG codon(s)?
		if ctg_pos.empty? then
			Helper.worked_or_throw_error(false, "No CTG in predicted gene")
		end
		results[:ref_seq_num] = @protein.ref_alignment.keys.size # total number of reference sequences

		ctg_pos_mapped, ref_alignment_cols, ref_codons = map_ctg_pos(ctg_pos)

		if ctg_pos_mapped.compact.any? then

			# 1) preference for hydrophob/ polar residues?
			ref_alignment_cols.each_with_index do |col, ind|

				# handle unmatched CTG positions
				if ctg_pos_mapped[ind].nil? then 
					results[:pred_ctg_unaligned].push ctg_pos[ind]
					next
				end
				aa_freq, aa_num = word_frequency(col)
				is_significant, prob_transl = predict_translation(aa_freq)
				results[:ref_chem][ctg_pos[ind]] = {aa_comp: aa_freq, aa_num: aa_num}
				if is_significant then
					results[:ref_chem][ctg_pos[ind]][:is_significant] = true
					results[:ref_chem][ctg_pos[ind]][:transl] = prob_transl
				else
					results[:ref_chem][ctg_pos[ind]][:is_significant] = false
				end
			end

			# 2) CTGs in reference data at CTG positions in predicted sequence?
			results[:ref_ctg], results[:no_ref_ctg] = ref_ctg_usage(ref_codons, ref_alignment_cols, ctg_pos)
			
		else 
			Helper.worked_or_throw_error(false, "No CTG position in predicted gene aligned with reference sequences")
		end
		return results
	end

	def map_ctg_pos(ctg_pos)
		ref_pos = []
		ref_cols = []
		ref_codons = []

		ctg_pos.each do |pos|

			pred_apos = Helper::Sequence.sequence_pos2alignment_pos(pos, @predseq_aligned)
			if @refseq_aligned[pred_apos] == "-" then
				ref_pos << nil
				ref_cols << []
				ref_codons << []
				next
			end

			ref_spos = Helper::Sequence.alignment_pos2sequence_pos(pred_apos, @refseq_aligned)
			if ref_spos.nil? || @protein.ref_alignment[@protein.ref_key].nil? then
				ref_pos << nil
				ref_cols << []
				ref_codons << []
				next
			end
			ref_cymopos = Helper::Sequence.sequence_pos2alignment_pos(ref_spos, @protein.ref_alignment[@protein.ref_key])

			col = []
			codons = []
			@protein.ref_alignment.each do |cymo_header, cymo_prot|
				
				if cymo_prot[ref_cymopos] == "-" then
					col << nil
					codons << nil
				else
					# save codon and translation
					# get codon translation from alignment to get translation regardless gene status (also if incomplete)
					
					spos = Helper::Sequence.alignment_pos2sequence_pos(ref_cymopos, cymo_prot) # position in sequece without gaps
					this_gene = nil
					if @protein.ref_genes[cymo_header] then
						# gene does exist (so it is not a partial)
						this_gene = @protein.ref_genes[cymo_header]["gene"]
					end

					this_aa = nil
					this_codon = nil
					begin
						# this_aa = Helper::Sequence.extract_translation(this_gene)[spos]
						this_aa = cymo_prot[ref_cymopos]
						this_codon = Helper::Sequence.split_cdna_into_codons(this_gene, spos)
					rescue
						# add 1 amino acid but 3 codons!
						this_aa = "X" if ! this_aa
						this_codon = "XXX" if ! this_codon
					ensure
						col << this_aa
						codons << this_codon
					end
				end

			end

			ref_pos << ref_cymopos
			ref_cols << col.compact
			ref_codons << codons.compact

		end

		return ref_pos, ref_cols, ref_codons
	end

	def word_frequency(arr)
		res = Hash.new(0)
		arr.each { |a| res[a] += 1 }
		res.delete(nil)
		sum = res.inject(0) {|s, (_,val)| s + val}.to_f
		sum = sum - res["-"]
		norm_res = res.each {|k,v| res[k] = v/sum}
		return norm_res, sum.to_i
	end

	# determine most probable translation scheme based on amino acid frequency
	# determination algorithm:
	# max of hydrophobic/ polar/ other amino acids
	# -> standard/ alternative/ indiscriminative position
	def predict_translation(aa_freq)
		aa_freq.delete("-")
		if aa_freq.empty? then
			return false
		end

		num_aas = aa_freq.inject(0) {|sum, (_,val)| sum + val} # total number of amino acids
		pol_aas = POLAR_AAS.collect{|aa| aa_freq[aa]}.compact.sum # polar amino acids
		hyd_aas = HYDORPHOBIC_AAS.collect{|aa| aa_freq[aa]}.compact.sum # hydrophobic amino acids
		pct_pol = pol_aas/num_aas.to_f
		pct_hyd = hyd_aas/num_aas.to_f
		pct_other = (num_aas - pol_aas - hyd_aas )/ num_aas.to_f

		# maximum of values gives predicted translation
		# requirements for discriminative position: majority of amino acids is hydrophobic or polar
		# requirements of leucine position: majority of amino acids is hydrophobic
		# requirements of serine position: majority of amino acids is polar
		# ... indiscriminative: majority of amino acids neither hydrophobic nor polar
		max_val = [pct_hyd, pct_pol, pct_other].max
		if max_val == pct_pol then 
			is_discrim = true
			transl = "S"
		end
		if max_val == pct_hyd then 
			is_discrim = true
			transl = "L"
		end
		if max_val == pct_other then 
			is_discrim = false
			transl = ""
		end

		return is_discrim, transl
	end

	def ref_ctg_usage(codons, aas, ctg_pos)
		counts = {} # usage of CTG codons in reference data, per CTG position in prediction
		no_ctgs = [] # predicted CTG positiions with no CTG codons in reference data

		codons.each_with_index do |codons_thiscol, ind|
			indices = codons_thiscol.find_each_index("CTG")
			ctg_num = indices.size
			if ctg_num > 0 then

				counts[ctg_pos[ind]] = {ctg_usage: Hash.new(), ctg_num: ctg_num} 
				counts[ctg_pos[ind]][:ctg_usage]["S"] = 0
				counts[ctg_pos[ind]][:ctg_usage]["L"] = 0

				ctg_translations = indices.collect {|x| aas[ind][x]}
				n_ser = ctg_translations.count("S")
				n_leu = ctg_translations.count("L")

				pct_ser = n_ser / ctg_num.to_f
				pct_leu = n_leu / ctg_num.to_f

				counts[ctg_pos[ind]][:ctg_usage]["S"] = pct_ser
				counts[ctg_pos[ind]][:ctg_usage]["L"] = pct_leu

				if (pct_ser + pct_leu) == 0.0 then
					counts[ctg_pos[ind]][:transl] = "?"
					next
				end

				counts[ctg_pos[ind]][:is_significant], counts[ctg_pos[ind]][:transl] = predict_translation(counts[ctg_pos[ind]][:ctg_usage])
			else
				no_ctgs.push(ctg_pos[ind])
			end # if ctg_num
		end # codons.each_with_index
		return counts, no_ctgs
	end


	def save(results)
		# "fail-save" save only if new values are not nil 
		results[:ref_species] = @protein.ref_species if @protein.ref_species
		results[:ref_prot] = @protein.ref_seq if @protein.ref_seq
		results[:ref_seq_num] = @protein.ref_alignment.keys.size if @protein.ref_alignment.keys
		results[:n_hits] = @n_blast_hits if @n_blast_hits
		results[:hit_shown] = @this_hit_nr if @this_hit_nr
		results[:pred_prot] = @pred_seq if @pred_seq

		if @pred_cug then
			results.merge!(@pred_cug)
		end
	end

	def save_message(results)
		results[:message] |= [@err_msg] if ! @err_msg.blank?
	end

end