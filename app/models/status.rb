class Status

	def initialize()
	end

	# update data statistics for this query 
	def self.update(pred_data, stats)

		if pred_data[:pred_prot] && ! pred_data[:pred_prot].blank? then
			# number of proteins with predicted protein sequence
			stats[:n_prots] += 1
		end

		if pred_data[:pred_prot] && pred_data[:ctg_pos] && pred_data[:ctg_pos].empty? then 
			# number of proteins without any ctg
			stats[:n_prots_wo_ctg] += 1
		end

		if pred_data[:ref_species] then
			# one more protein was tested
			stats[:total_n_prots] += 1
		end

		if pred_data[:ctg_pos] then
			# number of CTGs in predicted prot
			stats[:pred_ctg] += pred_data[:ctg_pos].size
		end

		if pred_data[:pred_ctg_unaligned] then
			stats[:pred_ctg_unaligned] += pred_data[:pred_ctg_unaligned].size
		end

		if pred_data[:ref_chem] then 
			# number of ctg codons discriminative (alternative/ standard codon usage) and not discriminative
			n_ser, n_leu, n_unknown = count_unknown_and_suggested_transl( pred_data[:ref_chem] )
			stats[:ref_chem_leu] += n_leu
			stats[:ref_chem_ser] += n_ser
			stats[:ref_chem_unknown] += n_unknown
		end

		if pred_data[:ref_ctg] then 
			# number of ctg codons discriminative (alternative/ standard codon usage) and not discriminative
			n_ser, n_leu, n_unknown = count_unknown_and_suggested_transl( pred_data[:ref_ctg] )
			stats[:ref_ctg_leu] += n_leu
			stats[:ref_ctg_ser] += n_ser
			stats[:ref_ctg_unknown] += n_unknown
		end

		if pred_data[:no_ref_ctg] then 
			stats[:no_ref_ctg] += pred_data[:no_ref_ctg].size
		end

		if pred_data[:ref_chem] && pred_data[:ref_ctg] then 
			# number of ctg codons discriminative for _both_ alternative and standard by chemical prop and reference ctg!
			n_contr = count_contradicting_pos_chem_ctg(pred_data[:ref_chem], pred_data[:ref_ctg])
			stats[:contradiction_ref_chem_ref_ctg] += n_contr
		end

		if pred_data[:score] then 
			stats[:trna_score] = pred_data[:score]
		end
		if pred_data[:blast_hits] then 
			n_ser, n_leu, n_unknown = count_unknown_and_suggested_transl( pred_data[:blast_hits] )
			stats[:trna_leu] += n_leu
			stats[:trna_ser] += n_ser
			stats[:trna_unknown] += n_unknown
		end

		# if pred_data[:ref_chem]&& pred_data[:ref_ctg] then
		# 	# number of CTGs discriminative for alternative usage/ standard usage
		# 	arr = suggested_transl(pred_data[:ref_chem], pred_data[:ref_ctg])
		# 	stats[:ser]  += arr[0] # predicted translation is serine (alternative usage)
		# 	stats[:leu] += arr[1] # predicted translation is leucine (standard usage)
		# end
		# if pred_data[:ref_ctg] then
		# 	n_conserved_pos = pred_data[:ref_ctg].keys.size
		# 	# number of CTG positions mapped to a CTG position in reference data -> is already conserved!
		# 	stats[:conserved_pos] += n_conserved_pos

		# 	# number of proteins which contain conserved CTG positions
		# 	stats[:prot_conserved_pos] += 1 if n_conserved_pos > 0
		# end

		return stats
	end

	# read stats from file, updates them and writes them back to file
	def self.update_file(pred_data, file)
		if Helper.does_file_exist(file) then 
			stats = read(file)
		else
			stats = Hash.new(0)
		end
		stats = update(pred_data, stats)
		save(file, stats)
	end

	# save statistics to file
	def self.save(file, stats)
		str = stats.map {|obj| obj.join(":") }.join("\n")
		File.open(file, 'w') { |f| f.write(str) }
	end

	# read statistics from file into Hash
	def self.read(file)
		str = File.read(file)
		stats = Hash.new(0)

		# parse data and put them into the hash
		str.scan(/([\w_]+):(\d+)/).each do |key, val|
			key = key.to_sym # convert all keys to symbols
			val = val.to_i # convert all numbers to integers
			stats[key] = val
		end
		return stats
	end

	# returns number of positions leading to suggest serine or leucine or no suggestion 
	def self.count_unknown_and_suggested_transl(ref_dat)
		ser_pos, leu_pos = ctg_pos_by_suggested_transl(ref_dat)
		unknown_pos = ctg_pos_by_unknown_transl(ref_dat)
		return ser_pos.size, leu_pos.size, unknown_pos.size
	end

	# parse ref_chem and ref_ctg for ctg positions leading to alternative, standard, or unknown
	def self.ctg_pos_by_suggested_transl(ref_dat)
		pos_ser = ref_dat.collect {|k,v| k if v[:transl] == "S" && v[:is_significant] }
		pos_leu = ref_dat.collect {|k,v| k if v[:transl] == "L" && v[:is_significant] }
		return pos_ser.compact, pos_leu.compact
	end

	# count number of ctg positions without an clear assignment
	# all ctg positions which suggest neither leucine nor serine translation
	def self.ctg_pos_by_unknown_transl(ref_dat)
		pos_indiscrim = ref_dat.collect {|k,v| k if ! v[:is_significant] }
		return pos_indiscrim.compact
	end
	# def self.ctg_pos_by_unknown_transl(pos_ser, pos_leu, all_ctgpos)
	# 	return all_ctgpos - (pos_ser | pos_leu)
	# end

	# count number of ctg positions where amino acid chemistry and ctg translation lead to conflicting results
	def self.count_contradicting_pos_chem_ctg(ref_dat_chem, ref_dat_ctg)
		ser_pos_chem, leu_pos_chem = ctg_pos_by_suggested_transl(ref_dat_chem)
		ser_pos_ctg, leu_pos_ctg = ctg_pos_by_suggested_transl(ref_dat_ctg)

		return (ser_pos_chem & leu_pos_ctg).size + (leu_pos_chem & ser_pos_ctg).size
	end

	# def self.suggested_transl(ref_chem, ref_ctg)

	# 	# counts for aa distribution in reference data
	# 	pos_ser = ref_chem.collect{|k,v| k if v[:transl] == "S"}
	# 	pos_leu = ref_chem.collect{|k,v| k if v[:transl] == "L"}

	# 	# add counts for ctg usage in reference data
	# 	pos_ser |= pos_ser + ref_ctg.collect{|k,v| k if v[:transl] == "S"} # set union
	# 	pos_leu |= ref_ctg.collect{|k,v| k if v[:transl] == "L"} # set union

	# 	pos_ser.compact!
	# 	pos_leu.compact!

	# 	alt_usage = (pos_ser - pos_leu).size # set difference
	# 	std_usage = (pos_leu - pos_ser).size # set difference
	# 	# strange = (pos_ser & pos_leu).size # set intersect

	# 	return [alt_usage, std_usage]

	# end

end