class Status

	def initialize()
	end

	# update data statistics for this query 
	def self.update(pred_data, stats)

		if pred_data[:pred_prot] then
			stats[:n_prots] += 1
		end

		if pred_data[:ref_chem]&& pred_data[:ref_ctg] then
			# number of CTGs discriminative for alternative usage/ standard usage
			arr = suggested_transl(pred_data[:ref_chem], pred_data[:ref_ctg])
			stats[:ser]  += arr[0] # predicted translation is serine (alternative usage)
			stats[:leu] += arr[1] # predicted translation is leucine (standard usage)
		end

		if pred_data[:ctg_pos] then
			# number of CTGs in predicted prot
			stats[:pred_ctg] += pred_data[:ctg_pos].size
		end

		if pred_data[:ref_ctg] then
			n_conserved_pos = pred_data[:ref_ctg].keys.size
			# number of CTG positions mapped to a CTG position in reference data -> is already conserved!
			stats[:conserved_pos] += n_conserved_pos

			# number of proteins which contain conserved CTG positions
			stats[:prot_conserved_pos] += 1 if n_conserved_pos > 0
		end

		return stats
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
		stats = Hash[str.scan(/([\w_]+):(\d+)/)] # parse data and put them right into the hash
		stats.each{ |key,val| stats[key] = val.to_i } # convert all numbers (string representation) to integers

		return stats
	end

	def self.suggested_transl(ref_chem, ref_ctg)

		# counts for aa distribution in reference data
		pos_ser = ref_chem.collect{|k,v| k if v[:transl] == "S"}
		pos_leu = ref_chem.collect{|k,v| k if v[:transl] == "L"}

		# add counts for ctg usage in reference data
		pos_ser |= pos_ser + ref_ctg.collect{|k,v| k if v[:transl] == "S"} # set union
		pos_leu |= ref_ctg.collect{|k,v| k if v[:transl] == "L"} # set union

		pos_ser.compact!
		pos_leu.compact!

		alt_usage = (pos_ser - pos_leu).size # set difference
		std_usage = (pos_leu - pos_ser).size # set difference
		# strange = (pos_ser & pos_leu).size # set intersect

		return [alt_usage, std_usage]

	end

end