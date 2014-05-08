class ProteinFamily
	attr_reader :prot, :prot_basename, 
		:ref_alignment_file, :ref_alignment, :ref_genes, :ref_prfl_file

	@@ref_data_path = BASE_PATH_PROTEIN # overwrite if reference data without Ca_b are used

	def initialize(prot, data)
		@prot = prot 
		@prot_basename = get_filesave_name 
		@ref_prfl_file = get_fname("prfl")
		@ref_alignment_file = get_fname("refalignment")
		@ref_alignment, @ref_genes = get_refdata(data)
	end

	def get_filesave_name
		@prot.gsub(" ", "-").downcase
	end

	def get_refdata(data)
		seqs = Hash[*data["alignment"].split("\n")]
		seqs.keys.each {|k| seqs[ k.sub(">", "") ] = seqs.delete(k)}
		genes = filter_genes(data["genes"])
		return seqs, genes
	end

	def filter_genes(all_genes)
		# remove all incomple gene structures
		all_genes.delete_if { |_,g| g["completeness"] != "complete" ||
			(YAML::load(g["gene"])[0]["status"] && YAML::load(g["gene"])[0]["status"] == "incomplete") }
		return all_genes
	end

	def get_fname(type)
		file = ""
		case type
		when "prfl"
			file = File.join(@@ref_data_path, "#{@prot_basename}.prfl")
		when "refalignment"
			file = File.join(@@ref_data_path, "#{@prot_basename}.fasta")
		end
		Helper.file_exist_or_die(file)
		return file
	end

end