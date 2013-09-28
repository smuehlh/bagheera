class RefProtein < ProteinFamily
	attr_reader  :ref_seq, :ref_key, :ref_species, 
		:ref_file

	def initialize(prot, data, output_file)
		super(prot, data)
		@file_basename = output_file
		@ref_key, @ref_seq = get_refprot 
		@ref_file = write_refprot
		@ref_species = get_refspecies 
	end

	def get_refprot
		# use a protein of the species specified as augustus model as reference
		matched_sp_abbr = ""
		case ProgCall.augustus_species
		when 'candida_albicans'
			if ProteinFamily.class_variable_get(:@@ref_data_path).include?("without_ca_b") then
				matched_sp_abbr = "Ca_a"
			else
				matched_sp_abbr = "Ca_b"
			end
		when 'candida_guilliermondii'
			matched_sp_abbr = "Mrg"
		when 'candida_tropicalis'
			matched_sp_abbr = "Ct_a"
		when 'debaryomyces_hansenii'
			matched_sp_abbr = "Deh"
		when 'eremothecium_gossypii'
			matched_sp_abbr = "Erg"
		when 'kluyveromyces_lactis'
			matched_sp_abbr = "Kl"
		when 'lodderomyces_elongisporus'
			matched_sp_abbr = "Loe"
		when 'pichia_stipitis'
			matched_sp_abbr = "Shs"
		when 'saccharomyces_cerevisiae_S288C'
			matched_sp_abbr = "Sc_c"
		when 'saccharomyces_cerevisiae_rm11-1a_1'
			matched_sp_abbr = "Sc_b"
		when 'yarrowia_lipolytica'
			matched_sp_abbr = "Yl"
		end

		# find a gene, if possible with complete gene sequence
		key = ""
		matched_genes = @ref_genes.keys.select do |key| 
			key =~ /#{matched_sp_abbr}/ && @ref_genes[key]["completeness"] =~ /[complete|partial]/
		end
		if matched_genes.size == 1 then
			key = matched_genes[0]
		elsif matched_genes.any? 
			# complete gene sequence if possible
			key = matched_genes.find do |key|
				@ref_genes[key]["completeness"] == "complete"
			end
			key = matched_genes[0] if key.nil? || key.empty?

		elsif matched_genes.empty?
			key = @ref_genes.keys[0]
		end
		seq = Helper::Sequence.extract_translation(@ref_genes[key]["gene"])

		return key, seq
	end

	def write_refprot
		file = File.join(@file_basename, "#{@prot_basename}-refseq.fasta")
		File.open(file, 'w') do |f| 
			f.write(Helper::Sequence.str2fasta(@ref_key, @ref_seq))
		end
		return file
	end

	def get_refspecies
		@ref_genes[@ref_key]["species"]
	end
end