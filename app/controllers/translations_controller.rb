class TranslationsController < ApplicationController

	require 'helper.rb'

	Codon_table = {
		"TTT" => "F", "TCT" => "S", "TAT" => "Y", "TGT" => "C", 
		"TTC" => "F", "TCC" => "S", "TAC" => "Y", "TGC" => "C",
		"TTA" => "L", "TCA" => "S", "TAA" => "*", "TGA" => "*",
		"TTG" => "L", "TCG" => "S", "TAG" => "*", "TGG" => "W",
		"CTT" => "L", "CCT" => "P", "CAT" => "H", "CGT" => "R",
		"CTC" => "L", "CCC" => "P", "CAC" => "H", "CGC" => "R",
		"CTA" => "L", "CCA" => "P", "CAA" => "Q", "CGA" => "R",
		"CTG" => "L", "CCG" => "P", "CAG" => "Q", "CGG" => "R",
		"ATT" => "I", "ACT" => "T", "AAT" => "N", "AGT" => "S",
		"ATC" => "I", "ACC" => "T", "AAC" => "N", "AGC" => "S",
		"ATA" => "I", "ACA" => "T", "AAA" => "K", "AGA" => "R", 
		"ATG" => "M", "ACG" => "T", "AAG" => "K", "AGG" => "R", 
		"GTT" => "V", "GCT" => "A", "GAT" => "D", "GGT" => "G",
		"GTC" => "V", "GCC" => "A", "GAC" => "D", "GGC" => "G",
		"GTA" => "V", "GCA" => "A", "GAA" => "E", "GGA" => "G",
		"GTG" => "V", "GCG" => "A", "GAG" => "E", "GGG" => "G" }
	Unknown_codon = "X"

	def start
		# renders the start-html page
	end

	# translates mRNA sequence into protein sequence
	# last nucleotides not belonging to any codon (e.g. if length mRNA % 3 != 0) are not translated
	# input: params[:mrna_seq] [String] The mRNA sequence to be translated, will be transformed to uppercase and DNA seq
	# input: params[:codonusage] [String] "L" or "S" depending on how CTG codon should be translated
	# accessible params in view: @protein_seq [String] Translated protein sequence
	# accessible params in view: @ctg_pos [Array] CTG positions in @protein_seq
	# accessible params in view: @error [String] Error method, or empty if no error occured
	# accessible params in view: @is_from_transl_mrna [False] Indicates that results should be rendered as trans_mrna results

	def transl_mrna
		@is_from_transl_mrna = true
		@protein_seq = ""
		@ctg_pos = []
		@error = ""
		
		if params[:mrna_seq].empty? then 
			@error = "No sequence provided."
		end

		mrna = Helper::Sequence.validate_dna_seq(params[:mrna_seq])

		transl_ctg_as = params[:codonusage]
		translation_length = mrna.length
		translation_length -= translation_length % 3
		
		# translate each codon into the correct amino acid, and keep track of CTG codons
		@protein_seq = mrna[0...translation_length].gsub(/.{3}/).with_index do |codon, ind| 
			aa = Codon_table[codon] 
			aa ||= Unknown_codon # set aa to unknown codon if it has no value
			if codon == "CTG" then 
				# collect CTG position
				@ctg_pos.push(ind)
				if transl_ctg_as == "S" then 
					# correct translation if it should be "S" instead of default "L"
					aa = "S"
				end
			end
			aa # needed in order to collect it
		end

		# remove trailing stop codon (its unusual to report it...)
		@protein_seq = @protein_seq.gsub(/\*$/, '')

	rescue NoMethodError, TypeError, NameError, RuntimeError, Errno::ENOENT
		@error = "Cannot translate mRNA into protein."

	ensure
		render :show_protein_seq, formats: [:js]

	end

	# translates protein sequence into mRNA and checks if CTG codons were correctly translated in protein sequence
	# input: params[:protein_seq] [String] protein sequence to check
	# input: params[:species] [String] species abbreviation of the species to use for gene prediction
	# accessible params in view: @protein_seq [String] (Correctly) Translated protein sequence
	# accessible params in view: @ctg_pos [Array] CTG positions in @protein_seq
	# accessible params in view: @ctg_pos_wrongly_transl [Array] CTG positions in @protein_seq wrongly transl. in input seq
	# accessible params in view: @error [String] Error method, or empty if no error occured
	# accessible params in view: @is_from_transl_mrna [False] Indicates that results should be rendered as protein-transl results
	def transl_protein

		require 'webScipio.rb'

		@is_from_transl_mrna = false
		@protein_seq = "" # correct translated
		@ctg_pos = []
		@ctg_pos_wrongly_transl = [] # positions which were wrongly translated
		@error = ""

		if params[:protein_seq].empty? then 
			@error = "No sequence provided."
		end
		@protein_seq = Helper::Sequence.validate_seq(params[:protein_seq])
		fasta_formatted_protein_seq = Helper::Sequence.str2fasta(">uploaded_seq", @protein_seq)
		species_name = params[:species] || ""
		is_use_alternative_codon_usage = view_context.is_species_with_alternative_codon_usage(species_name)

		if params[:scipio_relaxed] then 
			# don't use default if "use relaxed" option is set
			is_use_scipio_default_params = false
		else
			# otherwise use default (standard)
			is_use_scipio_default_params = true
		end

		# call webscipio
		yaml, @error = WebScipio.call_webscipio( species_name, fasta_formatted_protein_seq, is_use_scipio_default_params )

		if @error.blank? then 

			cdna = Helper::Sequence.extract_cdna(yaml)
			cdna.scan(/.{1,3}/).each_with_index do |codon, ind|
				if codon == "CTG" then 
					@ctg_pos.push(ind)

					transl_input = @protein_seq[ind] # translation in input sequence
					if is_use_alternative_codon_usage then 
						# translation according to codon usage of species
						transl_correct = "S"
					else
						transl_correct = "L"
					end

					if transl_input != transl_correct then 
						@ctg_pos_wrongly_transl.push(ind)
						@protein_seq[ind] = transl_correct
					end					

				end
			end

		end

	rescue NoMethodError, TypeError, NameError, RuntimeError, Errno::ENOENT
		@error = "An error occured. Cannot check CTG translation in protein."

	ensure
		render :show_protein_seq, formats: [:js]

	end

	# Creates a temporary file containing the protein sequence for download
	def download
		file_basename = File.dirname(session[:file][:path]) 
		file = File.join(file_basename, Helper.get_tmp_file)
		fasta = Helper::Sequence.str2fasta(">translated protein", params[:seq])
		File.open(file, 'w') {|f| f.write(fasta)}

		send_file file, :x_sendfile=>true

	rescue RuntimeError, Errno::ENOENT => exc
		# sufficient to initialize @error for template :show_tree, because if this is not empty, only error message will be shown
		@error = "Cannot prepare file for download."
		render :show_protein_seq, formats: [:js]
	end

end
