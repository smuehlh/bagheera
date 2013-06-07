module PredictionsHelper
	Line_width = 80
	Space = "&nbsp;"

	# fill in alignment options from pairAlign, dialign and Augustus (species) into selection_tag
	# @return [Array] Available algorithms (for alignment) methods in selection_tag format
	# @return [String] Preselected algorithm
	# @return [Array] Available algorithm options (for alignment) in selection_tag format
	# @return [String] Preselected configuration
	# @return [Array] Available species (for gene prediction) in selection_tag format
	# @return [String] Preselected species
	def fill_selection_tags()
		arr_algo = [['MAFFT (recommended)', 'mafft'],['Needleman-Wunsch', 'nw'], ['Gotoh', 'gotoh'], ['Smith-Waterman', 'sw'], ['Longest Common Subsequence', 'lcs']]
		pre_select_algo = 'mafft'
		arr_conf = [['End gap free', 'tttt'], ['Standard global alignment', 'ffff']]
		pre_select_conf = 'tttt'
		arr_species = [
			['Candida albicans', 'candida_albicans'],
			['Candida guilliermondii', 'candida_guilliermondii'],
			['Candida tropicalis', 'candida_tropicalis'],
			['Debaryomyces hansenii', 'debaryomyces_hansenii'],
			['Eremothecium gossypii', 'eremothecium_gossypii'],
			['Kluyveromyces lactis', 'kluyveromyces_lactis'],
			['Lodderomyces elongisporus', 'lodderomyces_elongisporus'],
			['Pichia stipitis', 'pichia_stipitis'],
			['Saccharomyces cerevisiae S288C', 'saccharomyces_cerevisiae_S288C'],
			['Saccharomyces cerevisiae RM11-1a', 'saccharomyces_cerevisiae_rm11-1a_1'],
			['Yarrowia lipolytica', 'yarrowia_lipolytica']]

		# arr_species = [
		# 	['Homo sapiens', 'human'], 
		# 	['Drosophila melanogaster', 'fly'], 
		# 	['Arabidopsis thaliana', 'arabidopsis'],
		# 	['Brugia malayi', 'brugia'],
		# 	['Aedes aegypti', 'aedes'],
		# 	['Tribolium castaneum', 'tribolium'],
		# 	['Schistosoma mansoni', 'schistosoma'],
		# 	['Tetrahymena thermophila', 'tetrahymena'],
		# 	['Galdieria sulphuraria', 'galdieria'],
		# 	['Zea mays', 'maize'],
		# 	['Toxoplasma gondii', 'toxoplasma'],
		# 	['Caenorhabditis elegans', 'caenorhabditis'],
		# 	['Aspergillus fumigatus', 'aspergillus_fumigatus'],
		# 	['Aspergillus nidulans', 'aspergillus_nidulans'],
		# 	['Aspergillus oryzae', 'aspergillus_oryzae'],
		# 	['Aspergillus terreus', 'aspergillus_terreus'],
		# 	['Botrytis cinerea', 'botrytis_cinerea'],
		# 	['Candida albicans', 'candida_albicans'],
		# 	['Candida guilliermondii', 'candida_guilliermondii'],
		# 	['Candida tropicalis', 'candida_tropicalis'],
		# 	['Chaetomium globosum', 'chaetomium_globosum'],
		# 	['Coccidioides immitis', 'coccidioides_immitis'],
		# 	['Coprinus cinereus', 'coprinus_cinereus'],
		# 	['Cryptococcus neoformans gattii', 'cryptococcus_neoformans_gattii'],
		# 	['Cryptococcus neoformans neoformans B', 'cryptococcus_neoformans_neoformans_B'],
		# 	['Cryptococcus neoformans neoformans JEC21', 'cryptococcus_neoformans_neoformans_JEC21'],
		# 	['Debaryomyces hansenii', 'debaryomyces_hansenii'],
		# 	['Encephalitozoon cuniculi', 'encephalitozoon_cuniculi_GB'],
		# 	['Eremothecium gossypii', 'eremothecium_gossypii'],
		# 	['Fusarium graminearum', 'fusarium_graminearum'],
		# 	['Histoplasma capsulatum', 'histoplasma_capsulatum'],
		# 	['Kluyveromyces lactis', 'kluyveromyces_lactis'],
		# 	['Laccaria bicolor', 'laccaria_bicolor'],
		# 	['Petromyzon marinus', 'lamprey'],
		# 	['Leishmania tarentolae', 'leishmania_tarentolae'],
		# 	['Lodderomyces elongisporus', 'lodderomyces_elongisporus'],
		# 	['Magnaporthe grisea', 'magnaporthe_grisea'],
		# 	['Neurospora crassa', 'neurospora_crassa'],
		# 	['Phanerochaete chrysosporium', 'phanerochaete_chrysosporium'],
		# 	['Pichia stipitis', 'pichia_stipitis'],
		# 	['Rhizopus oryzae', 'rhizopus_oryzae'],
		# 	['Saccharomyces cerevisiae S288C', 'saccharomyces_cerevisiae_S288C'],
		# 	['Saccharomyces cerevisiae RM11-1a', 'saccharomyces_cerevisiae_rm11-1a_1'],
		# 	['Schizosaccharomyces pombe', 'schizosaccharomyces_pombe'],
		# 	['Trichinella spiralis', 'trichinella'],
		# 	['Ustilago maydis', 'ustilago_maydis'],
		# 	['Yarrowia lipolytica', 'yarrowia_lipolytica'],
		# 	['Nasonia vitripennis', 'nasonia'],
		# 	['Solanum lycopersicum', 'tomato'],
		# 	['Chlamydomonas reinhardtii', 'chlamydomonas'],
		# 	['Amphimedon queenslandica', 'amphimedon'],
		# 	['Pneumocystis jirovecii', 'pneumocystis']]
		pre_select_species = 'candida_albicans'
		return arr_algo, pre_select_algo, arr_conf, pre_select_conf, arr_species, pre_select_species
	end

	# Format data for view
	# @param prot [String] predicted protein sequence
	# @param ctgpos [Array] Found CTG positions in the predicted sequence
	# @param chem_props [Hash] Hash of hashes; amino acids stats (from reference data) for each significant CTG position
	# @return [Array] Formatted protein sequence (including tags to highlight significant positions); 80 chars per element
	# @return [Array] Formatted CTG positions; elements correspond to the formatted seq
	# @return [Array] Formatted amino acid statistics for each significant position
	# @return [Array] Positions in other return array of non-significant positions -> display those in gray!
	def format_seq(prot, ctgpos, chem_props) 
		str, aa_pct, side_str, not_sign = [], [], [], []
		prot.gsub!("!", "") # it should really not contain any "!", but this may happen with dialign

		# mark all significant positions by "!" to find them unambigously
		seq = String.new(prot)
		ctgpos.each {|pos| seq[pos] = "!"}
		ind = 0 # index of ctgpos

		seq.scan(/.{1,#{Line_width}}/).each do |part|
			# indices = part.enum_for(:scan,/!/).map { Regexp.last_match.begin(0) }
			# indices.each do |i|
			side = []
			part.count("!").times do |i|
				pos = ctgpos[ind]
				part.sub!("!", content_tag(:span, prot[ctgpos[ind]], :class => "highlight_aa", :title => "Position " + ruby2human_counting(pos).to_s ))
				if chem_props.has_key?(pos) then
					aa_pct << format_aa_stats(chem_props[pos][:aa_comp])
					not_sign << ind if ! chem_props[pos][:is_significant]
				end

				side << ruby2human_counting(pos)
				ind = ind + 1
			end
			str << part
			side_str << side.join(", ")
		end
		return str, side_str, aa_pct, not_sign
	end

	# Puts data into a html-table
	# @param col1 [Array] Data for column 1 (e.g. position)
	# @param col2 [Array] Data for column 2 (e.g. sequence data)
	# @param col3 [Array] Data for column 3; optional
	# @param table_spec [Hash] Specifies table properties like column headings
	# @param has_table_tag [Boolean] specifies if table has surrounding <table>-Tags; Default: true
	# @return [String] html-formatted table
	def draw_table(col1, col2, *col3, table_spec)

		content_tag(:table, :class => table_spec[:table_class]) do

			content_tag(:tr) do
				content_tag(:th, table_spec[:th_left]) +
				if col3.any? then
					content_tag(:th, table_spec[:th_middle])
				end +
				content_tag(:th, table_spec[:th_right])
			end +
			if  ! table_spec[:col_left].empty? || 
				! table_spec[:col_right].empty? || 
				(col3.any? && ! table_spec[:col_middle].empty?) 
				then
				content_tag(:colgroup) do
					content_tag(:col, "", :class => table_spec[:col_left]) +
					if col3.any? then
						content_tag(:col, "", :class => table_spec[:col_middle]) 
					end +
					content_tag(:col, "", :class => table_spec[:col_right])
				end
			end +
			col1.each_with_index.collect do |data, ind|
				content_tag(:tr) do
					col1_safe = data.html_safe? ? data : data.html_safe
					col2_safe = col2[ind].html_safe? ? col2[ind] : col2[ind].html_safe
					content_tag(:td, col1_safe) + 
					content_tag(:td, col2_safe) +
					if col3.any? then
						col3.flatten!
						col3_safe = col3[ind].html_safe? ? col3[ind] : col3[ind].html_safe
						content_tag(:td, col3_safe)
					end
				end
			end.join.html_safe

		end
	end

	# no table-tag is printed, useful to concatenate tables
	def draw_table2(col1, col2, *col3, table_spec)
		content_tag(:thead) do
			content_tag(:tr) do
				content_tag(:th, table_spec[:th_left]) +
				if col3.any? then
					content_tag(:th, table_spec[:th_middle])
				end +
				content_tag(:th, table_spec[:th_right])
			end
		end +

		content_tag(:tbody) do
			content_tag(:i, "") + # need this tag for the following if/else!
			if  ! table_spec[:col_left].empty? || 
				! table_spec[:col_right].empty? || 
				(col3.any? && ! table_spec[:col_middle].empty?) 
				then
				content_tag(:colgroup) do
					content_tag(:col, "", :class => table_spec[:col_left]) +
					if col3.any? then
						content_tag(:col, "", :class => table_spec[:col_middle]) 
					end +
					content_tag(:col, "", :class => table_spec[:col_right])
				end
			end +

			col1.each_with_index.collect do |data, ind|
				# set class for table row -> display not significant data in grey
				if table_spec.has_key?(:display_grey) then
					if table_spec[:display_grey].include?(ind) then
						klass = "grey"
						change_col = "this.style.color='red';"
					else
						klass = ""
					end
				else
					klass = ""
				end
				content_tag(:tr, :class => klass) do
					col1_safe = data.html_safe? ? data : data.html_safe
					col2_safe = col2[ind].html_safe? ? col2[ind] : col2[ind].html_safe
					content_tag(:td, col1_safe) + 
					content_tag(:td, col2_safe) +
					if col3.any? then
						col3.flatten!
						col3_safe = col3[ind].html_safe? ? col3[ind] : col3[ind].html_safe
						content_tag(:td, col3_safe)
					end
				end
			end.join.html_safe
		end
	end

	# formatting amino acid statistics
	# @param stats [Array] Amino acid statistics for one CTG-pos
	# @param is_reject_low [Boolean] show only stats > 5% or all ; not mandatory; default = true
	# @return [String] Format: "L: 100%", sorted by percent; only > 5%
	def format_aa_stats(stats, is_reject_low=true)
		res = []
		has_rejected = false 
		stats.sort_by {|k, v|v}.reverse.each do |aa, freq|
			if freq >= 0.05 && is_reject_low then
				res << (aa + ": " + (freq*100).round.to_s + "%")
			elsif freq < 0.05 && is_reject_low 
				has_rejected = true
			end
		end
		res << "Others: < 5%" if has_rejected
		# check if any freq is < 0.05 than add "Others: , 5%"
		return res.join(", ")
	end

	# format statistics over CTG usage in reference data
	# @param [Array] Hash; CTG usage in reference data in percentages
	# @return [Array] Formatted statistics
	def format_ctg_stats(data)
		res = []
		data.each do |_, val|
			str = "S: " << 
				(val[:ctg_usage]["S"]*100).round.to_s << 
				"%, L: " << 
				(val[:ctg_usage]["L"]*100).round.to_s << "%"
			res << str
		end
		return res
	end

	# converts ruby counting to human counting: add 1 to each value
	# can handle both Arrays and Fixnums
	# @param lnum [Array or Fixnum] Ruby counted numbers
	# @return [Array or Fixnum] Human counted numbers
	def ruby2human_counting(num)
		if num.kind_of?(Array) then
			numplus = num.collect do |n|
				next if n.blank? || n.nil?
				n += 1
			end
			return numplus
		elsif num.kind_of?(Fixnum) 
			return num + 1
		end
	end


	def suggested_transl(ref_chem, ref_ctg)
		res = []

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
		strange = (pos_ser & pos_leu).size # set intersect

		# only add to results if values are not zero
		if alt_usage > 0 then
			res << pluralize(alt_usage, ' CTG position', ' CTG positions') + " suggest alternative codon usage."
		end
		if std_usage > 0 then
			res << pluralize(std_usage, ' CTG position', ' CTG positions') + " suggest standard codon usage."
		end
		if strange > 0 then
			res << pluralize(strange, ' CTG position', ' CTG positions') +
				" leads to contrary results." 
				# " suggest contrary codon usage based on distribution of amino acids and CTG usage in reference data."
		end

		return res.join("</br>").html_safe
	end

end