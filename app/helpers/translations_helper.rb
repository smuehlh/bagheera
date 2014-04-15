module TranslationsHelper
	attr_accessor :all_species_by_codon_usage

	@@all_species_by_codon_usage = {}

	require 'helper.rb'

	def fill_codon_usage_selection_tag
		options = [
			["Standard", "L"], 
			["Alternative yeast nuclear", "S"]
		]
		pre_selected = "L"
		return options, pre_selected
	end

	# provides species grouped by their translation scheme
	# @returns Hash [Species Name]
	# @returns String [Species name] - used as preselected value
	def fill_and_group_species_selection_tag()

		grouped_options = { "Standard codon usage" => [], "Alternative codon usage" => [] }
		pre_selected = "Saccharomyces cerevisiae S288c" # assuming that Saccharomyces is always under reference species
		seen_species = {} # hash containing all species already seen

		ref_data = Helper.load_ref_data

		ref_data.each do |key, data|
			genes = data["genes"]

			genes.each do |abbr, genedat|
				sp_abbr = abbr.match(/
				([A-Z][_a-z]+) # match species abbreviation -- everything until the 2. uppercase char
				[A-Z][\w]+ # protein class and variant -- everything after 2. uppercase char 
				/x)[1]

				if ! seen_species.has_key?(sp_abbr) then 
					seen_species[sp_abbr] = "" # add key to mark this species abbreviation as seen 
					begin 
						name = genedat["species"]
						codonusage = genedat["codonusage"].to_s || ""
					rescue
						next
					end
	
					# codonusage will be NCBI tag for standard (1) or alternative (12) usage
					if codonusage == "1" then 
						grouped_options["Standard codon usage"].push( [name] )
					elsif codonusage == "12" then 
						grouped_options["Alternative codon usage"].push( [name] )
					else
						if ! grouped_options.has_key?("Unknown codon usage") then 
							grouped_options["Unknown codon usage"] = []
						end
						grouped_options["Unknown codon usage"].push( [name] ) 	
					end
				end
			end
		end
		@@all_species_by_codon_usage = grouped_options

		return grouped_options, pre_selected
	end

	# Format protein sequence (with highlighted CTG positions) as table
	# @param protein_seq [String] protein sequence
	# @param ctg_pos [Array] Found CTG positions in the sequence
	# @return [Array] Formatted protein sequence (including tags to highlight significant positions); 80 chars per element
	# @return [Array] Formatted CTG positions; elements correspond to the formatted seq
	def pretty_format_protein_seq(protein_seq, ctg_pos)

		# call method from PredictionsHelper
		seq_string, ctg_string, dummy = format_seq(protein_seq, ctg_pos, {})

        return seq_string, ctg_string
	end

	# Check for a given species if it uses alternative codon usage
	# @param species_name [String] species name as submitted by select tag
	# @return [Boolean] True if species has alternative codon usage, false otherwise
	def is_species_with_alternative_codon_usage(species_name)
		@@all_species_by_codon_usage["Alternative codon usage"].flatten.include?(species_name)
	end

end
