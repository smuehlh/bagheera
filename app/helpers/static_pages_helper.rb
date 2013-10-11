module StaticPagesHelper

	require 'helper.rb'
	require 'proteinFamily.rb'

	def species_table
		coldata = refdata2species
		content_tag(:thead, escape: false) do
			content_tag(:tr) do 
				content_tag(:th, "Abbreviation") +
				content_tag(:th, "Species (teleomorph") +
				content_tag(:th, "Anamorph") +
				content_tag(:th, "Alternative names")
			end
		end +
		content_tag(:tbody) do
			coldata.keys.sort.each.collect do |abbr|
				data = coldata[abbr]
				content_tag(:tr) do
					col1_safe = abbr.html_safe? ? abbr : abbr.html_safe
					col2_safe = data[0]
					col3_safe = data[1]
					col4_safe = data[2]
					content_tag(:td, col1_safe) + 
					content_tag(:td, col2_safe) + 
					content_tag(:td, col3_safe) +
					content_tag(:td, col4_safe)
				end
			end.join.html_safe
		end

	end

	def refdata2species
		ref_data = Helper.load_ref_data

		spabbr2xxx = {}

		ref_data.each do |key, data|
			genes = data["genes"]

			genes.each do |abbr, genedat|
				sp_abbr = abbr.match(/
				([A-Z][_a-z]+) # match species abbreviation -- everything until the 2. uppercase char
				[A-Z][\w]+ # protein class and variant -- everything after 2. uppercase char 
				/x)[1]

				if ! spabbr2xxx.has_key?(sp_abbr) then
					begin
						anamorph = genedat["species_anamorph"] || ""
						altname = genedat["species_altname"] || ""
						name = genedat["species"]
					rescue
						anamorph = "?" if ! anamorph
						altname = "?" if ! altname
						name = "?" if ! name
					ensure
					spabbr2xxx[sp_abbr] = [name.html_safe, anamorph.html_safe, altname.html_safe]
					end
				end
			end
		end

		return spabbr2xxx

	end
end
