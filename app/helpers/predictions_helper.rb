module PredictionsHelper

	# return pretty formatted ctg-usage
	# input:
		# predictd protein sequence,
		# CTG positions
		# hash listing counts for each amino acid
	# output:
		# str [div at top]: sequence, Line_width chars per line, significant positions marked, mouse over
		# aa_table [div bottom]
			# aa_pos [div bottom-left]
			# aa_pct [div bottom-right]

	def format_output(prot, ctgpos, aa_stats)
		str, aa_pos, aa_pct = [], [], []

		# mark all significant positions by "!" to find them unambigously
		seq = String.new(prot)
		ctgpos.each {|pos| seq[pos] = "!"}
		ind = 0 # index of ctgpos, aa_stats

		seq.scan(/.{1,#{Line_width}}/).each do |part|
			# indices = part.enum_for(:scan,/!/).map { Regexp.last_match.begin(0) }
			# indices.each do |i|
			part.count("!").times do |i|
				pos = ctgpos[ind]
				part.sub!("!", "<span title=\"Position: #{pos}\" class=highlight_aa>#{prot[ctgpos[ind]]}</span>")
				aa_pos << ("Position: " + pos.to_s)
				aa_pct << format_aa_stats(aa_stats[ind])
				ind = ind + 1
			end
			str << part
		end
		return str.join("<br />"), aa_pos.join("<br />"), aa_pct.join("<br />")
	end

	# return pretty formatted amino acid statistics
	# input: 
		# statistics for one ctg position
	# output: 
		# string, convert count to percent
		# and sort amino acids by percent
	def format_aa_stats(stats)
		str = []
		num_aas = stats.inject(0) {|sum, (_,val)| sum + val}
		stats.sort_by {|k, v|v}.reverse.each do |aa, count|
			pct = (100 * count/num_aas.to_f).round
			if pct >= 5 then
				str << (aa + ": " + pct.to_s + "%")
			end
		end
		return str.join(", ")
	end

	# return pretty formatted ctg-usage
	# input: 
		# predicted protein sequence, 
		# positions of CTGs, 
		# percent of polar amino acids in reference alignment, 
		# percent of hydrophob amino acids in reference alignment
	# output:
		# str [left div]: sequence, 80 chars per line, discriminative positions marked
		# pct_str [right div]: 
			# 1 discrim. pos: percent polar | percent hydro
			# >=1 discrim. pos: marker for line of optional_line
		# optional_line [div underneath]: line \n pct_polar \n pct_hydro
	Line_width = 80
	Space = "&nbsp;"

	def format_output_old(prot, ctgpos, pct_polar, pct_hydro)
		data = []
		data << prot.split("")
		str, pct_str = "", ""

		# replace all discrim. positions by marker "!" and 
		ctgpos.each {|pos| data[0][pos] = "!"}
		ind = 0 # index of ctgpos, pct_polar, pct_hydrophob array

		(data[0].join("").scan /.{1,#{Line_width}}/).each do |part|
			str += part
			indices = part.enum_for(:scan,/!/).map { Regexp.last_match.begin(0) }
			if indices.length == 0 then
				# ensure that pctages are written next to appropriate line
				pct_str += "<br />"
				str += "<br />"
			# elsif indices.length == 1 	
			# 	pct_str += "<span title=\"polar [%]\">" + (pct_polar[ind] * 100).round.to_s +
			# 		"</span> | <span title=\"hydrophobic [%]\">" +
			# 		(pct_hydro[ind] * 100).round.to_s + "</span><br />"
			# 	# replace marker by actual amino acid
			# 	str.gsub!("!", "<span class=highlight_aa>#{prot[ctgpos[ind]]}</span>")
			# 	ind += 1	
			else
				indices.each do |i|
					pct_str += "<span title=\"polar [%]\">" + (pct_polar[ind] * 100).round.to_s +
						"</span> | <span title=\"hydrophobic [%]\">" +
						(pct_hydro[ind] * 100).round.to_s + "</span><br />"
					# replace marker by actual amino acid
					str.gsub!("!", "<span class=highlight_aa>#{prot[ctgpos[ind]]}</span>")
					ind += 1	
					str += "<br />"		
				end
			end
		end
		return str, pct_str
	end
end
