module PredictionsHelper
	Line_width = 80
	Space = "&nbsp;"

	# fill in alignment options from pairAlign into selection_tag
	def fill_selection_tags()
		arr_algo = [['Needleman-Wunsch', 'nw'], ['Gotoh', 'gotoh'], ['Smith-Waterman', 'sw'], ['Longest Common Subsequence', 'lcs'], ['DIALIGN2', 'dialign']]
		pre_select_algo = 'gotoh'
		arr_conf = [['End gap free', 'tttt'], ['Standard global alignment', 'ffff']]
		pre_select_conf = 'tttt'
		return arr_algo, pre_select_algo, arr_conf, pre_select_conf
	end

	# return pretty formatted ctg-usage
	# input:
		# predictd protein sequence,
		# CTG positions
		# [hash listing counts for each amino acid, total number of amino acids]
	# output:
		# array with each element containing one sequence part of Line_width, with siginifcant positions highlighted
		# array containing positions (as occuring in respecitive part of the sequence-array)
		# array containing amino acid statistics for each significant position

	def format_seq(prot, ctgpos, aa_stats)
		str, aa_pct, side_str = [], [], []

		# mark all significant positions by "!" to find them unambigously
		seq = String.new(prot)
		ctgpos.each {|pos| seq[pos] = "!"}
		ind = 0 # index of ctgpos, aa_stats

		seq.scan(/.{1,#{Line_width}}/).each do |part|
			# indices = part.enum_for(:scan,/!/).map { Regexp.last_match.begin(0) }
			# indices.each do |i|
			side = []
			part.count("!").times do |i|
				pos = ctgpos[ind]
				part.sub!("!", content_tag(:span, prot[ctgpos[ind]], :class => "highlight_aa", :title => "Position #{pos}"))
				aa_pct << format_aa_stats(aa_stats[ind][0])
				side << pos
				ind = ind + 1
			end
			str << part
			side_str << side.join(", ")
		end
		return str, side_str, aa_pct
	end

	# puts data into a html-table
	# input: 
		# array with data for column 1 (e.g. position)
		# array with data for column 2 (e.g. sequence data)
		# OPTIONAL: array with data for column 3 (e.g. additional data)
		# hash specifing classes for the table
	# output:
		# html-formatted tab
	def draw_table(col1, col2, *col3, table_spec)
		content_tag(:table, :class => table_spec[:table_class]) do
			content_tag(:tr) do
				content_tag(:th, table_spec[:th_left]) +
				if col3.any? then
					content_tag(:th, table_spec[:th_middle])
				end +
				content_tag(:th, table_spec[:th_right])
			end +
			content_tag(:colgroup) do
				content_tag(:col, "", :class => table_spec[:col_left]) +
				if col3.any? then
					content_tag(:col, "", :class => table_spec[:col_middle]) 
				end +
				content_tag(:col, "", :class => table_spec[:th_right])
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

	# return pretty formatted amino acid statistics
	# input: 
		# statistics for one ctg position
	# output: 
		# string, sorted by occurence of amino acids
		# and sort amino acids by percent
	def format_aa_stats(stats)
		str = []
		stats.sort_by {|k, v|v}.reverse.each do |aa, freq|
			if freq >= 0.05 then
				str << (aa + ": " + (freq*100).round.to_s + "%")
			end
		end
		return str.join(", ")
	end

	# return formatted ctg usage in reference data -statistics
	# input:
		# array of hashes with percentage of CTGs encoding Ser and Leu for each position
	# output:
		# array of strings
	def format_ctg_stats(data)
		data.map do |ele| 
			sum = (ele.values.sum * 100).round
			if sum > 0 then
				sum.to_s << "% total (" <<
				(ele[:ser]*100).round.to_s << "% encoding Ser, " <<
				(ele[:leu]*100).round.to_s << "% encoding Leu)"
			else
				sum.to_s << "%"
			end
		end
	end


end

	# def format_output_bkp(prot, ctgpos, aa_stats)
	# 	str, aa_pos, aa_pct, side_str = [], [], [], []

	# 	# mark all significant positions by "!" to find them unambigously
	# 	seq = String.new(prot)
	# 	ctgpos.each {|pos| seq[pos] = "!"}
	# 	ind = 0 # index of ctgpos, aa_stats

	# 	seq.scan(/.{1,#{Line_width}}/).each do |part|
	# 		# indices = part.enum_for(:scan,/!/).map { Regexp.last_match.begin(0) }
	# 		# indices.each do |i|
	# 		side = []
	# 		part.count("!").times do |i|
	# 			pos = ctgpos[ind]
	# 			part.sub!("!", "<span title=\"Position: #{pos}\" class=highlight_aa>#{prot[ctgpos[ind]]}</span>")
	# 			aa_pos << ("Position: " + pos.to_s)
	# 			aa_pct << format_aa_stats(aa_stats[ind])
	# 			side << pos
	# 			ind = ind + 1
	# 		end
	# 		str << part
	# 		side_str << side.join(", ")
	# 	end
	# 	return str.join("<br />"), aa_pos.join("<br />"), aa_pct.join("<br />"), side_str.join("<br />")
	# end
