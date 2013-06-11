# this class is for evaluation purposes
	# call me like this: 
	# rails c 
	# > Evaluation.eval_ref_data
class Evaluation < PredictionsController
	extend PredictionsHelper

	# evaluate reference data
	# get statistics about ALL CTG positions in ref data - takes some time!
	# for all proteins
		# for all genes
			# - CTG positions
			# - amino acid usage of other genes at this position
			# - CTGs in other genes at this position

	def Evaluation.eval_ref_data
 		
		fh = File.new("/tmp/cug/statistics_about_ref_data.txt", "w")
		ref_data, fatal_error = PredictionsController.new.load_ref_data

		if fatal_error.empty? then
			ref_data.keys.sort_by { |key| key.to_s.naturalized }.each do |prot|
			puts prot	
			# statistics about all CTG-Positions in reference data
			results = {ref_seq_num: "", ctg_pos: [], ref_chem: {}, ref_ctg: {}}
			pos2key = {} # map from ctg position to a sequece containing this ctg

			# get all CTG positions in all reference sequences
			ref_data[prot]["genes"].keys.each do |k|
				dna_seq = PredictionsController.new.extract_gene_seq(ref_data[prot]["genes"][k]["gene"])
				codons = dna_seq.scan(/.{1,3}/)
				ctg_pos = codons.each_with_index.map{ |ele, ind| (ele == "CTG") ? ind : nil }.compact

				# store ctg position
				results[:ctg_pos] |= ctg_pos # set union

				# store sequence key if it is a "new" ctg position
				ctg_pos.each {|pos| pos2key[pos] = k if ! pos2key.has_key?(pos)}

			end

			# extract ref_alignment 
			begin
			ref_alignment = Hash[*ref_data[prot]["alignment"].split("\n")]
			ref_alignment.keys.each {|k| ref_alignment[ k.sub(">", "") ] = ref_alignment.delete(k)}
			rescue => e
				puts "---"
				puts "ERROR (#{prot}):"
				puts e
				puts "---"
				next
			end

			results[:ref_seq_num] = ref_alignment.keys.size

			# get list of all amino acids and codons at all CTG -Positions
			ref_alignment_cols = []
			ref_codons = []

			results[:ctg_pos].each do |pos|
				# all amino acids and all codons for a given CTG position
				col = []
				codons = []
				
 				# get reference key of a sequence containing this CTG position
 				ref_key = pos2key[pos]
				ref_cymopos = PredictionsController.new.sequence_pos2alignment_pos(pos, ref_alignment[ref_key])

				ref_alignment.each do |cymo_header, cymo_prot|

					# collect cymo column at CTG position
					col << cymo_prot[ref_cymopos]

					# collect codons at CTG position
					if cymo_prot[ref_cymopos] == "-" then
						# add nil as placeholder if its a gap
						codons << nil
					else
						# actually collect codon
						spos = PredictionsController.new.alignment_pos2sequence_pos(ref_cymopos, cymo_prot) # position in sequece without gaps
						codons << PredictionsController.new.get_codon(ref_data[prot]["genes"][cymo_header]["gene"], spos)
					end
				end

				# store results
				ref_alignment_cols << col
				ref_codons << codons.compact # if only nil- placeholders, it will be an empty array
			end

			# amino acid usage 
			ref_alignment_cols.each_with_index do |col, ind|
				# check for each column (= each CTG position) the amino acid distribution in reference data
				aa_freq, aa_num = PredictionsController.new.word_frequency(col)
				is_significant, prob_transl = PredictionsController.new.predict_translation(aa_freq)

				# possibly, here the mapping on the CTG position is wrong (results[:ctg_pos][ind])
				results[:ref_chem][results[:ctg_pos][ind]] = {transl: prob_transl, aa_comp: aa_freq, aa_num: aa_num}

			end

			# CTG usage
			results[:ref_ctg] = PredictionsController.new.ref_ctg_usage(ref_codons, ref_alignment_cols, results[:ctg_pos])

			# save statistics to file
			fh.puts prot
			fh.puts "#{results[:ref_seq_num]} sequences."
			fh.puts "#{results[:ctg_pos].size} CTG positions."
			fh.puts ""
			fh.puts "CTG position\tDistribution of amino acids\tNumber of amino acids\tCTG usage\tNumber of CTG codons"
			results[:ctg_pos].sort.each do |pos|
				# check if for all ctg positsions info about aa distribution and ctg usage exists
				if results[:ref_chem].has_key?(pos) then
					# call method from PredictionsHelper (formerly with: view_context.format_aa_stats, without extend module)
					aa_stat = format_aa_stats(results[:ref_chem][pos][:aa_comp])
					n_aas = results[:ref_chem][pos][:aa_num]
				else 
					aa_stat = "Not discriminative."
					n_aas = "-"
				end
				if results[:ref_ctg].has_key?(pos) then
					ctg_stat = format_aa_stats(results[:ref_ctg][pos][:ctg_usage])
					n_ctgs = results[:ref_ctg][pos][:ctg_num]
				else
					ctg_stat = "Not discriminative."
					n_ctgs = "-"
				end

				# tab separated list of all statistics for this position
				# pos + 1 to convert between ruby and human counting! 
				fh.puts "#{ruby2human_counting(pos)}\t#{aa_stat}\t#{n_aas}\t#{ctg_stat}\t#{n_ctgs}"
			end
			fh.puts ""


			end # ref_data.each
		end # if fatal_error.empty?
		fh.close
	end

	# parsing file containing all statistics into single csv -files for excel
	def Evaluation.parse_all_stats 
		file_in = "/tmp/cug/statistics_about_ref_data.txt"

		# out-file 1: Verteilung CTG Positionen pro Protein
		out_ctg_prot = "/fab8/smuehlh/data/cugusage/ctg_prot.csv"
		fh_ctg_prot = File.new("/fab8/smuehlh/data/cugusage/ctg_prot.csv", "w")
		fh_ctg_prot.puts "Protein,# CTGs"

		# # out-file 2: Verteilung CTG Positionen pro Organismus
		# out_ctg_org = "/fab8/smuehlh/data/cugusage/ctg_org.csv"
		# fh_ctg_org = File.new("/fab8/smuehlh/data/cugusage/ctg_org.csv", "w")
		# fh_ctg_org.puts "Organismus,# CTGs"

		prot = ""
		IO.foreach(file_in) do |line|
			debugger
			if ! line.include?(".") && ! line.include?("\t") then
				# its an protein
				prot = line.chomp
			elsif line.include?("CTG positions.") 
				# its the number of CTG positions
				n_ctg = line.match(/\d+/)[0]
				fh_ctg_prot.put
				s prot + "," + n_ctg
			end
		end

		fh_ctg_prot.close
		# fh_ctg_org.close

	end
end
