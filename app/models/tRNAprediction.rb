class TRNAprediction

	attr_reader :err_msg

	def initialize(path_to_file, path_to_ref_data)
		@file = path_to_file
		@base_path = File.dirname(path_to_file)
		@ref_data_path = path_to_ref_data
		@anticodon = "CAG"

		@tRNA_seq = ""
		@tRNA_sec_struct = ""
		@anticodon_pos = ""
		@score = 0 # score from tRNAscan
		@blast_results = {}
		@f_blast_alignment = ""
		@err_msg = ""
	end

	def predict

		@err_msg = catch(:problem) {

			# 1) scan in genome for CAG tRNA
			run_trna_scan

			# 2) extract tRNA seq
			@tRNA_seq, @tRNA_sec_struct, @anticodon_pos, @score = parse_trna_scan

			# 3) blast against CAG tRNAS with known identity
			blast

			# 4) parse blast hits: how many leu and ser trnas?
			@blast_results = parse_and_analyse_blasthits
			@f_blast_alignment = parse_blast_alignment

			# 5) align found tRNA with reference
			align

			# repeat what should be catched by err_msg
			""
		}

	end

	def run_trna_scan
		is_success = ProgCall.trnascan(@file)
		Helper.worked_or_throw_error(is_success, "tRNA scan failed.")
	end

	def run_blast(db)
		@f_blasthits = File.join(@base_path, "trna.blast")
		is_success = ProgCall.blast( @f_blasthits, @f_seq, {:genome_db => db, :program => "blastn", :output_format => "4"} )
		Helper.worked_or_throw_error(is_success, "BLAST failed.")
	end

	def run_mafft(f_in_pred, f_in_ref)
		is_success = ProgCall.mafft(@f_align, f_in_pred, f_in_ref, {:input => "nuc"})
		Helper.worked_or_throw_error(is_success, "Aligning sequences (MAFFT) failed.")
	end

	def blast
		# 1) setup blast-db
		db_name = File.join(@base_path, "trna_db")
		ProgCall.create_blast_db(@ref_data_path, db_name)

		# 2) blast against blast-db
		run_blast(db_name)

	end

	def align
		# output file
		@f_align = File.join(@base_path, "trna-align-aligned.fasta")
		# input file: @f_seq
		# input file (ref data)
		f_ref = @ref_data_path.gsub(".fasta", "_aligned.fasta")
		run_mafft(@f_seq, f_ref)
	end

	def parse_trna_scan
		ext_trna_scan = "sta.ss"

		f_out = File.join( @base_path,  File.basename(@file, ".*") + ext_trna_scan)
		is_success = Helper.does_file_exist(f_out)
		Helper.worked_or_throw_error(is_success, "tRNA scan failed.")

		# output = File.read(f_out)
		is_collect_seq = false
		is_collect_struct = false
		seq = ""
		struct = ""
		anticodon_pos = ""
		score = ""
		intron_start = nil
		intron_stop = nil


# gi|92090993|gb|CM000311.1|.trna6 (58090-58009)  Length: 82 bp
# Type: Leu       Anticodon: CAG at 34-36 (58057-58055)   Score: 53.85
# Possible intron: 38-63 (1101298-1101323)
#          *    |    *    |    *    |    *    |    *    |    *    |    *    |    *    | 
# Seq: GATACGATGGCCGAGTGGTtAAGGCGAAGGATGCAGGTTCCTTTGGGCATTGCCCGCGCAGGTTCGAACCCTGCTCGTGTCG
# Str: >>>>>>>..>>>..........<<<.>>>>>.......<<<<<.>>>>...<<<<..>>>>>.......<<<<<<<<<<<<.

		IO.foreach(f_out) do |line|
			if line.match(/Anticodon: #{@anticodon}/) then 
				is_collect_seq = true
				is_collect_struct = true
				anticodon_pos = line.match(/at ([\d-]+) /)[1]
				score = line.match(/Score: ([.\d]+)/)[1]
			end
			if is_collect_seq && line.match(/Possible intron/) then 
				matches = line.match(/(\d+)-(\d+)\s/)
				intron_start = matches[1].to_i - 1 # convert human to ruby counting
				intron_stop = matches[2].to_i -1 # convert human to ruby counting
			end
			if is_collect_seq && line.match(/^Seq:/) then 
				seq = line.match(/Seq: (.*)/)[1]
				seq = seq.upcase
				# remove intron
				if intron_start && intron_stop then 
					seq.slice!(intron_start..intron_stop)
				end
				is_collect_seq = false
			end
			if is_collect_struct && line.match(/^Str:/) then 
				struct = line.match(/Str: (.*)/)[1]
				# remove intron
				if intron_start && intron_stop then 
					struct.slice!(intron_start..intron_stop)
				end
				is_collect_struct = false
			end
			if ! seq.blank? && ! struct.blank? then 
				# important: stop search
				break
			end
		end

		if seq == "" then 
			Helper.worked_or_throw_error(false, "No tRNA_#{@anticodon} found.")
		else
			# save sequence as fasta
			@f_seq = File.join(@base_path, "trna.fasta")
			fasta = Helper::Sequence.str2fasta("tRNA, anticodon CAG", seq)
			File.open(@f_seq, 'w') {|f| f.write(fasta)}
		end
		return seq, struct, anticodon_pos, score
	# rescue 
	# 	throw :problem, "tRNA scan failed."
	end

	# parse seq_id and e-values out of all blasthits
	def parse_blasthits
		seq_ids, e_values, scores = [], [], []
		file_lines = File.read(@f_blasthits)
		score_lines = file_lines.match(/Sequences producing significant alignments[^\n]*\n\n(.*)\n\n\nQuery_1/m)[1]

		score_lines.lines.to_a.each do |line|
			line.chomp!
			fields = line.split(/\s+/)
			seq_ids.push( fields[0].sub(/^[^\|]+\|/,"") ) # fasta header, everything after first "|"
			scores.push( fields[1] ) # bit score
			e_values.push( fields[2] ) # e-value
		end

		return seq_ids, e_values, scores

	# rescue 
	# 	throw :problem, "tRNA scan failed."
	end

	# parse blast output and writes alignment to file
	def parse_blast_alignment
		f_alignment = File.join(@base_path, "trna-blast-aligned.fasta")

		fasta = ""
		file_lines = File.read(@f_blasthits)
		alignment_lines = file_lines.match(/(Query_1.*)\n\nLambda/m)[1]

		fasta_parts = {}
		# each seq might be splitted over multiple lines
		# sequences might contain blanks
		start_query_seq = nil
		stop_query_seq = nil
		alignment_lines.lines.to_a.each do |line|
			line.chomp!
			next if line == ""
			if line.include?("Query") then
				fields = line.split(/\s+/) 

				# seq_positions = line.match(/(\d+)\s+[-A-X]+\s+(\d+)/)
				# start_query_seq = seq_positions[1].to_i
				# stop_query_seq = seq_positions[2].to_i

				# range_start = fields[1].to_i
				# range_stop = fields[3].to_i
				# start_query_seq = line.upcase.index(@tRNA_seq[range_start-1..range_stop-1])
				# stop_query_seq = line.upcase.rindex(@tRNA_seq[range_stop-1])
			end
			fields = line.split(/\s+/)
			header = fields[0]

			seq = line.match(/\d+\s+([-A-X]+)\s+\d+/)[1]
			# seq = line[start_query_seq..stop_query_seq]
			seq = seq.gsub(" ", "-")
			if fasta_parts.has_key?(header) then 
				fasta_parts[header] += seq 
			else
				fasta_parts[header] = seq
			end
		end

		# concatenate sequences to fasta
		fasta_parts.each do |header, seq|
			this_fasta = Helper::Sequence.str2fasta(header, seq)
			fasta += this_fasta + "\n"
		end

		File.open(f_alignment, 'w') { |f| f.write(fasta) }
		return f_alignment 
	# rescue 
	# 	throw :problem, "tRNA scan failed."
	end

	# map blast hits to the closes translation
	def parse_and_analyse_blasthits
		seq_ids, e_values, scores = parse_blasthits
		res = {}

		seq_ids.each_with_index do |seq_id, ind|
			res[seq_id] = {
				e_val: e_values[ind],
				score: scores[ind],
				transl: convert_seqid_to_tRNA_identity(seq_id),
				is_significant: true
			}
			if res[seq_id][:transl] == "?" then 
				res[seq_id][:is_significant] = false
			end
		end

		return res
	# rescue 
	# 	throw :problem, "tRNA scan failed."
	end 

	# parse translation of fasta header from reference data
	def convert_seqid_to_tRNA_identity(seq_id)
		transl = "?"
		if seq_id.match(/Ser/i) then 
			transl = "S"
		elsif seq_id.match(/Leu/i)
			transl = "L"
		end

		return transl
	end

	def save(results)
		results[:seq] = @tRNA_seq if ! @tRNA_seq.blank?
		results[:struct] = @tRNA_sec_struct if ! @tRNA_sec_struct.blank?

		results[:anticodon] = @anticodon if ! @anticodon.blank?
		results[:anticodon_pos] = @anticodon_pos if ! @anticodon_pos.blank?
		results[:score] = @score

		results[:blast_hits] = @blast_results if @blast_results.any?

		results[:message] = @err_msg if ! @err_msg.blank?

		return results
	end

end