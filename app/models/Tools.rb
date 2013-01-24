class Tools

end

class NeedlemanWunsch
	attr_accessor :sequence, :blosum62

	# gap function: g(k) = Gap_open + Gap_extend * k
	Gap_open = -11
	Gap_extend = -1
	Inf = 1.0/0.0 # infinity
	
	def initialize(sequence)
	    @sequence 				= sequence 		# sequence to be compared
    	@blosum62 				= get_blosum62  # hash containing blosum62 matrix
    end

	def get_blosum62
		hash = {}
		file = File.read("#{Rails.root}/lib/blosum62")
		aas = nil
		rows = file.lines.to_a
		rows.each do |row|
			next if row.length <= 1
			next if row[0].chr == "#"

			unless aas then
				aas = row.split(" ").to_a
				next
			end

			cells_of_row = row.split(" ").to_a
			aa = cells_of_row[0]
			cells_of_row.delete_at(0)
			aas.each_with_index do |bb, i|
				hash[aa.to_s + bb.to_s] = cells_of_row[i]
			end
		end
		return hash
 	end

# new for affine linear gap costs
# FIXME
# initialization not working ...
 	def gotoh(reference, is_endgap_free)
 		cols = sequence.length + 1 
 		rows = reference.length + 1
 		score = Array.new(rows){Array.new(cols,0)}
 		insert = del = Array.new(rows){Array.new(cols,-Inf)} # insert: costs to insert a char in score; del: costs to delete char in score
 		# for each cell: use optimum from those values: 
		opt_in = opt_del = [-Inf, -Inf]
		opt_score = [-Inf, -Inf, -Inf]
		# return values
		ref = ''
		seq = ''
		alignment_score = 0

 		# compute score (and initialize first row and col of score)
 		# score[0] = Array.new(cols, gap_function(1)) # initialize first row
 		(rows).times do |i|
			# score[i][0] = gap_function(1) # initialize first column
 			(cols).times do |j|
 				if j > 0 then
 				puts "writing insertmat for i = #{i}, j = #{j}"			
	 				opt_in[0] = score[i][j-1] + gap_function(1)
 				# opt_in[1] = (j >= 2) ? insert[i][j-1] + Gap_extend : -Inf 
	 				opt_in[1] = insert[i][j-1] + Gap_extend #if (j >= 1) 
	 				insert[i][j] = opt_in.max
	 			end
				
				if i > 0 then
					puts "writing delmat for i = #{i}, j = #{j}"
					opt_del[0] = score[i-1][j] + gap_function(1)
				# opt_del[1] = (i >= 2) ? del[i-1][j] + Gap_extend : -Inf
					opt_del[1] = del[i-1][j] + Gap_extend #if (i >= 1)
					del[i][j] = opt_del.max
				end

				opt_score[0] = score[i-1][j-1] + aa_score(reference[i-1], sequence[j-1])
				opt_score[1] = insert[i][j]
				opt_score[2] = del[i][j]
				score[i][j] = opt_score.max
 			end
 		end

 		# backtracing
 		# 1) find starting position for trace back
		if is_endgap_free then
			# end gap free alignment: starting at pos with max. score
			max_row = score[rows-1][0..cols-1].each_with_index.max # first element: val, second: index
			max_col = score.inject([]) {|res, val| res << val.last}.each_with_index.max  
			if max_row[0] > max_col[0] then
				i = reference.length
				j = max_row[1]
				ref = '-' * (sequence.length - j)
				seq = sequence[j..sequence.length]
			else
				i = max_col[1]
				j = sequence.length
				seq = '-' * (reference.length - i) 
				ref = reference[i..reference.length]
			end
		else
			# "normal" alignment: starting always at lower right corner
			i = reference.length
			j = sequence.length
		end	
		# 2) actual backtracing
		max_score = score[i][j]

		while (i > 0 && j > 0) do
			this_score = score[i][j]
			this_del = del[i][j]
			this_insert = insert[i][j]

			score_diag = score[i-1][j-1] + aa_score(reference[i-1], sequence[j-1])
			score_left = score[i][j-1] + gap_function(1)
			score_up = score[i-1][j] + gap_function(1)

			use_ins = (this_insert == this_score) ? true : false
			use_del = (this_del == this_score) ? true : false

			ins_left = insert[i][j-1] + Gap_extend
			del_up = del[i-1][j] + Gap_extend
debugger
			if (this_score == score_diag) then
				ref = reference[i-1] + ref
				seq = sequence[j-1] + seq
				i -= 1
				j -= 1
				alignment_score += aa_score(reference[i-1], sequence[j-1])
			elsif ( this_insert == score_left && ! use_ins)
				ref = reference[i-1] + ref 
				seq = '-' + seq 
				i -= 1
				alignment_score += gap_function(1)
			elsif (this_insert == ins_left && use_ins)
				ref = reference[i-1] + ref 
				seq = '-' + seq 
				i -= 1
				alignment_score += Gap_extend
			elsif ( this_del == score_up && ! use_del)
				ref = '-' + ref 
				seq = sequence[j-1] + seq 
				j -= 1
				alignment_score += gap_function(1)
			elsif (this_del == del_up && use_del)
				ref = '-' + ref 
				seq = sequence[j-1] + seq 
				j -= 1
				alignment_score += Gap_extend
			else
				puts "Ups! Unknown case ..."
			end
		end

		while (i > 0) do 
			ref = reference[i-1] + ref 
			seq = '-' + seq 
			i -= 1
		end

		while (j > 0) do
			ref = '-' + ref 
			seq = sequence[j-1] + seq 
			j -= 1
		end

		return seq, ref, alignment_score/max_score.to_f
 	end

# new for gotoh
 	
 	def gap_function(gap_length)
		return (Gap_open + gap_length * Gap_extend)
	end


	# needleman wunsch implementation: compare @sequence with reference
	# optional: use free end-gaps
	# optional: use affine linear gap costs
	def needleman_wunsch(reference, is_endgap_free, is_const)
		rows = reference.length + 1
		cols = sequence.length + 1
		arr = Array.new(rows){Array.new(cols,0)}
		choice = [0, 0, 0]
		ref = ''
		seq = ''
		score, score_diag, score_up, score_left = nil

		reference.length.times do |i|
			sequence.length.times do |j|
				choice[0] = arr[i][j] + aa_score(reference[i], sequence[j])
				choice[1] = arr[i][j+1] + gap_cost(arr[i][j+1], arr[i][j], is_const)
				choice[2] = arr[i+1][j] + gap_cost(arr[i+1][j], arr[i][j], is_const)
				arr[i+1][j+1] = choice.max
			end
		end

		if is_endgap_free then
			# end gap free alignment: starting at pos with max. score
			max_row = arr[rows-1][0..cols-1].each_with_index.max # first element: val, second: index
			max_col = arr.inject([]) {|res, val| res << val.last}.each_with_index.max  
			if max_row[0] > max_col[0] then
				i = reference.length
				j = max_row[1]
				ref = '-' * (sequence.length - max_row[1])
				seq = sequence[max_row[1]..sequence.length]
			else
				i = max_col[1]
				j = sequence.length
				seq = '-' * (reference.length - max_col[1]) 
				ref = reference[max_col[1]..reference.length]
			end
		else
			# "normal" needleman wunsch alignment: starting always at lower right corner
			i = reference.length
			j = sequence.length
		end		
 
		max_score = arr[i][j]
		alignment_score = 0
		while (i > 0 && j > 0) do
			score = arr[i][j]
			score_diag = arr[i-1][j-1]
			score_up = arr[i][j-1]
			score_left = arr[i-1][j]

			if ( score == (score_diag + aa_score(reference[i-1], sequence[j-1])) ) then
				ref = reference[i-1] + ref
				seq = sequence[j-1] + seq 
				i -= 1
				j -= 1
				alignment_score += aa_score(reference[i-1], sequence[j-1])
			elsif ( score == (score_left + gap_cost(score_left, arr[i-2][j], is_const)) )
				ref = reference[i-1] + ref 
				seq = '-' + seq 
				i -= 1
				alignment_score += gap_cost(score_left, arr[i-2][j], is_const)
			elsif ( score == (score_up + gap_cost(score_up, arr[i][j-2], is_const)) )
				ref = '-' + ref 
				seq = sequence[j-1] + seq 
				j -= 1
				alignment_score += gap_cost(score_up, attr_accessor[i][j-2], is_const)
			else 
				puts "ERROR in Needleman-Wunsch"
				exit
			end
		end

		while (i > 0) do 
			ref = reference[i-1] + ref 
			seq = '-' + seq 
			i -= 1
		end

		while (j > 0) do
			ref = '-' + ref 
			seq = sequence[j-1] + seq 
			j -= 1
		end
		return seq, ref, alignment_score/max_score.abs.to_f
	end

	def gap_cost(score_pos, score_last_pos, is_const)
		# constant gap costs
		gap = -5
		return gap if is_const 
		# use affine linear gap costs
		gap_start = -7
		gap_extend = -1 
		if (score_pos == score_last_pos + gap_extend || score_pos == score_last_pos + gap_start) then
			# extend existing gap
			return gap_extend
		else
			# start a new gap
			return gap_start
		end
	end

	def aa_score(a, b)
		blosum62[(a.chr + b.chr).upcase].to_i
	end
end # class NeedlemanWunsch