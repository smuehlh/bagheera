module ProgCall
	extend self

	attr_accessor :augustus_config_path, :augustus_species, :blosum62_path, :genome_db, :blast_filtering

	@augustus_species = "candida_albicans"
	@augustus_config_path = "/usr/local/bin/augustus/config/"

	@blast_filtering = "F"

	@blosum62_path = Rails.root.join('lib', 'blosum62').to_s
	@genome_db = ""


	def create_blast_db(f_in)
		# -i [INPUT: GENOME FILE] -p [PROTEIN] FALSE -n [DB NAME] -o [CREATE INDEX OVER SEQID] TRUE
		is_success = system(FORMATDB, "-i", f_in, "-p", "F", "-n", @genome_db, "-o", "T")
		# no output needed, so system is sufficient
		Helper.worked_or_die(is_success, "Cannot create BLAST database from genome data.")
	end

	def blast(f_out, f_ref)
		# -p [PROGRAM] protein query against nt-db -d [DATABASE] 
		# -i [FILE] -m8 [OUTPUT FORMAT] -F [FILTERING] -s [SMITH-WATERMAN ALIGNMENT] T 
		stdin, stdout_err, wait_thr = Open3.popen2e(BLASTALL, 
			"-p", "tblastn", "-d", @genome_db, "-i", f_ref, "-m8", "-F", @blast_filtering, "-s", "T")
		stdin.close
		output = stdout_err.read
		stdout_err.close

		if ! wait_thr.value.success? || output.include?("ERROR") || output.empty? then
			return false
		end

		File.open(f_out, 'w') {|f| f.write(output)}
		return true 
	end

	def fastacmd(f_out, seq_id, start, stop, strand)
		# -d [DB] -p [PROTEIN] false -s [SEARCH STRING] -L [START,STOP]
		# add 1000 nucleotides to start/stop 
		stdin, stdout_err, wait_thr = Open3.popen2e(FASTACMD, "-d", @genome_db, "-p", "F", "-s", seq_id, 
		"-L", "#{start},#{stop}", "-S", strand.to_s, "-o", f_out)
		stdin.close
		output = stdout_err.read
		stdout_err.close
		if ! wait_thr.value.success? then
			return false
		end
		return true,output
	end

	def augustus(f_out, f_in, *prfl_file)
		# --species [REFERENCE SPEC] QUERY --genemodel=exactlyone [predict exactly one complete gene] --codingseq=on [output also coding sequence]
		# --strand=forward [use only forward strand for gene prediction; this works as fastacmd always translates into forward strand]
		# --proteinprofile profilefile [use protein profile about alignment as foundation for gene prediciton]
		# redirection: add stderr to stdout (screen)
		if prfl_file.any? then
			stdin, stdout_err, wait_thr = Open3.popen2e(AUGUSTUS, "--AUGUSTUS_CONFIG_PATH=#{@augustus_config_path}", 
				"--species=#{@augustus_species}", "--genemodel=exactlyone", "--strand=forward", "--codingseq=on", 
				"--proteinprofile=#{prfl_file[0]}", f_in)
		else
			stdin, stdout_err, wait_thr = Open3.popen2e(AUGUSTUS, "--AUGUSTUS_CONFIG_PATH=#{@augustus_config_path}", 
				"--species=#{@augustus_species}", "--genemodel=exactlyone", "--strand=forward", "--codingseq=on", f_in)
		end
		stdin.close
		output = stdout_err.read
		stdout_err.close
		if (! wait_thr.value.success?) || output.include?("ERROR") || ! output.include?("coding sequence") then
			return false, wait_thr.value.termsig
		end
		File.open(f_out, 'w') {|f| f.write(output)}
		return true, 0
	end

	def mafft(f_out, f_in_pred, f_in_ref)
		# --addfragments: add non-full length sequence to existing MSA
		# --amino: input is protein
		# --anysymbol: replace unusal symbols by "X"
		# --quiet: output only alignment
		stdin, stdout_err, wait_thr = Open3.popen2e(MAFFT, "--addfragments", f_in_pred, "--amino", "--anysymbol", "--quiet", f_in_ref)
		stdin.close
		output = stdout_err.read
		stdout_err.close
		if ! wait_thr.value.success? then
			return false, ""
		end
		File.open(f_out, 'w') {|f| f.write(output)}
		return true
	end

	def pairalign(f_out, f_in, conf, meth)
		# -c tttt: initialize first row & column with zeros, search last row & column for maximum
		# 	=> end gap free alignment
		# --matrx matrix file
		# implicit options: gap open penalty: -11, gap extension penalty: -1
		# 					protein sequences
		# don't use system, as pair_align is quite verbose and has an easy parsable success-output
		stdin, stdout_err, wait_thr = Open3.popen2e(PAIR_ALIGN, "--seq", f_in, "--matrix", @blosum62_path, 
			"--config", conf, "--method", meth, "--outfile", f_out)
		stdin.close
		output = stdout_err.read
		stdout_err.close
		if ! wait_thr.value.success? || ! output.include?("Alignment score") then
			return false
		end
		return true
	end

	def gblocks(f_in)
		stdout = IO.popen([GBLOCKS, f_in, "-b5=h", "-p=n"])
		output = stdout.read
		stdout.close
		# return value is always false, so parse output to find out if it was successful
		if output.include?("selected block(s)") then
			return true
		else
			return false
		end
	end

	def fasttree(f_out, f_in)
		is_success = system(FastTree, "-out", f_out, f_in)
	end

	# def config(configfile)

	# 	@augustus_species = session[:augustus][:species]
	# 	@blast_filtering = session[:blast][:use_low_comp_filter]

	# end
end