module ProgCall
	extend self

	attr_accessor :augustus_config_path, :augustus_species, :blosum62_path, :genome_db, :blast_filtering

	@augustus_species = "candida_albicans"
	@augustus_config_path = "/usr/local/bin/augustus/config/"

	# @blast_filtering = "F"
	@blast_filtering = "no"
	@blast_filtering_masking = "false"

	@blosum62_path = Rails.root.join('lib', 'blosum62').to_s
	@genome_db = ""

	@is_tRNA_scan_general_model = true

	def create_blast_db(f_in, genome_db=@genome_db)

		path_to_program = File.join( New_blast_path, "makeblastdb")
		# -dbtype nucl [nucleotide] -out [DB NAME] -in [INPUT:FASTA] -parse_seqids [create index over seqid]
		is_success = system(path_to_program, 
			"-dbtype", "nucl", "-out", genome_db, "-in", f_in, "-parse_seqids")

		# no output needed, so system is sufficient
		Helper.worked_or_die(is_success, "Cannot create BLAST database from genome data.")
	end

	def blast(f_out, f_ref, options={})
		defaults = {
			:genome_db => @genome_db,
			:program => "tblastn",
			:output_format => "6" 
		}
		options = defaults.merge(options)

		path_to_program = File.join( New_blast_path, options[:program] )
		# -db [DB NAME] -query [INPUT FILE] -outfmt 6 [output in tabular format] -seg [yes|no] [filter query]
		# 	-soft_masking [true|false] [filter locations as soft mask]
		args_for_program = ["-db", options[:genome_db], "-query", f_ref, 
			"-outfmt", options[:output_format], "-soft_masking", @blast_filtering_masking]

		if options[:program] == "tblastn" then 
			# option "-use_sw_tback" only supported by tblastn!
			# FIXME - it does not work!!!
			# args_for_program.push("-use_sw_tback") # FIXME - causes blast to crash

			args_for_program.push("-seg", @blast_filtering)
		end

		stdin, stdout_err, wait_thr = Open3.popen2e(path_to_program, *args_for_program) 

		stdin.close
		output = stdout_err.read
		stdout_err.close

		if ! wait_thr.value.success? || output.include?("ERROR") || output.empty? then
			return false
		end

		File.open(f_out, 'w') {|f| f.write(output)}
		return true 
	end

	def fastacmd(f_out, seq_id, start, stop, strand, genome_db=@genome_db)

		path_to_program = File.join( New_blast_path, "blastdbcmd" )
		# -db [DB NAME] -dbtype nucl [nucleotide] -entry [SEARCH STRING] -out [OUT FILE] -range [start-stop in fasta] -strand [DNA strand]
		stdin, stdout_err, wait_thr = Open3.popen2e(path_to_program, 
			"-db", genome_db, "-dbtype", "nucl", "-entry", seq_id, "-out", f_out, 
			"-range", "#{start}-#{stop}", "-strand", strand)

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

	def mafft(f_out, f_in_pred, f_in_ref, options={})
		# --addfragments: add non-full length sequence to existing MSA
		# --amino: input is protein
		# --nuc: is nucleotide
		# --anysymbol: replace unusal symbols by "X"
		# --quiet: output only alignment

		defaults = {
			:input => "amino",
		}
		options = defaults.merge(options)

		stdin, stdout_err, wait_thr = Open3.popen2e(MAFFT, "--addfragments", f_in_pred, 
			"--#{options[:input]}", "--anysymbol", "--quiet", 
			f_in_ref)
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

	def trnascan(f_in)
		# -D: no search for pseudogenes, makes search faster
		# -f#: output sequence and its secondary structure
		# -Q: do not prompt before overwriting an existing file
		
		if @is_tRNA_scan_general_model then 
			# -G: general model, as opposed to a specific one for eukaryotes
			is_success = system(TRNASCAN, "-G", "-D", "-f#", "-Q", f_in )
		else
			# eukaryotic model
			is_success = system(TRNASCAN, "-D", "-f#", "-Q", f_in )
		end
		return is_success
	end

	# def config(configfile)

	# 	@augustus_species = session[:augustus][:species]
	# 	@blast_filtering = session[:blast][:use_low_comp_filter]

	# end

	def set_blast_filtering(is_filtering) 
		if is_filtering then 
			# @blast_filtering = "yes" #"T"
			@blast_filtering = "yes"
			@blast_filtering_masking = "true"
		else
			@blast_filtering = "no" # "F"
			@blast_filtering_masking = "false"
		end
	end

	def set_tRNAscan_model(is_general_model)
		@is_tRNA_scan_general_model = is_general_model
	end

	# how to call the old blast version
	# formatdb
		# is_success = system(FORMATDB, "-i", f_in, "-p", "F", "-n", genome_db, "-o", "T")
	# blast
		# stdin, stdout_err, wait_thr = Open3.popen2e(BLASTALL, 
		# "-p", options[:program], "-d", options[:genome_db], "-i", f_ref, "-m8", "-F", @blast_filtering, "-s", "T")
		# caution: @blast_filtering must be "F" or "m S"; no @blast_filtering_masking
		# c
	# fastacmd
		# stdin, stdout_err, wait_thr = Open3.popen2e(FASTACMD, "-d", genome_db, "-p", "F", "-s", seq_id, 
		# 	"-L", "#{start},#{stop}", "-S", strand.to_s, "-o", f_out)
		# caution: strand must be "1" or "2"
end