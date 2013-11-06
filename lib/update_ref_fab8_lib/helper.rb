module Helper
	extend self

	# general error handling methods
	def abort(msg="Fatal: An error occured")
		$stderr.puts msg
		exit 1
	end

	def worked_or_die(is_success, msg)
		if ! is_success then
			abort "Fatal: #{msg}"
		end
	end

	# def worked_or_throw_error(is_success, msg)
	# 	if ! is_success then
	# 		throw :problem, msg
	# 	end
	# 	is_success
	# end

	# general file handling methods
	def file_exist_or_die(path)
		if ! FileTest.file?(path) then
			abort "Fatal: File #{path} does not exist."
		end
	end

	def does_file_exist(path)
		if ! FileTest.file?(path) then
			return false
		end
		true
	end

	# def dir_exist_or_die(path)
	def does_dir_exist(path)
		if ! FileTest.directory?(path) then
			# abort "Fatal: Directory #{path} does not exist."
			return false
		end
		true
	end

	def mkdir_or_die(path)
		if ! FileTest.directory?(path) then
			begin
				Dir.mkdir(path)
			rescue
				abort "Fatal: Cannot create directory #{path}"
			end
		end
	end

	def move_or_die(f_scr, f_dest)
		FileUtils.mv(f_scr, f_dest)
	rescue Errno::ENOENT => exc
		abort "Fatal: Cannot move files"
	end

	def del_file_or_dir(path)
		FileUtils.rm_rf path # this way, it can handle both string and array
	rescue Errno::ENOENT => exc
		warn "Cannot delete #{path}"
	end

	def run_mafft(f_in, f_frag, is_store_to_f_in=false)
		stdin, stdout_err, wait_thr = "", "", ""
		if f_frag.empty? then
			# align normally
			stdin, stdout_err, wait_thr = Open3.popen2e("/usr/local/bin/mafft", 
				"--amino", "--anysymbol", "--quiet", f_in)
		else
			# add fragment
			stdin, stdout_err, wait_thr = Open3.popen2e("/usr/local/bin/mafft", 
				"--addfragments", f_frag, "--amino", "--anysymbol", "--quiet", f_in )
		end
		stdin.close
		output = stdout_err.read
		stdout_err.close
		if ! wait_thr.value.success? then
			return false
		end

		if is_store_to_f_in then
			# save mafft - output to f_in
			File.open(f_in, 'w') {|f| f.write(output)}
		end

		true

	end

	def calc_protein_profile(f_in, f_out)
		# test if the enviromental variable is set
		if ! ENV.has_key?("AUGUSTUS_CONFIG_PATH") then
			# not set, so set it now
			ENV["AUGUSTUS_CONFIG_PATH"] = "/usr/local/bin/augustus/config/"
		end

		# calculate profiles
		output = system "/usr/local/bin/augustus/scripts/msa2prfl.pl", f_in, :out => f_out, :err => '/tmp/cug/err.log'
	end

	module Sequence
		extend self

		def fasta2str(fasta)
			headers = []
			seqs = []
			# get every record: everything between two ">"
			fasta.scan(/^>([^>]*)/m).flatten.each do |rec|
				rec.chomp!
				nl = rec.index("\n") # first match: separating header from seq
				header = rec[0..nl-1]
				seq = rec[nl+1..-1]
				seq.gsub!(/\r?\n/,'') # removing all remaining \n from seq 
				headers.push(header)
				seqs.push(seq)
			end
			return headers, seqs
		end


		def str2fasta(header, seq)
			fasta = header.include?(">") ? header << "\n" : ">" << header << "\n"
			fasta += seq
		end

		def save_alignment(file, fasta_hash)
			fasta = fasta_hash.map{|k,v| ">#{k}\n#{v}"}.join("\n") + "\n"
			File.open(file, 'w'){ |f| f.write(fasta) }
		end
	end
end 