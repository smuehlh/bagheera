class CronjobController < ApplicationController
# debug me by calling "script/rails runner 'CronjobController.prepare_ref_data'" on command line and including a "debugger" anywhere
	def self.prepare_ref_data

		# wait for wget to be finished, but don't wait forever! (max. 10 minutes)
		start = Time.now
		puts "Starting execution: prepare reference data. Wait until they are fetched!"
		begin
			fin = Time.now
			puts fin
			sleep(30) 
		end while (! FileTest.file?(BASE_PATH + "new_" + REF_DATA) || fin - start >= 600)
		if fin - start >= 600 then
			puts "ERROR with wget. Aborting."
			exit
		else
			puts "wget finished: #{fin}"
			# waiting another 30 seconds to make sure the entire file is written
			sleep(30)
		end

		# create directory for newly created reference data
		path_new_data = BASE_PATH + "new/"
		puts "Saving all new data to temporaray file #{path_new_data}"
		if File.exists?(path_new_data) then
			FileUtils.rm_r(Dir.glob(path_new_data + "*"), :verbose => true)
		else
			FileUtils.mkdir(path_new_data, :mode => 0775)
		end

		puts ("Loading reference data file ... ")
		# load reference data
		data, errors = load_ref_data
		puts errors.join(" ")
puts data.keys
		if errors.empty? then
			puts "Successfully loaded reference data"
			puts "Separating myosin, actin and kinesin by class"

			# separate actin, myosin and kinesin by class
			myo_data = separate_by_class(data, "Myosin heavy chain")
puts "Myo done"
			actin_data = separate_by_class(data, "Actin related protein")
puts "Actin done"
			kinesin_data = separate_by_class(data, "Kinesin")
puts "Kin done"
			# delete unseparated data 
			puts "Deleting unseparated data"

			data.delete("Myosin heavy chain")
			data.delete("Actin related protein")
			data.delete("Kinesin")
			# add separated (do this after deletion of unseparated !!! as old key is used as "default class")
			puts "Add separated data to reference data"

			data.merge!(myo_data)
			data.merge!(actin_data)
			data.merge!(kinesin_data)

			# delete kinesin and myosin orphans, and myosin class 17 (extremly unusual for saccharomycetes!)
			data.delete("Myosin heavy chain")
			data.delete("Myosin heavy chain Class 17")
			data.delete("Kinesin")

			# save updated reference data
			puts "Saving updated reference data to file"
			fh = File.new(path_new_data + REF_DATA, "w")
			fh.puts JSON.dump(data)
			fh.close

			# ensure all alignments are of same length
			puts "Ensure reference data are aligned"
			puts "Writing alignments ... "
			data.each do |prot, prot_data|

				data_file = path_new_data + prot.gsub(" ", "-").downcase + ".fasta"
				fasta = prot_data["alignment"]
				File.open(data_file, 'w'){|f| f.write(fasta)}
				puts "\t file #{data_file}"

				# check and correct length of each aligned sequence
				is_true = ensure_length(data_file)
				if ! is_true then
					puts "\t #{prot}: Not all sequences were of same lenght. Now they are."
				end
				# check if mafft accepts them as aligned
				is_true = ensure_mafft_is_fine(data_file)
				if ! is_true then
					puts "\t #{prot}: Mafft thinks, there are not aligned. So align them with mafft ... "

					# align them with mafft
					is_true = run_mafft(data_file)
					if is_true then 
						puts " done."
					else
						puts " ... Could not align."
					end
				end 
				
				# calculate profile
				puts "Generating protein profiles ... "
				is_true = calc_protein_profiles(data_file)

			end # data.each
			puts "Finished."
			# set correct permissions

			begin
				# puts "Moving new file to place ... "
				# FileUtils.mv(BASE_PATH + "new_" + REF_DATA, BASE_PATH + REF_DATA)
				puts "Moving new reference data to place ..."
				FileUtils.rm Dir.glob(BASE_PATH + "*.fasta")
				FileUtils.mv Dir.glob(path_new_data + "*"), BASE_PATH, :verbose => true
				FileUtils.rm_rf path_new_data, :secure => true, :verbose => true
				FileUtils.rm_rf "/tmp/cug/new_alignment_gene_structure.json", :verbose => true
				puts "done."
			rescue => e
				puts "Could not move new files to place:"
				puts e
			end	
		else
			puts "Could not load/ prepare reference data"
		end # if errors.empty?
	end 

	# delete old data in /tmp/cug/ and /tmp/cymobase_alignment_cug*
	# data from old predictions and show_alignment requests
	def delete_old_data(days = 1)
		begin
			Dir.glob(BASE_PATH + "**/") do |dir|
				if dir == BASE_PATH then
					next
				elsif File.mtime(dir) <= days.day.ago
					FileUtils.rm_rf(dir, :secure=>true)
				end
			end

			# delete alignment file for show-alignment function
			FileUtils.rm Dir.glob(Dir::tmpdir + '/cymobase_alignment_cug*')
		rescue => e
			puts "Could not delete old data: #{e}"
		end
	end

	# calculate protein profiles for reference alignments (needed by augustus to enhance gene prediction)
	# @param file [String] reference alignment file
	# @return [Boolean] Status of profile calculation
	def self.calc_protein_profiles(file)
		# test if the enviromental variable is set
		if ! ENV.has_key?("AUGUSTUS_CONFIG_PATH") then
			# not set, so set it now
			ENV["AUGUSTUS_CONFIG_PATH"] = "/usr/local/bin/augustus/config/"
		end
		prfl_file = file.sub("fasta", "prfl")
		# calculate profiles
		output = system "/usr/local/bin/augustus/scripts/msa2prfl.pl", file, :out => prfl_file, :err => '/tmp/cug/err.log'
		return true
	end

	# load reference data from file alignment_gene_structure.json
	# @return [Hash] Reference data
	# @return [Array] Errors occured during file load
	def self.load_ref_data
		path = BASE_PATH + "new_" + REF_DATA
		if ! FileTest.file?(path) then
			puts "Fatal error: "
			errors = ["Fatal error.", "Cannot load reference data. Please contact us!"]
			return false, errors
		end
		return JSON.load(File.read(path)), []
	end

	# separate reference data for action, myosin and kinesin by class
	# @param data [Hash] reference data
	# @param prot [String] key in data - hash, those data to separate
	def self.separate_by_class(data, prot)
		new_data = {}
		genes = data[prot]["genes"].keys
		old_alignment = data[prot]["alignment"]
		headers, seqs = fasta2str(old_alignment)

begin
		genes.each do |name|
			if prot.include?("Myosin") then
				name =~ /[a-zA-Z_]+(Myo)?([0-9]+|Mhc)/
				klass = $2 if $2
			elsif prot.include?("Actin") 
				name =~ /[a-zA-Z_]+(Arp)([0-9]+)/
				klass = $2 if $2
			else
				name =~ /[a-zA-Z_]+(Kinesin)([0-9]+)/
				klass = $2 if $2
			end
		  
			if klass then
				# could extract class name $2
				key = prot + " Class " + $2.to_s # convert to key name
			else
				# could not extract a class name, simply use protein name as class name
				key = prot
			end

			if ! new_data.has_key?(key) then
				# initialize hash for this class
				new_data[key] = {}
				new_data[key]["genes"] = {}
				new_data[key]["alignment"] = ""
			end

			# add gene structure and aligned sequence to new data hash
			new_data[key]["genes"][name] = data[prot]["genes"][name]

			if headers.include?(name) then
				ind = headers.index(name)
				new_data[key]["alignment"] += str2fasta(name, seqs[ind], true) << "\n" # true: no line break after 80 chars
				# no else needed: even if no sequence exists for this key, the other functions can handle missing data
			end
		end
rescue => e
	puts e
end
		
	
		return new_data
	end

	def self.ensure_length(file)
		ret_value = true
		lines = IO.readlines(file)#.map(&:chomp) # chomps every line
		seq_lines = lines.reject{|line| line.match(/^>/)}
		len = seq_lines.collect {|line| line.size}

		if ! (len.min == len.max) then
			ret_value = false
			# add gaps to end of each line shorter than longest line
			max = len.max 
			lines.map! do |line|
				if line.match(/^>/) then
					line
				elsif line.size == max then
					line
				else
					line.chomp!
					line = line + "-" * (max - 1 - line.size) + "\n" # -1 as in max the "\n" is counted
				end
			end
		end

		File.open(file, 'w') {|f| f.write(lines.join(""))}
		return ret_value
	end

	def self.ensure_mafft_is_fine(file)
		# use first seq to align with, if it works	
		test_file = BASE_PATH + "test.fasta"
		headers, seqs = fasta2str(File.read(file))
		File.open(test_file, 'w') {|f| f.write(str2fasta(">firstseq", seqs[0], true))}

		stdin, stdout_err, wait_thr = Open3.popen2e("/usr/local/bin/mafft", "--addfragments", test_file, "--amino", "--anysymbol", file)
		stdin.close
		output = stdout_err.read
		stdout_err.close
		FileUtils.rm(test_file)
		if ! wait_thr.value.success? then
			return false
		end
		return true
	end

	def self.run_mafft(file)
		stdin, stdout_err, wait_thr = Open3.popen2e("/usr/local/bin/mafft", "--auto", "--amino", "--anysymbol", "--quiet", file)
		if ! wait_thr.value.success? then
			return false
		end
		stdin.close
		output = stdout_err.read
		stdout_err.close
		File.open(file, 'w') {|f| f.write(output)}
		return true
	end

	# convert header and sequence into fasta-format
	# @param header [String] fasta header
	# @param seq [String] sequence
	# @param no_split [Boolean] include line break each 80 chars?
	# @return [String] fasta formatted header and sequence
	def self.str2fasta(header, seq, no_split=false)
		fasta = header.include?(">") ? header << "\n" : ">" << header << "\n"
		if no_split then
			fasta += seq
		else
			fasta += seq.scan(/.{1,80}/).join("\n")
		end
		return fasta
	end

	# extract headers and sequences from a fasta-formatted string
	# can deal with multiple sequence fasta
	# @param fasta [String] fasta-formatted sequence (MSA also possible)
	# @return [Array] All fasta header found
	# @return [Array] All fasta sequences found
	def self.fasta2str(fasta)
		headers = []
		seqs = []
		# get every record: everything between two ">"
		fasta.scan(/^>([^>]*)/m).flatten.each do |rec|
			rec.chomp!
			nl = rec.index("\n") # first match: separating header from seq
			header = rec[0..nl-1]
			seq = rec[nl+1..-1]
			seq.gsub!(/\n/,'') # removing all remaining \n from seq 
			headers.push(header)
			seqs.push(seq)
		end
		return headers, seqs
	end

end