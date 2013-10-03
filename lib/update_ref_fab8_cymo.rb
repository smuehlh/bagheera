#!/usr/local/bin/ruby

### Update Bagheeras reference data on fab8
### Called monthly in crontab

# 00 04 1 * * ruby /fab8/server/bagheera/lib/update_ref_fab8_cymo.rb
# require 'ruby-debug'
require 'timeout'
require 'fileutils'
require 'json'
require 'open3'

Base_path = "/tmp/cug"
Without_ca_b = "without_ca_b/"
New_ref_data = "new_alignment_gene_structure.json"

Del_unusual_prots = ["Calcineurin", "Calmodulin", "Centrin", "Frequenin", 
	"Dynein Light Intermediate Chain",
	"Myosin heavy chain", "Myosin essential light chain", "Myosin regulatory light chain", "Myosin heavy chain Class 17",
	"Kinesin", "Kinesin Class 4", "Kinesin Class 16"]
Split_prot_families = ["Myosin heavy chain", "Actin related protein", "Kinesin", "Tubulin", "Capping Protein"]
Split_prot_abbrs = { "Myosin heavy chain" => "Myo", 
	"Actin related protein" => "Arp", 
	"Kinesin" => "Kinesin", 
	"Tubulin" => "Tub", 
	"Capping Protein" => "CAP" }

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

	def worked_or_throw_error(is_success, msg)
		if ! is_success then
			throw :problem, msg
		end
		is_success
	end

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

	def load_ref_data(path)
		file_exist_or_die(path)
		return JSON.load(File.read(path)) 
	end

	def save_ref_data(path, data)
		fh = File.new(path, "w")
		fh.puts JSON.dump(data)
		fh.close
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

	end
end 

class ProteinFamily
	attr_reader :prot, :prot_filesave_name, 
		:ref_alignment, :ref_genes

	def initialize(prot, data)
		@prot = prot 
		@prot_filesave_name = get_filesave_name 
		# @ref_prfl_file = get_fname("prfl")
		@ref_alignment, @ref_genes = get_refdata(data)
	end

	def get_filesave_name
		@prot.gsub(" ", "-").downcase
	end

	def get_refdata(data)
begin
		seqs = Hash[*data["alignment"].split("\n")]
rescue => exc
	Helper.abort("Alignment of protein #{@prot} is corrupt.")
end
		seqs.keys.each {|k| seqs[ k.sub(">", "") ] = seqs.delete(k)}
		genes = data["genes"]
		return seqs, genes
	end

	# separate by class & delete unseparated data & store back in ref_data
	def split_into_classes(ref_data)

		result = {} # store them back to ref_data, keys: class identifier
	
		# first: separate this obj into classes
		keys = @ref_alignment.keys & @ref_genes.keys
		keys.each do |name|

			prot_abbr = Split_prot_abbrs[@prot]

			klass = ""
			new_prot_key = ""

			matches = name.match(/#{prot_abbr}([0-9]+|Mhc)/)

			# Myosin protein has to abbreviations: Mhc for class 2 and Myo for all other classes
			if !  matches && name.include?("Mhc") then
				matches = ["", "Mhc"]
			end

			if matches && matches[1] then
				klass = matches[1]
				new_prot_key = @prot + " Class " + klass
			else
				new_prot_key = @prot
			end

			if ! result.has_key?(new_prot_key) then
				result[new_prot_key] = {}
				result[new_prot_key]["alignment"] = ""
				result[new_prot_key]["genes"] = {}
			end

			# save genes and alignments
			if @ref_genes.has_key?(name)
				result[new_prot_key]["genes"][name] = @ref_genes[name]
			end

			if @ref_alignment.has_key?(name)
				result[new_prot_key]["alignment"] += Helper::Sequence.str2fasta(name,@ref_alignment[name])
				result[new_prot_key]["alignment"] += "\n"
			end

		end

		# second: delete unseparated data (that is basiscally which lives in this obj) from reference data
		ref_data.delete(@prot)

		# third: save separated data into reference data
		ref_data.merge!(result)

	end

	def ensure_same_length

		seq_lines = @ref_alignment.values
		len = seq_lines.collect {|line| line.size}

		if ! (len.min == len.max) then
			# add gaps to end of each line shorter than longest line
			max = len.max 
			@ref_alignment.each do |key, val|
				
				if val.size != max then
					@ref_alignment[key] = val + "-" * (max - val.size)
				end

			end

		end

	end

	# method is based on /fab8/db_scripts/alignment/edit.rb
	def remove_common_gaps
  
		seqs = @ref_alignment.values
		max_len = seqs.collect do |s| s.length end.max

		#look for common gaps
		gaps = []
		max_len.times do |pos|
			if (seqs.select do |s| (!s[pos] || s[pos].chr == "-") end.length == seqs.length) then
				gaps << pos
			end
		end
		gaps = gaps.reverse

		#delete gaps
		@ref_alignment.each_pair do |k, v|
			v = v.split("")
			gaps.each do |pos| v.delete_at(pos) end
			v = v.join("")
			@ref_alignment[k] = v
		end
	end

	def save_alignment(file)
		fasta = ""
		@ref_alignment.each do |key, val|
			fasta += Helper::Sequence.str2fasta(key, val)
			fasta += "\n"
		end
		File.open(file, 'w'){ |f| f.write(fasta) }
	end

	def ensure_mafft_is_fine(file, test_file)
		# try to add first sequence again
		File.open(test_file, 'w') { |f| f.write( Helper::Sequence.str2fasta(">test",@ref_alignment.values[0]) ) }
		is_success = Helper.run_mafft(file, test_file)
		Helper.del_file_or_dir(test_file)

		# if this is not successful, mafft does not recognize file as aligned
		if ! is_success then
			# align reference sequences with mafft
			Helper.run_mafft(file,"",true) # true: save output
			# Helper.worked_or_throw_error(is_success, "Mafft failed: #{file}")
		end
	end
end


def callCymoAPI 
	wget_path = `which wget`.chomp
	system(wget_path, "--spider", "http://fab8:2001/api_cug_alignment/all", "-o", File.join(Base_path, "cron.log"))
end

# "main script"

# fetch reference data from cymobase

ref_data_path = File.join(Base_path,New_ref_data) # File loaded by wget

new_data_path = File.join(Base_path,"new") # Folder containing new files
new_ref_data_path = File.join(new_data_path, New_ref_data.sub("new_","")) # Wget-file in new location

new_data_path_no_cab = File.join(new_data_path, Without_ca_b) # Folder containing new files for example
new_ref_data_path_no_cab = File.join(new_data_path_no_cab, New_ref_data.sub("new_","")) # Wget-file for example

Helper.del_file_or_dir(ref_data_path)
Helper.del_file_or_dir(new_ref_data_path)
Helper.del_file_or_dir(new_data_path_no_cab)

max_secs = 60*10 # wait max. 10 minutes for cymoapi
puts "Start wget: #{Time.now}"
begin
	status = Timeout::timeout(max_secs) { callCymoAPI }
rescue Timeout::Error => exc
	$stderr.puts exc
	exit 1
end
Helper.worked_or_die(status && Helper.does_file_exist(ref_data_path), "wget failed")
puts "Finished wget: #{Time.now}"

# prepare reference data

Helper.mkdir_or_die(new_data_path) if ! Helper.does_dir_exist(new_data_path)
Helper.move_or_die(ref_data_path,new_ref_data_path)

Helper.mkdir_or_die(new_data_path_no_cab) if ! Helper.does_dir_exist(new_data_path_no_cab)

ref_data = Helper.load_ref_data(new_ref_data_path)

# Important: order does matter! first split into classes, than delete unusal proteins
# separate by class & delete unseparated data & store back in all_prot_data
Split_prot_families.each do |fam|
	all_prot_data = ref_data[fam]
	prot_fam_obj = ProteinFamily.new(fam, all_prot_data)
	prot_fam_obj.split_into_classes(ref_data)
end

# delete unusual proteins from ref_data
Del_unusual_prots.each do |prot| 
	ref_data.delete(prot)
end

ref_data.each do |prot, all_prot_data|
puts prot
	prot_fam_obj = ProteinFamily.new(prot,all_prot_data)

	# clean_up alignments: remove common gaps && ensure same lenght
	prot_fam_obj.ensure_same_length
	prot_fam_obj.remove_common_gaps

	# save alignments to file
	f_out = File.join(new_data_path, "#{prot_fam_obj.prot_filesave_name}.fasta")
	prot_fam_obj.save_alignment(f_out)

	# test this by calling mafft
	f_test = File.join(new_data_path,"test.fasta")
	prot_fam_obj.ensure_mafft_is_fine(f_out, f_test)

	# calculate protein profile
	f_prfl = f_out.sub("fasta", "prfl")
	Helper.calc_protein_profile(f_out, f_prfl)

	# save genes and alignments (without common gaps!) back to ref_data
	ref_data[prot]["genes"] = prot_fam_obj.ref_genes
	ref_data[prot]["alignment"] = prot_fam_obj.ref_alignment.map{|k,v| ">#{k}\n#{v}"}.join("\n") + "\n"
end

# save json

puts "Saving updated reference data to file"
Helper.save_ref_data(new_ref_data_path, ref_data)



# create set of reference data without Ca_b (used as example for webserver)
ref_data = Helper.load_ref_data(new_ref_data_path)
ref_data.each do |prot, all_prot_data|
	puts prot

	prot_fam_obj = ProteinFamily.new(prot,all_prot_data)

	# remove all Ca_b genes and alignment sequences
	prot_fam_obj.ref_alignment.delete_if { |k,v| k =~ /Ca_b/ }
	prot_fam_obj.ref_genes.delete_if { |k,v| k =~ /Ca_b/ }

	# calc again: remove_common_gaps, save files, calc profile
	prot_fam_obj.remove_common_gaps
	f_out = File.join(new_data_path_no_cab,"#{prot_fam_obj.prot_filesave_name}.fasta")
	prot_fam_obj.save_alignment(f_out)
	
	f_prfl = f_out.sub("fasta", "prfl")
	Helper.calc_protein_profile(f_out, f_prfl)

	# save genes and alignment to ref_data
	ref_data[prot]["genes"] = prot_fam_obj.ref_genes
	ref_data[prot]["alignment"] = prot_fam_obj.ref_alignment.map{|k,v| ">#{k}\n#{v}"}.join("\n") + "\n"
end

puts "Save updated reference data (without Ca_b) to file"
Helper.save_ref_data(new_ref_data_path_no_cab, ref_data)

# move everything to place 

Helper.del_file_or_dir( Dir.glob(File.join(Base_path, "*.fasta")) )
Helper.del_file_or_dir( Dir.glob(File.join(Base_path, "*.prfl")) )
Helper.del_file_or_dir( File.join(Base_path, Without_ca_b) )
FileUtils.mv Dir.glob(new_data_path + "/*"), Base_path, :force => true
Helper.del_file_or_dir(new_data_path)
