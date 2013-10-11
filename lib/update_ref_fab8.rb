#!/usr/local/bin/ruby

### Update Bagheeras reference data on fab8
### Called monthly in crontab

# 00 04 1 * * ruby update_ref_neu.rb > /tmp/cug/cron-ruby.log 2> /tmp/cug/err-ruby.log
require 'timeout'
require 'fileutils'
require 'json'
require 'open3'

Base_path = "/tmp/cug/new"
Without_ca_b_path = "/tmp/cug/new/without_ca_b"
Ref_data_path = "/tmp/cug/new/alignment_gene_structure.json"

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

class ReferenceData
	attr_reader :path, :data

	def initialize(file)
		@path = file
		@data = load_ref_data
	end

	# separate by class & delete unseparated data & store back in ref_data
	def split_prot_into_classes(prot)
		prot_fam_obj = ProteinFamily.new(prot,@data[prot])	

		# first: delete unseparated data
		del_prot(prot)

		# second: separate this obj into classes
		keys = prot_fam_obj.ref_alignment.keys & prot_fam_obj.ref_genes.keys
		keys.each do |name|

			prot_abbr = Split_prot_abbrs[prot]

			klass = ""
			new_prot_key = ""

			matches = name.match(/#{prot_abbr}([0-9]+|Mhc)/)

			# Myosin protein has two abbreviations: Mhc for class 2 and Myo for all other classes
			if ! matches && name.include?("Mhc") then
				matches = ["", "Mhc"]
			end

			if matches && matches[1] then
				klass = matches[1]
				new_prot_key = prot + " Class " + klass
			else
				new_prot_key = prot
			end

			if ! @data.has_key?(new_prot_key) then
				add_prot(new_prot_key)
			end

			# third: save genes and alignments
			if prot_fam_obj.ref_genes.has_key?(name)
				@data[new_prot_key]["genes"][name] = prot_fam_obj.ref_genes[name]
			end

			if prot_fam_obj.ref_alignment.has_key?(name)
				@data[new_prot_key]["alignment"] += Helper::Sequence.str2fasta(name,prot_fam_obj.ref_alignment[name])
				@data[new_prot_key]["alignment"] += "\n"
			end

		end
	end

	def add_prot(key)
		@data[key] = {}
		@data[key]["genes"] = {}
		@data[key]["alignment"] = ""
	end

	def del_prot(key)
		@data.delete(key)
	end

	def update_alignment(key, fasta_hash)
		@data[key]["alignment"] = fasta_hash.map{|k,v| ">#{k}\n#{v}"}.join("\n") + "\n"
	end

	def update_genes(key, genes_hash)
		@data[key]["genes"] = genes_hash
	end

	def create_protfam_obj(key)
		obj = ProteinFamily.new(key,@data[key])	
	end

	def load_ref_data
		Helper.file_exist_or_die(@path)
		return JSON.load(File.read(@path)) 
	end

	def save_ref_data(this_path=@path)
		fh = File.new(this_path, "w")
		fh.puts JSON.dump(data)
		fh.close
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

	def ensure_mafft_is_fine(file)
		# test if mafft recognizes alignment by adding first sequence again to alignment
		tmp_file = File.join(Base_path, "test.fasta")
		File.open(tmp_file, 'w') { |f| f.write( Helper::Sequence.str2fasta(">test",@ref_alignment.values[0]) ) }
		is_success = Helper.run_mafft(file, tmp_file)
		Helper.del_file_or_dir(tmp_file)

		# if this is not successful, mafft does not recognize file as aligned
		if ! is_success then
			puts "#{@prot}:: need to calc alignment again"
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
# prepare directory structure

Helper.del_file_or_dir(Base_path)
Helper.mkdir_or_die(Base_path)
Helper.mkdir_or_die(Without_ca_b_path)

# fetch reference data from cymobase
max_secs = 60*10 # wait max. 10 minutes for cymoapi
puts "Start wget: #{Time.now}"
begin
	status = Timeout::timeout(max_secs) { callCymoAPI }
rescue Timeout::Error => exc
	Helper.abort(exc)
end
Helper.worked_or_die(status && Helper.does_file_exist(Ref_data_path), "wget failed")

ref_data_obj = ReferenceData.new(Ref_data_path)

# split protein families into protein classes
Split_prot_families.each do |fam|
	ref_data_obj.split_prot_into_classes(fam)
end

# delete unusual proteins 
Del_unusual_prots.each do |prot| 
	ref_data_obj.del_prot(prot)
end
ref_data_obj.save_ref_data

# 'reload' ref_data to have new set of keys in ref_data_obj
ref_data_obj = ReferenceData.new(Ref_data_path)
ref_data_wo_cab_obj = ReferenceData.new(Ref_data_path)

prot_list = ref_data_obj.data.keys
prot_list.sort.each do |prot|
	puts prot
	prot_obj = ref_data_obj.create_protfam_obj(prot)

	# 1) adapt alignments: ensure same lenght of all seqs and remove common gaps
	prot_obj.ensure_same_length
	prot_obj.remove_common_gaps

	# 2) save alignments to file and to reference data

	f_out = File.join(Base_path, "#{prot_obj.prot_filesave_name}.fasta")	
	Helper::Sequence.save_alignment(f_out, prot_obj.ref_alignment)
	ref_data_obj.update_alignment(prot, prot_obj.ref_alignment)

	# 3) test if mafft can handle alignment files
	prot_obj.ensure_mafft_is_fine(f_out)

	# 4) precalculate protein profiles
	f_prfl = f_out.sub("fasta", "prfl")
	Helper.calc_protein_profile(f_out, f_prfl)

	## do exactly the same again for dataset without example species 'Ca_b'
	prot_obj = ""
	prot_obj = ref_data_wo_cab_obj.create_protfam_obj(prot)

	# 5) remove all 'Ca_b' genes from genes and alignments
	prot_obj.ref_alignment.delete_if { |k,v| k =~ /Ca_b/ }
	prot_obj.ref_genes.delete_if { |k,v| k =~ /Ca_b/ }

	# 6) save alignments to file and to reference data
	f_out = File.join(Without_ca_b_path, "#{prot_obj.prot_filesave_name}.fasta")	
	Helper::Sequence.save_alignment(f_out, prot_obj.ref_alignment)
	ref_data_wo_cab_obj.update_alignment(prot, prot_obj.ref_alignment)
	ref_data_wo_cab_obj.update_genes(prot, prot_obj.ref_genes)

	# 3) precalculate protein profiles
	f_prfl = f_out.sub("fasta", "prfl")
	Helper.calc_protein_profile(f_out, f_prfl)

end

puts "Save reference data to file"
ref_data_obj.save_ref_data
ref_data_wo_cab_obj.save_ref_data( File.join( Without_ca_b_path, File.basename(Ref_data_path) ) )

puts "Move everything to place"
final_base_path = Base_path.sub("/new", "")
final_wo_ca_b_path = Without_ca_b_path.sub("new/", "")

Helper.del_file_or_dir( Dir.glob(File.join(final_base_path, "*.*")) ) # delete old fasta, prfl, json and log files
Helper.del_file_or_dir( final_wo_ca_b_path )
FileUtils.mv Dir.glob(Base_path + "/*"), final_base_path, :force => true
Helper.del_file_or_dir(Base_path)

