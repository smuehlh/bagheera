module Helper
	extend self


	# general error handling methods
	def raise_runtime_error(msg="Fatal: An error occured")
		raise msg
	end

	def worked_or_die(is_success, msg)
		if ! is_success then
			raise_runtime_error "Fatal: #{msg}"
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
			raise_runtime_error "Fatal: File #{path} does not exist."
		end
	end

	def does_file_exist(path)
		if ! FileTest.file?(path) then
			return false
		end
		true
	end

	def dir_exist_or_die(path)
		if ! FileTest.directory?(path) then
			raise_runtime_error "Fatal: Directory #{path} does not exist."
		end
	end

	def mkdir_or_die(path)
		if ! FileTest.directory?(path) then
			begin
				# Dir.mkdir(path)
				FileUtils.mkdir_p(path)
			rescue
				raise_runtime_error "Fatal: Cannot create directory #{path}"
			end
		end
	end

	def move_or_copy_file(f_scr,f_dest,operation)
		file_exist_or_die(f_scr)
		case operation
		when "copy"
			FileUtils.cp(f_scr, f_dest)
		when "move"
			FileUtils.mv(f_scr, f_dest)
		end
	rescue
		raise_runtime_error "Error during setup. Please contact us."
	end

	def chmod(f_dest, mode)
		File.chmod(mode, f_dest)
	rescue
		raise_runtime_error "Error during setup. Please contact us."
	end

	# expecting a [params] file 
	def filesize_below_limit(file, max_size)
		if file.size > max_size then
			msg = []
			msg << "File must be less than #{(max_size/1024)/1024} MB"
			msg << "Please contact us to upload larger files."
			throw :error, msg
		end
		true
	end

	# load reference data from file alignment_gene_structure.json
	# @return [Hash] reference data
	# @return [Array] Errors occured during file load
	def load_ref_data
		path = File.join( ProteinFamily.class_variable_get(:@@ref_data_path), REF_DATA)
		file_exist_or_die(path)
		return JSON.load(File.read(path)) 
	end

	def get_tmp_file(extension="fasta")
		rand(1000000).to_s + "." + extension
	end

	def make_new_tmp_dir(base_path)
		id = rand(1000000000).to_s
		path = File.join(base_path, id)
		mkdir_or_die(path)
		return id
	end


	module Sequence
		extend self
		require 'yaml'

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


		def str2fasta(header, seq, no_split=false)
			fasta = header.include?(">") ? header << "\n" : ">" << header << "\n"
			if no_split then
				fasta += seq
			else
				fasta += seq.scan(/.{1,80}/).join("\n")
			end
			return fasta
		end

		# calculate position of residue in aligned sequence based on the position in unaligned sequence
		# spos: posion [int] in unaligned sequence
		# aseq: aligned sequence [string]
		def sequence_pos2alignment_pos(spos, aseq)
			pats = []
			aseq.gsub("-", "")[0..spos].split("").each {|chr| pats << ("-*" + chr)}
			pat = Regexp.new(pats.join)
			pat.match(aseq)[0].length - 1
		end

		# calculate position of residue in unaligned sequence based on the position in aligned sequence
		# apos: position [int] in aligned sequence
		# aseq: aligned sequence [string]
		def alignment_pos2sequence_pos(apos, aseq)
			aseq[0..apos].gsub("-", "").length - 1
		end

		def extract_translation(gene)
			g = Gene.new(YAML.load(gene))
			g.translation
		end

		def extract_cdna(gene)
			g = Gene.new(YAML.load(gene))
			g.cdna.upcase
		end

		def split_cdna_into_codons(gene, pos)
			dna_seq = extract_cdna(gene) # dna sequence
			codons = dna_seq.scan(/.{1,3}/)
			return codons[pos]
		end
	end

	class Gene
		require 'yaml'

		def initialize(blat_data)
			@blat_data = blat_data
			exon_number = 0
			@blat_data.each_with_index do |contig, contig_index|
			contig["number"] = contig_index.to_i + 1
			contig["matchings"].each do |matching|
			if matching["type"] == "exon" then
				exon_number += 1
				matching["number"] = exon_number.to_i
			elsif matching["type"] == "intron" ||  matching["type"] == "intron?" then
				matching["number"] = exon_number.to_i
			end
				matching["contig"] = contig_index.to_i + 1
			end
		end

		def contigs
			@blat_data
		end
		def matchings
			contigs.collect do |c| c["matchings"] end.flatten
		end
		def exons
			matchings.select do |m| m["type"] == "exon" || m["type"] == "exon_alternative" end
		end

		def cdna
			exons.select{|e| e["type"] != "exon_alternative" || e["alternative_type"] == "mutual_exclusive" || e["alternative_type"] == "duplicated_exon"}.collect do |e| e["seq"] end.join('').upcase
		end
		
		def translation
			exons.select{|ex| ex["type"] == "exon"}.collect do |m| m["translation"] end.join('')
		end

	end
    
end

end

class Array
    def sum
        self.inject{|sum,x| sum + x }
    end
	def find_each_index find
		found, index, q = -1, -1, []
		while found
			found = self[index+1..-1].index(find)
			if found
				index = index + found + 1
				q << index
			end
		end
		q
	end
end

class String
	def naturalized
		scan(/[^\d\.]+|[\d\.]+/).collect { |f| f.match(/\d+(\.\d+)?/) ? f.to_f : f }
	end
end
