require "benchmark"
require 'peach'

class PredictionsController < ApplicationController

	MAX_SIZE = 52428800
	BASE_PATH = Dir::tmpdir + "/cug/" # resulting in final filename /tmp/cug/uploaded_genomefile_#{id}
	REF_DATA = BASE_PATH + "alignment_gene_structure.json"
	SORT = "/usr/bin/sort"
	FORMATDB = "/usr/bin/formatdb"
	BLASTALL = "/usr/bin/blastall"
	FASTACMD = "/usr/bin/fastacmd"
	AUGUSTUS = "/usr/local/bin/augustus/src/augustus"
	PAIR_ALIGN =  Rails.root + "lib/pair_align" #"/usr/local/bin/pair_align"
	HYDORPHOBIC_AAS = ["V", "I", "L"]
	POLAR_AAS = ["S", "T"]

	def search
		delete_old_data # delete old uploaded genome files
		prepare_new_session # a fresh session
	end

	def upload_file
		`fromdos #{params[:uploaded_file].path}` # remove \r, if uploaded from windows
		# is the file too big? maybe jquery-check did not work, so better check again
		@content_error = check_filesize(params[:uploaded_file].size)
		# check content, only if file is not too big
		@content_error |= check_fasta(params[:uploaded_file]) if @content_error.blank?
		if @content_error.blank? then
			# rename and save file
			file_id = rand(1000000000).to_s
			file_name = File.basename(params[:uploaded_file].original_filename)
			file_path = BASE_PATH + "uploaded_genomefile_" + file_id + ".fasta"
			File.rename(params[:uploaded_file].path, file_path)
			File.open(file_path, 'wb') {|file| file.write(params[:uploaded_file].read)}
			session[:file] = { id: file_id, name: file_name, path: file_path }
		end
		render :upload_file_ajax
	end

	def load_example
		# do not check content, but check file size, just to be sure
		example_file = "Candida_albicans_WO_1.fasta"
		example_path = Rails.root + "spec/fixtures/files/" + example_file
		@content_error = check_filesize(File.size(example_path))
		if @content_error.blank? then
			file_id = rand(1000000000).to_s
			file_path = BASE_PATH + "uploaded_genomefile_" + file_id + ".fasta"
			File.copy_stream(example_path, file_path)
			@file = example_file
			session[:file] = { id: file_id, name: example_file, path: file_path }
		end
		render :upload_file_ajax
	end

# uncomment if alignment options should get an own submit-button instead of the big "Predict button"
	# def set_alignment_options
	# 	session[:align] = { algo: params[:algo], config: params[:config] }
	# 	render :set_options
	# end

	def predict_genes
		# general workflow

		# 1) extract reference proteins
		# 2) gene prediction foreach reference protein:
		# 2.1) BLAST
		# 2.2) AUGUSTUS
		# 3) Compare with reference data
		# 3.1) Compare with reference alignment
		# 3.2) Compare with reference genes

		# implementation
# TODO: benchmark entfernen
time = Benchmark.realtime do

		# add alignment options to session
		session[:align] = { algo: params[:algo], config: params[:config] }
		ref_data, @errors = load_ref_data
		@predicted_prots = {} # containing final results
		if @errors.empty? then
			# sucessfully loaded reference data file

			# for BLAST (2.1): convert uploaded genome into blast-db
			genome_db = session[:file][:path].gsub(".fasta", "_db")
				# -i [INPUT: GENOME FILE] -p [PROTEIN] FALSE -n [DB NAME] -o [CREATE INDEX OVER SEQID] TRUE
			output = `#{FORMATDB} -i #{session[:file][:path]} -p F -n #{genome_db} -o T 2>&1`
			if ! $?.success? || output.include?("ERROR") then
				@errors << "Gene prediction failed: Cannot create BLAST database"
				next
			end 
			# all other tasks need to be done for each ref protein
			ref_data.peach(10) do |prot, all_prot_data|
				### 1) extract reference proteins
				ref_prot_species, ref_prot_stat, ref_prot_seq, ref_prot_geneseq, ref_prot_key = "", "", "", "", ""

				sc_genes = all_prot_data["genes"].keys.select {|key| key =~ /Sc_/ && all_prot_data["genes"][key]["completeness"] =~ /[complete|partial]/}
				if sc_genes.any? then
					# Sc_b and complete gene structure
					ref_prot_key = sc_genes.find {|key| all_prot_data["genes"][key]["completeness"] == "complete" &&  all_prot_data["genes"][key] =~ /Sc_b/}
					# ... or ... other Sc and complete gene structure
					ref_prot_key = sc_genes.find {|key| all_prot_data["genes"][key]["completeness"] == "complete"} if ref_prot_key.blank?
					# ... or ... Sc_b and partial
					ref_prot_key = sc_genes.find {|key| key =~ /Sc_b/} if ref_prot_key.blank?
				end
				# ... or ... use first gene with complete gene structure
				ref_prot_key = all_prot_data["genes"].keys.find {|key| all_prot_data["genes"][key]["completeness"] =~ /complete/} if ref_prot_key.blank?
				# ... or ... or just any gene at all
				ref_prot_key = all_prot_data["genes"].keys.first if ref_prot_key.blank?
				# fill other ref_data-types
				ref_prot_species = all_prot_data["genes"][ref_prot_key]["species"]
				ref_prot_stat = all_prot_data["genes"][ref_prot_key]["completeness"]
				ref_prot_seq = all_prot_data["genes"][ref_prot_key]["gene"].match(/prot_seq: (.*?)\n/)[1]
				ref_prot_geneseq = extract_gene_seq(all_prot_data["genes"][ref_prot_key]["gene"])
				prot_basename = prot.gsub(" ", "-").downcase

				### 2) gene prediction foreach reference protein:

				### 2.1) BLAST
				# 2.1.1) convert uploaded genome into blast-db
				# 2.1.2) blast-query
				file_in = BASE_PATH + "blast_in_" + prot_basename + "_" + session[:file][:id] + ".fasta"
				# create file_in (protein file) for blast query
				fasta = str2fasta(prot_basename, ref_prot_seq)
				File.open(file_in, 'w'){|file| file.write(fasta)}
				# 2.1)
						# -p [PROGRAM] protein query against nt-db -d [DATABASE] -m8 [OUTPUT FORMAT] | sort by e-value
				output = `#{BLASTALL} -p tblastn -d #{genome_db} -i #{file_in} -m8 2>&1 | #{SORT} -g -k 11 | head -1`
				if ! $?.success? || output.include?("ERROR") then
					@errors << "#{prot}: Gene prediction (BLAST) failed"
					next
				end
				# parse blast output and add start/stop position to selected_ref_prots
				output_fields = output.chomp.split("\t") # best hit
				seq_id = output_fields[1] # sequence contig
				start = output_fields[8].to_i # sequence start
				stop = output_fields[9].to_i # sequence stop
				strand = 1 # set strand to plus for fastacmd
				if stop < start then
					# change start/ stop if prediction on minus strand 
					strand = 2 # set strand to minus for fastacmd
					tmp = start
					start = stop
					stop = tmp
				end
				### 2.2) AUGUSTUS
				# 2.2.1) get matching search sequence 
				file_in.gsub!("blast_in", "aug_in")
						# -d [DB] -p [PROTEIN] false -s [SEARCH STRING] -L [START,STOP]
						# add 1000 nucleotides to start/stop 
				output = `#{FASTACMD} -d #{genome_db} -p F -s \'#{seq_id}\' -L #{start-1000},#{stop+1000} -S #{strand} > #{file_in}` 
				if ! $?.success? || output.include?("ERROR") then
					@errors << "#{prot}: Gene prediction (BLAST-FASTACMD) failed"
					next
				end
				# 2.2.2) augustus
						# --species [REFERENCE SPEC] QUERY
	 			output = `#{AUGUSTUS} --AUGUSTUS_CONFIG_PATH=/usr/local/bin/augustus/config/ --species=saccharomyces_cerevisiae_S288C #{file_in}`
				if ! $?.success? || output.include?("ERROR") then
					@errors << "#{prot}: Gene prediction (AUGUSTUS) failed"
					next
				end
				# 2.2.3) parse augustus output
				# return longest predicted protein and according coding sequence
				pred_seq, pred_dnaseq = parse_augustus(output, fasta2str(File.read(file_in)))
				pred_seq.upcase!

				### 3.1) Compare with reference alignment
				@predicted_prots[prot] = {has_ctg: true, pred_prot: pred_seq, ctg_pos: [], ctg_transl: [], aa_comp: [], ctg_ref: []}
				# tracking of ctg_pos, ctg_transl, aa_comp, ctg_ref only for significant positions!

				# has predicted protein CTG codons?
				codons = pred_dnaseq.scan /.{1,3}/
				if !codons.include?("CTG")
					# nope, no CTG => nothing to predict
					@predicted_prots[prot][:has_ctg] = false
					next
				end
				ctg_pos = codons.each.with_index.map{ |ele, ind| (ele == "CTG") ? ind : nil }.compact

				### 3.1.1) map predicted gene onto reference protein
				file_in.gsub!("aug_in", "align_in")
				file_out = file_in.gsub("align_in", "align_out")

				fasta = str2fasta("ref", ref_prot_seq) << "\n" << str2fasta("pred", pred_seq)
				File.open(file_in, 'w'){|file| file.write(fasta)}
					# -c tttt: initialize first row & column with zeros, search last row & column for maximum
					# 	=> end gap free alignment
					# --matrx matrix file
					# implicit options: gap open penalty: -11, gap extension penalty: -1
					# 					protein sequences
				output = `#{PAIR_ALIGN} --seq #{file_in} --matrix #{Rails.root+"lib/blosum62"} --config #{session[:align][:config]} --method #{session[:align][:algo]} --outfile #{file_out}`
				if ! $?.success? || ! output.include?("Alignment score") then
					@errors << "#{prot}: Codon usage prediction (Seqan::pair_align) failed"
					next
				end
				ref_seq_aligned, pred_seq_aligned = parse_seqan(File.read(file_out))

				### 3.1.2a) map CTG positions in predicted protein onto reference alignment columns
				### 3.1.2b) compare with chemical properties of reference alignment columns 
				cymo_algnmnt = Hash[*ref_data[prot]["alignment"].split("\n")] # cymo-alignment
				ctg_pos.each do |pos|

					algnmnt_col_aa, algnmnt_col_pos = parse_alignment_by_ctgpos(cymo_algnmnt, pos, pred_seq_aligned, ref_seq_aligned, ref_prot_key)
					aa_freq = word_frequency(algnmnt_col_aa)
					is_significant, prob_transl = predict_translation(aa_freq)
					if is_significant then
						@predicted_prots[prot][:ctg_pos] << pos
						@predicted_prots[prot][:ctg_transl] << prob_transl
						@predicted_prots[prot][:aa_comp] << aa_freq
					end

					### 3.2) Compare with reference genes: count usage of CTG for each mapped L and S
					pct_ctg_S, pct_ctg_L = ref_ctg_usage(cymo_algnmnt, algnmnt_col_pos, ref_data[prot]["genes"])
					@predicted_prots[prot][:ctg_ref] << {ser: pct_ctg_S, leu: pct_ctg_L}
				end
			end # ref_data.peach do |prot, des|
		end # if @errors.empty?
end
puts "Time elapsed #{time*1000} milliseconds"
		render :predict_genes
	end

	def check_filesize(size)
		errors = []
		if size > MAX_SIZE then
			errors << "File must be less than 50 MB."
			errors << "Please contact us to upload larger files."
		end
		return errors
	end

	def check_fasta(file)
		errors = []
		is_fasta = true
		contains_ctg = false
		is_first_line = true
		last_nucleotides = "" # need to store 2 last nts of each line to check for ctg!
		# read file line by line, since it might be large
		IO.foreach(file.path) do |line|
			line.chomp!
			if is_first_line then
				# expect very first line to be fasta header
				if ! line.starts_with?('>') then
					is_fasta = false
					# test failed, stop reading file!
					break
				end
				is_first_line = false
			else
				# it's not the very first line, might be header, comment or sequence
				if line.starts_with?('>') then
					# a new header: clear last_nucleotides-buffer
					last_nucleotides = ""
					next
				elsif line.starts_with?(';')
					# comment line, simply skip
					next
				elsif line =~ /[^ACGTURYMKSWBDHVN]/i
					# a sequence line, but containing charaters not part of iupac nucleotide definition
					# test failed, stop reading file!
					is_fasta = false
					break
				else
					# a sequence line, valid content
					# ctg might be splitted into 2 lines
					contains_ctg = true if (last_nucleotides + line).upcase.include?("CTG")
					last_nucleotides = line[-2,2]
				end
			end
		end

		# file read, check what happend
		if ! is_fasta then
			errors << "Invalid input."
			errors << "Expected fasta formatted genome file."
		end
		if  is_fasta && (! contains_ctg) then
			errors << "Genome sequence contains no \'CTG\'."
			errors << "Cannot predict codon usage."
		end
		return errors
	end

	# expects header and fasta-sequence as string and converts both in fasta-format
	def str2fasta(header, str)
		fasta = header.include?(">") ? header << "\n" : ">" << header << "\n"
		# breaking string at every 80th character, joining by newline
		fasta += (str.scan /.{1,80}/).join("\n")
		return fasta
	end

	# expects fasta-formatted string (header + seq, containing \n) and returns the sequence
	def fasta2str(fasta)
		return fasta.split("\n")[1..-1].join("")
	end

	# extracting gene coding dna from gene entry in reference data
	def extract_gene_seq(gene)
		a = gene.scan(/\sseq:\s?(\w+)\n/).flatten
		# quick & dirty version: every second entry belongs to an exon
		a.values_at(*a.each_index.select(&:even?)).join("").upcase
	end

	def delete_old_data(days = 1)
		system("find #{BASE_PATH}uploaded_genomefile_* -mtime +#{days} -type f -delete 2&> /dev/null")
		system("find #{BASE_PATH}formatdb.log -mtime +#{days} -type f -delete 2&> /dev/null")
		system("find #{BASE_PATH}*.fasta -mtime +#{days} -type f -delete 2&> /dev/null")
	end

	def load_ref_data
		if ! FileTest.file?(REF_DATA) then
			errors = ["Sorry, cannot load reference data. Please contact us!"]
			return false, errors
		end
		return JSON.load(File.read(REF_DATA)), []
	end

	def prepare_new_session
		reset_session
	    session[:file] = {} 		# uploaded genome described by keys: :name, :id, :path
	    session[:align] = {}		# alignment options for seqan: pair_align
	    return true
	end

	def parse_augustus(output, dna_seq)
		# input: augustus output
		#        dna sequence used by augustus

		# output: predicted protein sequence
		#         coding sequence of predicted protein

		if output =~ /start gene g2/ then
			# more than one hit predicted, find longest one
			hits = output.scan(/\[(.*?)\]\n# end gene g(\d+)/m) # [[prot, hit], [prot, hit], ...]
			pred_prot, best_hit = hits.max_by{|ele| ele[0].length}
			# combine with right gene!!! gene_id "g1"
			cds_pos = output.scan(/AUGUSTUS\tCDS\t(\d+)\t(\d+)\t.*?gene_id "g#{best_hit}"/)
		else
			# only one hit predicted
			pred_prot = output.match(/\[(.*)\]/m)[1]
			cds_pos = output.scan(/AUGUSTUS\tCDS\t(\d+)\t(\d+)\t/) # [[exon1_start, exon1_stop], [exon2_start, exon2_stop], ...]
		end

		# prepare output
		pred_prot.gsub!("\n# ", "")
		pred_dna = cds_pos.inject("") {|res, exon| res += dna_seq[(exon[0].to_i - 1) .. (exon[1].to_i - 1)].upcase }
		return pred_prot, pred_dna
	end

	def parse_seqan(file)
		parts = file.split("\n>")
		return fasta2str(parts[0]), fasta2str(parts[1])
	end

	def alignment_pos2sequence_pos(apos, aseq)
		aseq[0..apos].gsub("-", "").length - 1
	end

	def sequence_pos2alignment_pos(spos, aseq)
		pats = []
		aseq.gsub("-", "")[0..spos].split("").each {|chr| pats << ("-*" + chr)}
		pat = Regexp.new(pats.join)
		pat.match(aseq)[0].length - 1
	end

	# get alignment column corresponding to ctg position in predicted protein
	def parse_alignment_by_ctgpos(algnmnt, pred_spos, pred_aseq, ref_aseq, ref_key)
		# 1) ctg pos -> pos in predicted sequence aligned with reference seq
		pred_apos = sequence_pos2alignment_pos(pred_spos, pred_aseq) # = pos_ref_seq_aligned
		# 2) pos aligned -> cymo alignment
		ref_spos = alignment_pos2sequence_pos(pred_apos, ref_aseq)
		ref_cymopos = sequence_pos2alignment_pos(ref_spos, algnmnt[">"<<ref_key])
		# 3) return column
		col = algnmnt.collect {|cymo_prot, cymo_seq| cymo_seq[ref_cymopos]}
		return col, ref_cymopos
	end

	# counting the frequency (absolute number) of words (amino acids...) in an array
	# returning hash of aa1: count; aa2: count
	def word_frequency(arr)
		res = Hash.new(0)
		arr.each { |a| res[a] += 1 }
		res.delete("-") # delete count for gaps in alignment!
		return res
	end

	# def predict_translation_old(aas)
	# 	# statistics
	# 	# use only amino acids, no gaps: exclude "-" from statistics
	# 	num_aas = aas.count{|ele| ele =~ /[A-Z]/} # total number of amino acids
	# 	pol_aas = aas.count{|ele| POLAR_AAS.include?(ele)} # polar amino acids
	# 	hyd_aas = aas.count{|ele| HYDORPHOBIC_AAS.include?(ele)} # hydrophobic amino acids
	# 	pct_pol = pol_aas/num_aas.to_f
	# 	pct_hyd = hyd_aas/num_aas.to_f
	# 	# requirements for a discriminative position:
	# 		# 1) occurence in more than half of sequences
	# 		# 2) the other usage should occure in less than quater or sequences
	# 	is_discrim = ([pct_hyd, pct_pol].max > 0.5 && [pct_hyd, pct_pol].min < 0.25) ? true : false 
	# 	# default translation = ""
	# 		# discriminative AND more hydrophobic aas: "L"
	# 		# discriminative AND more polar aas: "S"
	# 	transl = ""
	# 	if is_discrim then
	# 		transl = (pol_aas > hyd_aas) ? "S" : "L"
	# 	end
	# 	return is_discrim, transl, pct_pol.round(2), pct_hyd.round(2)
	# end

	def predict_translation(aa_freq)
		num_aas = aa_freq.inject(0){|res, obj| res += obj[1]} # total number of amino acids
		pol_aas = POLAR_AAS.collect{|aa| aa_freq[aa]}.sum # polar amino acids
		hyd_aas = HYDORPHOBIC_AAS.collect{|aa| aa_freq[aa]}.sum # hydrophobic amino acids
		pct_pol = pol_aas/num_aas.to_f
		pct_hyd = hyd_aas/num_aas.to_f
		# requirements for a discriminative position:
			# 1) occurence in more than half of sequences
			# 2) the other usage should occure in less than quater or sequences
		is_discrim = ([pct_hyd, pct_pol].max > 0.5 && [pct_hyd, pct_pol].min < 0.25) ? true : false 
		# default translation = ""
			# discriminative AND more hydrophobic aas: "L"
			# discriminative AND more polar aas: "S"
		transl = ""
		if is_discrim then
			transl = (pol_aas > hyd_aas) ? "S" : "L"
		end
		return is_discrim, transl
	end

	# check for a position in the cymo alignment, how many "S" and how many "L" are encoded by an CTG
	# number of other codons = 100 - pct_ctg
	# input: 
	#  		cymo alignment, column in cymo alignment, ref_data[prot]["genes"]	 
	def ref_ctg_usage(algnmnt, apos, genes)
		pct_ctg_S, pct_ctg_L = "", "" 
		# Serine
		ref_codons = []
		algnmnt.collect{|k,v| k if v[apos]=="S"}.compact.each do |ref_key|
			dna_seq = extract_gene_seq(genes[ref_key.gsub(">", "")]["gene"]) # access via "gene-name"
			codons_ref = dna_seq.scan /.{1,3}/
			spos = alignment_pos2sequence_pos(apos, algnmnt[ref_key]) # access via ">gene-name"
			ref_codons << codons_ref[spos]
		end
		pct_ctg_S = ref_codons.count("CTG")/ref_codons.length.to_f 

		# same for Leucine
		ref_codons = []
		algnmnt.collect{|k,v| k if v[apos]=="L"}.compact.each do |ref_key|
			dna_seq = extract_gene_seq(genes[ref_key.gsub(">", "")]["gene"]) # access via "gene-name"
			codons_ref = dna_seq.scan /.{1,3}/
			spos = alignment_pos2sequence_pos(apos, algnmnt[ref_key]) # access via ">gene-name"
			ref_codons << codons_ref[spos]
		end
		pct_ctg_L = ref_codons.count("CTG")/ref_codons.length.to_f
		return pct_ctg_S.round(2), pct_ctg_L.round(2)
	end
end
