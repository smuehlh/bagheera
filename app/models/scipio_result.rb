class ScipioResult

#  include DRbUndumped
  require 'yaml'
  require 'yaml2log'
  require 'pretty_aligner'
  attr_reader :blat_data, :name, :comment_lines


  TYPES_OF_DISCREPANCY = {
    1 => "mismatch",
    2 => "undetermined query",
    3 => "undetermined target",
    4 => "additional codon in target",
    5 => "unmatched query",
    6 => "frameshift (+1) target only",
    7 => "frameshift (+1) target/query",
    8 => "frameshift (+2) target only",
    9 => "frameshift (+2) target/query",
    10 => "frameshift (-2) target only",
    11 => "frameshift (-2) target/query",
    12 => "frameshift (-1) target only",
    13 => "frameshift (-1) target/query",
    14 => "stopcodon target/query",
    15 => "stopcodon, target only",
    16 => "stopcodon, undetermined query",
    17 => "additional stopcodon"
  }


  def ScipioResult.new_from_file(fn)
    
    return ScipioResult.new_from_string(File.read(fn))

  end

  
  def ScipioResult.new_from_string(str)

    comment_lines = []
    str.each_line do |line|
      comment_lines << line.rstrip if line.match(/\A#+/)
    end

    return ScipioResult.new(YAML.load(str), "Unnamed Sequence", comment_lines)

  end
  

  def initialize(blat_data, name = "Unnamed Sequence", comment_lines = [])
    @comment_lines = comment_lines
    @blat_data = blat_data
    @name = name

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

  end


  def <=>(other)
    self.name <=> other.name
  end


  def contigs
    @blat_data
  end


  def is_tandem_gene?
    contigs.any?{|contig| contig["tandem_gene_number"]}
  end

  def has_tandem_genes?
    contigs.any?{|contig| contig["tandem_genes"] && !contig["tandem_genes"].empty?}
  end

  def has_alternative_exons?
    matchings.any? {|match| match["type"] == "exon_alternative"}
  end

  def has_duplicated_exons?
    matchings.any? {|m| m["type"] == "exon_alternative" && m["alternative_type"] == "duplicated_exon"}
  end

  def alternative_exons
    matchings.select {|m| m["type"] == "exon_alternative"}
  end

  def has_mutually_exclusive_exons?
    matchings.any? {|match| match["type"] == "exon_alternative" && match["alternative_type"] == "mutual_exclusive"}
  end


  def mutually_exclusive_exons
    matchings.select {|m| m["type"] == "exon_alternative" && m["alternative_type"] == "mutual_exclusive"}
  end

  def duplicated_exons
    matchings.select {|m| m["type"] == "exon_alternative" && m["alternative_type"] == "duplicated_exon"}
  end

  def exons_original
    exons.select{|exon| exon["type"] == "exon"}
  end


  def exons_have_alternatives
    exons_to_return = []
    exons_original.each_with_index {|exon, i| exons_to_return << exon if exon_has_alternatives?(i+1)}
    return exons_to_return
  end


  def exons_have_mutually_exclusives
    exons_to_return = []
    exons_original.each_with_index {|exon, i| exons_to_return << exon if exon_has_mutually_exclusives?(i+1)}
    return exons_to_return
  end


  def tandem_gene_exons
    return self.exons.select{|exon| exon["tandem_gene"]}
  end


  def tandem_genes

    tandem_genes = []

    tandem_gene_number = 1
    loop do
      exons_of_tandem_gene = self.tandem_gene_exons.select{|tandem_gene_exon| tandem_gene_exon["tandem_gene"].to_i == tandem_gene_number}
      break if exons_of_tandem_gene.empty?
      tandem_genes << exons_of_tandem_gene
      tandem_gene_number += 1
    end

    return tandem_genes
  end


  def ScipioResult.tandem_gene_score(exons_of_tandem_gene, protein_translation_length)

    number_of_aa_found = 0

    exons_of_tandem_gene.each do |exon_of_tandem_gene|
      (0...exon_of_tandem_gene["exon_al"].length).each do |aa_index|
        number_of_aa_found += 1 if exon_of_tandem_gene["exon_al"][aa_index].chr != "-" &&
          exon_of_tandem_gene["cand_al"][aa_index].chr != "-"
      end
    end

    score = number_of_aa_found.to_f / protein_translation_length.to_f * 100.0

    return score

  end


  def tandem_gene_score(exons_of_tandem_gene)

    return ScipioResult.tandem_gene_score(exons_of_tandem_gene, self.translation.length)

  end

  
  def tandem_gene_score_by_number(tandem_gene_number)

    exons_of_tandem_gene = self.tandem_gene_exons.select{|exon| exon["tandem_gene"].to_i == tandem_gene_number.to_i}
    return ScipioResult.tandem_gene_score(exons_of_tandem_gene, self.translation.length)

  end

  


  def yaml
    yamlstring = ""
    yamlstring += @comment_lines.join("\n") + "\n" + "\n" unless @comment_lines.empty?
    yamlstring += @blat_data.to_yaml

    return yamlstring
  end

  
  def query_length_aa
    contigs[0]["prot_len"]
  end
    
  def query_length_nuc
    query_length_aa * 3
  end
    
  def exon_length
#    puts "DEBUG: exon_length: #{@name}: #{exons.select{|exon| exon["seq"].length != (exon["dna_end"]-exon["dna_start"]).abs}.collect{|exon| exon["number"]}.inspect}"
    exons.sum{|exon| exon["seq"] ? exon["seq"].length : 0}
  end
    
    
  ## return merged ranges of the nucleotides that were found
  ## e.g. [66..90, 165..204]
  def matched_nucl_ranges
    ranges = []
    ranges = exons.collect { |e| (e["nucl_start"].to_i..e["nucl_end"].to_i) if e["nucl_start"].to_i <= e["nucl_end"].to_i}
    ranges_merged = []
    ranges.each do |r|
      if ranges_merged[-1] && (m = ranges_merged[-1] + r) then
        ranges_merged[-1] = m
      else
        ranges_merged << r
      end
    end
    return ranges_merged
  end
    
  def mismatches
    exons.collect do |e| e["mismatchlist"] end.flatten.uniq.compact
  end
    
  def seqshifts
    exons.collect do |e|
      if (ss = e["seqshifts"]) then
        ss
        #ss.collect do |s| (s["nucl_start"]..s["nucl_end"]).to_a end
      else
        nil
      end
    end.flatten.uniq.compact
  end
    
  def matched_nucls
    matched_nucl_ranges.collect do |r| r.to_a end.flatten.uniq
  end
    
  def unmatched_nucls
    (1..query_length_nuc).to_a - matched_nucls
  end
    
  def nucl_match_ratio
    (matched_nucls.length.to_f - (mismatches.length.to_f + seqshifts.length.to_f)) / query_length_nuc.to_f
  end

  def gaps
    matchings.select do |match|match["type"] == "gap" end
  end

  def gap_ratio
    (matched_nucls.length.to_f / query_length_nuc.to_f)
  end
    
  def gap_length
    length = 0
    matchings.each do |match|
      if (match["type"] == "gap")
        #  l = match["dna_end"] - match["dna_start"]
        #   puts "l:#{l}"
        #length += (match["dna_end"] - match["dna_start"])
        length += match["seq"].length
      end
    end
    #    puts "gap-length:#{length}"
    return length
  end
    
  def gap_count
    return gaps.size
  end

  def contig_gap_count
    gap_contig_lengths.select{|g| g > 0}.count
  end
    
  def intron_putative_length
    length = 0
    matchings.each do |match|
      if (match["type"] == "intron?")
        length += (match["dna_end"]- match["dna_start"])
      end
    end
    #   puts length
    return length
  end
    
  def contiggap_length
    gap_contig_lengths.sum
  end
    
  def nucl_match_ratio_color
    percentage = (nucl_match_ratio * 100).round
    color =
      case percentage
    when 100
      "green"
    when 95..99
      "yellow"
    else
      "red"
    end
    return color
  end


  def complete?
    #   puts "complete?"
    #    puts log.scan(/^status\s.+incomplete/).empty?
    #    log.scan(/^status\s.+incomplete/).empty?
    return !yaml.match(/incomplete/)
  end


  def perfect_match?
    unmatched_nucls.empty? && mismatches.empty? && seqshifts.empty?
  end

  def perfect_match_but_mismatches_seqshifts?
    unmatched_nucls.empty?
  end
      
  def contig_names
    contigs.collect do |c| c["target"] end
  end


  def contig_strands
    contigs.collect do |c| c["strand"] end
  end

  def contig_lengths
    contigs.collect do |c| (c["dna_end"] || 0) - (c["dna_start"] || 0) end
  end

    
  def contig(name)
    contigs.select do |c| c["target"] == name end
  end
    
  def downstream_gaps
    gaps = Array.new
    c = contigs.collect do |c| c["downstream_gap"] end
    c.each do |downstream_gap|
      if (downstream_gap != nil)
        gaps << downstream_gap.length
      else
        gaps << 0
      end
    end
    return gaps
  end
    
  def upstream_gaps
    gaps = Array.new
    c = contigs.collect do |c| c["upstream_gap"] end
    c.each do |upstream_gap|
      if (upstream_gap != nil)
        gaps << upstream_gap.length
      else
        gaps << 0
      end
    end
    return gaps
  end
                
  def matchings
    contigs.collect do |c| c["matchings"] end.flatten
  end

  def introns
    matchings.select do |m| m["type"] == "intron" || m["type"] == "intron?" end
  end

   def uncertain_introns
    matchings.select do |m| m["type"] == "intron?" end
  end

  def exons
    matchings.select do |m| m["type"] == "exon" || m["type"] == "exon_alternative" end
  end
    
  def cdna
    exons.select{|e| e["type"] != "exon_alternative" || e["alternative_type"] == "mutual_exclusive" || e["alternative_type"] == "duplicated_exon"}.collect do |e| e["seq"] end.join('').upcase
  end

  def genomic_dna
    matchings.select{|e| e["type"] != "exon_alternative" || e["alternative_type"] == "mutual_exclusive" || e["alternative_type"] == "duplicated_exon"}.collect do |m| 
      m["seq"]
    end.join('')
  end
    
  def upstream_dna
    upstream_dna = contigs.collect do |c| c["upstream"] end.compact[0]
    return upstream_dna
  end
    
  def downstream_dna
    downstream_dna = contigs.collect do |c| c["downstream"] end.compact[0]
    return downstream_dna
  end
    
  def translation
    exons.select{|ex| ex["type"] == "exon"}.collect do |m| m["translation"] end.join('')
  end
    
  def scipio_yaml

    yamlstring = ""
    yamlstring += @comment_lines.join("\n") + "\n" + "\n" unless @comment_lines.empty?

    yamlstring += ZAML::dump({@name => @blat_data})
    #yamlstring += YAML::dump({@name => @blat_data})

    return yamlstring
  end

  #returns array with uniq putative exon numbers
  def count_alternative_exons
    matchings.select do |m| m["type"] == "exon_alternative" end.collect do |m| m["exon"] end.uniq
  end
    
  def log

    return @logfile if @logfile

    duplicated_exons = self.exons_original.select{|e| e["alternative_type"] == "duplicated_exon"}
    duplicated_exons.each do |dup_exon|
      dup_exon["type"] = "duplicated_exon"
    end

    @logfile = Yaml2Log.logfile(scipio_yaml)

    duplicated_exons.each do |dup_exon|
      dup_exon["type"] = "exon"
    end


    # Add alignment for alternative exons
    pretty_aligner = PrettyAligner.new
   
    blat_data_with_cdna = self.blat_data.select{|hit_data| hit_data["cdna_seq"]}

    if blat_data_with_cdna.size > 0 then

      alt_exons = blat_data_with_cdna.collect do |hit_data|
        hit_data["matchings"].select do |m|
          (m["type"] == "exon_alternative" && m["alternative_type"] != "mutual_exclusive" && m["alternative_type"] != "duplicated_exon")
        end
      end
      alt_exons.flatten!

      gaps = blat_data_with_cdna.collect do |hit_data|
        hit_data["matchings"].select do |m|
          m["type"] == "gap"
        end
      end

      gaps.flatten!
      gap_sequences = gaps.collect do |gap|
        gap["cdnaseq"]
      end


#      index_first_sequence = @logfile.index(/upstream\n|exon\n|intron\n|exon_alternative\n/)
#
#      puts "DEBUG: #{blat_data_with_cdna[0]["matchings"].sort_by{|m| m["dna_start"]}[0]["cdna_start"].inspect}"
#
#      cdna_upstream_end = blat_data_with_cdna[0]["matchings"].sort_by{|m| m["dna_start"]}[0]["cdna_start"].to_i
#      cdna_upstream_seq = blat_data_with_cdna[0]["cdna_seq"][0, cdna_upstream_end]
#
#      puts "DEBUG: cdna_upstream_seq: #{cdna_upstream_seq}"
#
#      @logfile.insert(index_first_sequence, "cdna_upstream_seq\n#{cdna_upstream_seq}\n\n")
#
##      puts "DEBUG: index_first_sequence: #{index_first_sequence}"
#

#      alt_exons.each {|e| puts e["seq"]}

      
      alt_exon_index = -1
      @logfile = @logfile.gsub(/exon_alternative.+?(?=intron|exon|exon_alternative|gap|upstream|downstream|\Z)/m) {|matched_text|
        "exon_alternative\n#{pretty_aligner.get_alignment_cdna(alt_exons[alt_exon_index+=1])}"
      }

      i = -1
      @logfile = @logfile.gsub(/gap.+?(?=intron|exon|exon_alternative|gap|upstream|downstream|\Z)/m) {|matched_text|
        "#{matched_text}\ngap_cdna\n#{gap_sequences[i+=1]}\n"
      }

    else

      dup_exon_index = -1
      @logfile = @logfile.gsub(/duplicated_exon(.|\n)+?(intron|gap|downstream|\Z)/) {|alt_exon|
        alt_exon_match = alt_exon.match(/duplicated_exon(.|\n)+?(intron|gap|downstream|\Z)/)
        "exon\n#{pretty_aligner.get_alignment(duplicated_exons[dup_exon_index+=1])}#{alt_exon_match[2]}"
      }

      alt_exons = self.alternative_exons.select{|alt_exon| alt_exon["alternative_type"] == "mutual_exclusive" || alt_exon["alternative_type"] == "duplicated_exon"}
      alt_exon_index = -1
      @logfile = @logfile.gsub(/exon_alternative(.|\n)+?(intron|gap|downstream|\Z)/) {|alt_exon|
        alt_exon_match = alt_exon.match(/exon_alternative(.|\n)+?(intron|gap|downstream|\Z)/)
        "exon_alternative\n#{pretty_aligner.get_alignment(alt_exons[alt_exon_index+=1])}#{alt_exon_match[2]}"
      }
      
    end

    return @logfile
  end
    
  def contiginfo(i)
    
    return @info if @info

    @info = Hash.new

    @info["identity"]=log.scan(/^identity.*?(\d+\.\d+%)/)[i]
    @info["score"]=log.scan(/^score.*?(\d+\.\d+)/)[i]
    #@info["check"]=log.scan(/^check\s*?(\w+)/)[i]
      
    return @info
  end
    
  #bad name
  def diskrep

    return @partial_hit, @mismatches_and_sequence_shifts if @partial_hit && @mismatches_and_sequence_shifts

    @mismatches_and_sequence_shifts = nil
    @partial_hit = nil
    #partial hit
    #or ### ?

    firstsplit = log.split(/###+/)
    firstsplit.each do |fspl|
      #    puts "splitted"
      # puts spl.to_s
      
      secondsplit = fspl.split(/---+/)
      secondsplit.each do |sspl|
        #     puts "second"
        if sspl.include?("partial hit")
          @partial_hit = sspl
          c = @partial_hit.scan(/\(#(\d*)\)/)
          #puts c.class
          c.each do |contignumber|
            # puts "loop"
            #  puts contignumber[0].to_i - 1
            # puts contigs.class
            index = contignumber[0].to_i - 1
            reason = contigs[index]["reason"]
            #puts reason
            @partial_hit = @partial_hit.sub(/^\(##{contignumber}\)/,"Contig: #{contignumber} \n #{reason} \n")
          end
        end
        if sspl.include?("mismatches and sequence shifts")
          @mismatches_and_sequence_shifts = sspl
        end
      end
    end
    #p = log.scan(/(^partial hits.*?)#}-{3,}|\#{3,}/)
    #puts contigs.size
      
    #discrepancy list
    #d = log.scan(/(^mismatches and sequence shifts.*?)###/m)
    #puts "runs here"
    return @partial_hit, @mismatches_and_sequence_shifts
  end


  def get_discrepancies_and_undetermined()
    discrep_partial_hits, discrep_mismatches_and_sequence_shifts = self.diskrep

    discrepancy_lines = []
    discrep_mismatches_and_sequence_shifts.split("\n").each_with_index do |line, line_index|
      if line_index <= 1 then
        discrepancy_lines << line
      else
        if line.match(/(\A\s+?)([0-9]+?)(:)/) then
          discrepancy_lines << line
        else
          discrepancy_lines[-1] += line
        end
      end
    end if discrep_mismatches_and_sequence_shifts


    discrepancies = {}
    undetermined = []

    discrepancy_lines.each_with_index do |line, line_index|

      # Ignore first and second lines
      next if line_index <= 1
      
      line.gsub!(/\s+|\n/, " ")

      ScipioResult::TYPES_OF_DISCREPANCY.each_pair do |type_number, type|
        escaped_type = type.gsub(/[\^\[\]\(\)\{\}\.\*\+\-\|\?\\\>\<]/){|match| "\\" + match}
        if line.match(/\A.+?[^0-9]#{type_number},.+?#{escaped_type}/) then
          aa = line.match(/(\A\s+?)([0-9]+?)(:)/)[2].to_i
          discrepancies[aa] = [] unless discrepancies[aa]
          discrepancies[aa] << type_number

          if type_number == 2 || type_number == 3 || type_number == 16 then
            undetermined << aa
          end

        end
      end

    end if discrep_mismatches_and_sequence_shifts

    return [discrepancies, undetermined]
  end

    
  def gff

    yaml_file = Tempfile.new("file_for_yaml2gff.yaml")
    #yaml_file_path = "#{TMP_DIR}/file_for_yaml2gff.yaml_#{Time.now.strftime("%Y%m%d_%H%M")}_#{rand(1000000000)}"

    begin

      exons.each do |exon|
        exon["undeterminedlist"] = [] unless exon["undeterminedlist"]
        exon["mismatchlist"] = [] unless exon["mismatchlist"]
      end
      
      yamlstring = scipio_yaml
      yamlstring = yamlstring.gsub(/exon_alternative/, 'exon')

      yaml_file_stream = yaml_file.open
      yaml_file_stream.puts yamlstring
      yaml_file.close

      gff_string = `#{YAML2GFF} #{yaml_file.path}`

    ensure
   
      yaml_file.close
      yaml_file.unlink
      FileUtils.remove_file(yaml_file.path, true)

    end

    return gff_string
    
  end

  def translation_fasta
    fasta = ">#{@name}_scipio_translation\n"
    fasta.concat(translation)
    return fasta
  end
    
  def fasta
    fasta = String.new
    exons.each_with_index do |exon,i|
      fasta.concat(">#{exon["type"]}_#{i+1}\n")
      #puts e.select("seq")[0].value
      fasta.concat(exon["seq"].upcase)
      fasta.concat("\n")
    end
    introns.each_with_index do |intron,j|
      fasta.concat(">intron_#{j+1}\n")
      fasta.concat(intron["seq"].upcase)
      fasta.concat("\n")
    end
    fasta.concat(">cdna\n")
    fasta.concat(self.cdna)
    fasta.concat("\n")
    fasta.concat(">genomic_dna\n")
    fasta.concat(self.genomic_dna)
    fasta.concat("\n")
    fasta.concat(">upstream\n")
    fasta.concat(self.upstream_dna.upcase) if self.upstream_dna
    fasta.concat("\n")
    fasta.concat(">downstream\n")
    fasta.concat(self.downstream_dna.upcase) if self.downstream_dna
    fasta.concat("\n")
    return fasta
  end
    
  #exonfasta
  def exonfasta
    fasta = String.new
    exons.each_with_index do |exon,i|
      fasta.concat(">exon_#{i+1}\n")
      #puts e.select("seq")[0].value
      fasta.concat(exon["seq"].upcase)
      fasta.concat("\n")
    end
    return fasta
  end

  def exon_has_alternatives?(exonnumber)
    has_alternatives = false
    # puts exonnumber
    exons.each do |ex|
      #  p ex
      if (ex["type"] == "exon_alternative") && (ex["exon"] == exonnumber || (ex["exon_tuple"] && ex["exon_tuple"].include?(exonnumber)))
        has_alternatives = true
      end
    end
    return has_alternatives
    
  end

  def exon_has_mutually_exclusives?(exonnumber)
    return mutually_exclusive_exons.any? do |ex| ex["exon"] == exonnumber end
  end

  def exon_has_duplicated_exons?(exonnumber)
    return duplicated_exons.any? do |ex| ex["exon"] == exonnumber end
  end

  def exon_get_alternatives(exonnumber)
    exons.find_all {|ex| (ex["type"] == "exon_alternative") && (ex["exon"] == exonnumber) }
  end

  def exon_get_mutually_exclusives(exonnumber)
    exons.find_all {|ex| (ex["type"] == "exon_alternative" && ex["alternative_type"] == "mutual_exclusive") && (ex["exon"] == exonnumber) }
  end
  
  def alt_exon_name(number)
    all =  exons.find_all { |exon| exon["type"] == "exon_alternative"  }.length
    before = alternative_exons.find_all { |exon| exon["dna_start"] < alternative_exons[number]["dna_start"] &&
        exon["exon"] == alternative_exons[number]["exon"]}.length
    in_casette = alternative_exons[number]["exon"].to_s + ('a' .. 'zzz').to_a[before + 1]
    # puts "#{alternative_exons[number - 1]["exon"]}" + in_casette
    #return in_casette
  end

  def exon_translation_fasta
    fasta = String.new
    #alt_counter = "b"
    exon_counter = 1
    exons.find_all{|ex| ex["type"] == "exon" || (ex["type"] == "exon_alternative" && ex["alternative_type"] != "mutual_exclusive")}.each do |exon|
      if exon_has_alternatives?(exon_counter)
        #   puts "Has alternatives"
        alt_counter = "a"
        #  puts exon_get_alternatives(exon_counter)
        exon_get_alternatives(exon_counter).find_all {|before| before["dna_start"] < exon["dna_start"] }.each do |before_alt|
          #   p before_alt
          fasta.concat(">alternative_exon_#{exon_counter}_#{alt_counter}\n")
          fasta.concat(before_alt["translation"].upcase)
          fasta.concat("\n")
          alt_counter = alt_counter.succ
        end
        fasta.concat(">#{exon["type"]}_#{exon_counter}_#{alt_counter}\n")
        fasta.concat(exon["translation"].upcase)
        fasta.concat("\n")
        alt_counter =  alt_counter.succ
        exon_get_alternatives(exon_counter).find_all {|after| after["dna_start"] > exon["dna_start"] }.each do |after_alt|
          #    p after_alt
          fasta.concat(">alternative_exon_#{exon_counter}_#{alt_counter}\n")
          fasta.concat(after_alt["translation"].upcase)
          fasta.concat("\n")
          alt_counter =  alt_counter.succ
        end
      else
        fasta.concat(">#{exon["type"]}_#{exon_counter}\n")
        fasta.concat(exon["translation"].upcase)
        fasta.concat("\n")
      end
      exon_counter += 1
    end

    return fasta
 
  end

  def exon_translation
    exon_translation = String.new
    exon_counter = 1
    exons.find_all {|ex| ex["type"] == "exon" }.each do |exon|
      if exon_has_alternatives?(exon_counter)
        alt_counter = "a"
        exon_get_alternatives(exon_counter).find_all {|before| before["dna_start"] < exon["dna_start"] }.each do |before_alt|
          exon_translation.concat(before_alt["translation"].upcase)
          alt_counter = alt_counter.succ
        end
        exon_translation.concat(exon["translation"].upcase)
        alt_counter =  alt_counter.succ
        exon_get_alternatives(exon_counter).find_all {|after| after["dna_start"] > exon["dna_start"] }.each do |after_alt|
          exon_translation.concat(after_alt["translation"].upcase)
          alt_counter =  alt_counter.succ
        end
      else
        exon_translation.concat(exon["translation"].upcase)
      end
      exon_counter += 1
    end

    return exon_translation

  end

  #intronfasta
  def intronfasta
    fasta = String.new
    introns.each_with_index do |intron,j|
      fasta.concat(">intron_#{j+1}\n")
      fasta.concat(intron["seq"].upcase)
      fasta.concat("\n")
    end
    return fasta
  end
    
  #codingdnafasta
  def cdnafasta
    fasta = String.new
    fasta.concat(">cdna\n")
    fasta.concat(self.cdna)
    fasta.concat("\n")
  end
    
  #genomicdnafast
  def genomicfasta
    fasta = String.new
    fasta.concat(">genomic_dna\n")
    fasta.concat(self.genomic_dna)
    fasta.concat("\n")
  end
    
  def gap_contig_lengths
    lens = []
    len = 0

    #debugger
    # puts"contigs:#{contigs.length}"
    #leading gap?
    e1 = contigs[0]["matchings"].clone.sort_by{|m| m["dna_start"]}[0]
    if e1["type"] == "exon" || e1["type"] == "exon_alternative" then
      if e1["gen_start"] then
        len = e1["gen_start"] - 0
    else
        len = e1["nucl_start"] - 0
      end
    end
    lens << len

    contigs[0..-2].each_index do |i|
      len = 0
      #  puts "i:#{i}"
      e1 = contigs[i]["matchings"].clone.sort_by{|m| m["dna_end"]}[-1]
      e2 = contigs[i + 1]["matchings"].clone.sort_by{|m| m["dna_start"]}[0]
      if (e1["type"] == "exon" || e1["type"] == "exon_alternative") && (e2["type"] == "exon" || e2["type"] == "exon_alternative") then
        len = e2["nucl_start"] - e1["nucl_end"]
      end
      lens << len
    end

    len = 0
    e1 = contigs[-1]["matchings"].clone.sort_by{|m| m["dna_end"]}[-1]
    if (e1["type"] == "exon" || e1["type"] == "exon_alternative") then
      if contigs[-1]["gen_len"] then
        len = contigs[-1]["gen_len"] - e1["gen_end"] if e1["gen_end"]
    else
        len = contigs[-1]["prot_len"]*3 - e1["nucl_end"] if contigs[-1]["prot_len"]
      end
    end
    lens << len #if len > 0

    return lens
  end


#   def gap_contig_lengths
#    lens = []
#    len = 0
#    # puts"contigs:#{contigs.length}"
#    #leading gap?
#    if (e1 = contigs[0]["matchings"][0])["type"] ==  "exon" then
#      len = e1["nucl_start"] - 0
#    end
#    lens << len
#    contigs[0..-2].each_index do |i|
#      #  puts "i:#{i}"
#      if (e1 = contigs[i]["matchings"][-1])["type"] ==  "exon" && (e2 = contigs[i + 1]["matchings"][0])["type"] ==  "exon" then
#        len =  e2["nucl_start"] - e1["nucl_end"]
#      end
#      lens << len #if len > 0
#    end
#    #puts "lens: #{lens}"
#    return lens
#  end

    
  def test_contig_length
    "Test contig length"
    contigs.each_with_index do |c,i|
      #   puts "contig_length #{i}:#{contig_lengths[i]}"
      seq_sum = 0#gap_contig_lengths[i]
      c["matchings"].each do |m|
        s = m["seq"]
        if(s != nil)
          #      puts m["type"]
          seq_sum += s.length
        end
      end
      # puts "seq_sum:#{seq_sum}"
    end
  end
    
  def reading_frames
    frames = [0]
    exons.each do |exon|
      # puts "exon"
      # puts exon["seq"].length / 3.0
      # puts exon["translation"].length
    end
  end
    
  def scale_int_ex_length(in_ex_ratio)
    scaled = exon_length + gap_contig_lengths.inject {|sum, n| sum + n }  + (upstream_gaps.inject {|sum, n| sum + n }  + intron_length  + gap_length + downstream_gaps.inject {|sum, n| sum + n } ) / in_ex_ratio
      
    return scaled
  end
    
  def intron_length
#    puts "DEBUG: intron_length: #{@name}: #{introns.select{|intron| intron["seq"].length != (intron["dna_end"]-intron["dna_start"]).abs}.collect{|intron| [intron["number"], intron["seq"].length - (intron["dna_end"]-intron["dna_start"]).abs]}.inspect}"
    introns.sum{|intron| intron["seq"] ? intron["seq"].length : 0}
 end


  def scipio_version
    version = ""
    @comment_lines.each do |line|
      if match = line.match(/(###\s)(.+)(\soutput)/) then
        version = match[2]
        break
      end
    end
    return version
  end


  def query_file
    file = ""
    @comment_lines.each do |line|
      if match = line.match(/(#\squery\sfile\s+)(.+)/) then
        file = match[2]
        break
      end
    end
    return file
  end


  def target_file
    file = ""
    @comment_lines.each do |line|
      if match = line.match(/(#\starget\sfile\s+)(.+)/) then
        file = match[2]
        break
      end
    end
    return file
  end


  def blat_output_file
  file = ""
    @comment_lines.each do |line|
      if match = line.match(/(#\sBLAT\soutput\s+)(.+)/) then
        file = match[2]
        break
      end
    end
    return file
  end

  
  def timestamp
    timestamp_string = ""
    @comment_lines.each do |line|
      if match = line.match(/(#\sTimestamp\s+)(.+)/) then
        timestamp_string = match[2]
        break
      end
    end

    timestamp = nil

    unless timestamp_string == ""

      begin
        timestamp = Time.parse(timestamp_string)
      rescue ArgumentError
      end
      
    end

    return timestamp
  end


  def make_color(exon_nr, type)
    @colors = HueRange.without_red
    @altcolors = Array.new
    @exoncolors = Hash.new

    self.contigs.each do |contig|
      contig["matchings"].each do |m|
        if m["type"] == "exon_alternative" then
          if m["alternative_type"] == "mutual_exclusive" || m["alternative_type"] == "duplicated_exon" then
            fillcolor = (self.count_alternative_exons.index(m["exon"]) + 1).to_f / self.count_alternative_exons.size.to_f
          else
            fillcolor =  (m["dna_start"] - contig["dna_start"]).to_f / (contig["dna_end"] - contig["dna_start"]).to_f
          end

          @altcolors.push(fillcolor)
          @exoncolors[m["exon"]] = @altcolors.last
          if m["exon_tuple"] then
            m["exon_tuple"].each do |number|
              @exoncolors[number] = @altcolors.last unless @exoncolors[number]
            end
          end

        end
      end
    end
    # p @altcolors
    #p @exoncolors
    #p @opacities
    #puts "exon:#{exon_nr} #{exon_has_alternatives(exon_nr +1,result)}"
    case  type
    when "exon"
      if self.exon_has_alternatives?(exon_nr) then
        #   puts exon_nr
        #  puts @exoncolors[exon_nr - 1]
        # p @colors
        color = @colors.map(@exoncolors[exon_nr])
      else
        color = "black"
      end
    when "alt_exon"
      #debugger
      color = @colors.map(@altcolors[exon_nr-1])
    end
    # puts "color:#{color}"
    return color
  end

  
  def to_s
    status = nil
    if perfect_match? then
      status = "complete"
    elsif perfect_match_but_mismatches_seqshifts?
      status = "complete_mm_fs"
    else
      status = "incomplete"
    end
    return name.ljust(25) + status.ljust(25) + (((gap_ratio * 1000).round / 10.0).to_s + "%").ljust(25) + (((nucl_match_ratio * 1000).round / 10.0).to_s + "%").ljust(25)  + contigs.length.to_s
  end
    
end

# Injection to handle DNA stretches
class Range
  
  # Returns true if the ranges overlap
  def overlap?(other)
    ol = false
    #|------|
    #   X------|
    # if self.member?(other.begin) || self.member?(other.end) then ol = true end
    #  |------|
    #|------X
    if other.member?(self.begin) || other.member?(self.end) then ol = true end
    return ol
  end
  
  # Compare by begin
  def <=>(other)
    self.begin <=> other.begin
  end
  
  # Returns a new range if the ranges overlap or are adjacent.
  # The resulting range is reverse if both input ranges are reverse.
  # Returns false if the input ranges have a gap between them.
  def +(other)
    if self.grow(1).overlap?(other) then
      # determin if concat is reverse
      rev = (self.begin > self.end && other.begin > other.end)
      b = [self.begin, self.end, other.begin, other.end].min
      e = [self.begin, self.end, other.begin, other.end].max
      if rev then b, e = e, b end
      ret = (b..e)
    else ret = false end
    return ret
  end

  # Expands the edges of the range by 'by'
  def grow(by)
    if    self.begin < self.end then g = ((self.begin - by)..(self.end + by))
    elsif self.begin > self.end then g = ((self.begin + by)..(self.end - by))
    else g = ((self.begin - by)..(self.begin + by)) end
    return g
  end
  
  
end
