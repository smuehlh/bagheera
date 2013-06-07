
class PrettyAligner

  def initialize
    
  end

  def get_alignment(alt_exon)

    @alt_exon = alt_exon
    
    @pretty_cand = ""
    @pretty_al = ""
    @fill_line = ""

    @pretty_dna = @alt_exon["seq"].upcase
    @dna_count = @alt_exon["dna_start"]
    @aa_count = @alt_exon["nucl_start"] / 3

    @frame = @alt_exon["reverse"] ? (@alt_exon["seq"].length+@alt_exon["nucl_end"]) % 3 : @alt_exon["nucl_start"] % 3

    @gaps = []

    similarity_hash = SubstitutionMatrix.similarity()

    gap_pos = 0
    k=0

    @alt_exon["cand_al"].each_byte do |a|
      @pretty_cand.concat(" " + a.chr + " ")
      if a.chr == "-"
        @fill_line.concat(" - ")
        @gaps.push gap_pos
      else
        if a.chr != @alt_exon["exon_al"][k].chr then
          if similarity_hash[a.chr + @alt_exon["exon_al"][k].chr] then
            @fill_line.concat(" + ")
          else
            @fill_line.concat(" X ")
          end
        else
          @fill_line.concat(" | ")
        end
      end
      gap_pos += 1
      k +=1
    end

    @gaps.each do |gap|
      case @frame
      when 0
        if @pretty_dna.size < gap * 3 then
          @pretty_dna.concat("---")
        else
          @pretty_dna.insert(gap * 3,"---")
        end
      when 1
        if @pretty_dna.size < gap * 3 - 1 then
          @pretty_dna.concat("---")
        else
          @pretty_dna.insert(gap * 3 - 1,"---")
        end
      when 2
        if @pretty_dna.size < gap * 3 + 1 then
          @pretty_dna.concat("---")
        else
          @pretty_dna.insert(gap * 3 + 1,"---")
        end
      end
    end

    @alt_exon["exon_al"].each_byte{|a| @pretty_al.concat(" " + a.chr + " ") }

    @pretty_dna.lstrip!
    @pretty_cand.lstrip!
    @fill_line.lstrip!
    @pretty_al.lstrip!

    if @frame == 1 then
      @pretty_dna = " " + @pretty_dna
      @pretty_cand = " " + @pretty_cand
      @fill_line = " " + @fill_line
      @pretty_al = " " + @pretty_al
    elsif @frame == 2 then
      @pretty_dna = "  " + @pretty_dna
      @pretty_cand = "    " + @pretty_cand
      @fill_line = "    " + @fill_line
      @pretty_al = "    " + @pretty_al
    else
      @pretty_dna = "" + @pretty_dna
      @pretty_cand = " " + @pretty_cand
      @fill_line = " " + @fill_line
      @pretty_al = " " + @pretty_al
    end

    #changes to_array!!!
    @pretty_dna = cut_dna_lines(@pretty_dna, @dna_count)
    @pretty_cand = cut_aa_lines(@pretty_cand,nil)
    @fill_line = cut_aa_lines(@fill_line,nil)
    @pretty_al = cut_aa_lines(@pretty_al,@aa_count)

    al = String.new

    # Fixing bug if @pretty_al is shorter than @pretty_dna
    @pretty_cand << " "*(@pretty_dna.size-@pretty_cand.size) if @pretty_dna.size > @pretty_cand.size
    @fill_line << " "*(@pretty_dna.size-@fill_line.size) if @pretty_dna.size > @fill_line.size
    @pretty_al << " "*(@pretty_dna.size-@pretty_al.size) if @pretty_dna.size > @pretty_al.size
 
    @pretty_dna.each_with_index do |line,i|
      al << "  " + line + "\n"
      al << "  " + (@pretty_cand[i] || " ") + "\n"
      al << "  " + (@fill_line[i] || " ") + "\n"
      al << "  " + (@pretty_al[i] || " ") + "\n"
      al << "\n"
    end
    return "#{al}"
  end

  
  def get_alignment_cdna(exon)

    return "" unless exon && exon["seq"]

    exon_dna_line = exon["seq"]
    exon_trans_line = String.new
    fill_line = String.new
    cdna_trans_line = String.new

    exon_trans_sequence = exon["translation"].to_s.upcase


    if exon["cdnaseq"] then
      cdna_line = exon["cdnaseq"].to_s.upcase
    else
      cdna_line = "N" * (exon["nucl_end"] - exon["nucl_start"]).abs
    end

    total_shift_genomic_dna = 0
    total_shift_cdna = 0

    exon["seqshifts"].each do |seqshift|

      cdna_nucl_diff = (seqshift["nucl_end"] - seqshift["nucl_start"]).abs
      genomic_dna_diff = (seqshift["dna_end"] - seqshift["dna_start"]).abs

      shift = cdna_nucl_diff - genomic_dna_diff

      if shift > 0
        position_in_exon = (seqshift["dna_start"] - exon["dna_start"]).abs + total_shift_genomic_dna
        exon_dna_line.insert(position_in_exon, "-"*shift) if exon_dna_line.size > position_in_exon
        total_shift_genomic_dna += shift.abs
      else
        position_in_exon = (seqshift["nucl_start"] - exon["cdna_start"]).abs + total_shift_cdna
        cdna_line.insert(position_in_exon, "-"*-shift) if cdna_line.size > position_in_exon
        total_shift_cdna += shift.abs
      end

   end

    if exon["protseq"] then
      cdna_trans_sequence = exon["protseq"].to_s.upcase
    else
      cdna_trans_sequence = "X" * ((exon["nucl_end"] - exon["nucl_start"]).abs/3)
#      cdna_trans_sequence = Bio::Sequence::AA.new(hit_data["prot_seq"]).subseq(exon["prot_start"]+1, exon["prot_end"]).to_s.upcase
    end



    genomic_dna_offset = exon["dna_start"]

    # Strand "minus"
    if exon["nucl_start"] < 0 || exon["nucl_end"] < 0 then
      # TODO:
      #cdna_offset = cdna_len + exon["nucl_end"]
      cdna_offset = 0
    else
      cdna_offset = exon["nucl_start"]
    end

    
    exon_trans_sequence.each_byte{|a| exon_trans_line.concat(" " + a.chr + " ") }
    exon_trans_sequence.each_byte{|a| fill_line.concat(" | ") }
    cdna_trans_sequence.each_byte{|a| cdna_trans_line.concat(" " + a.chr + " ") }

    #frameshift
    frame = exon["frameshift"]

    gaps = []
    index = -1

    exon_dna_line.each_char do |n|

      index += 1

      if n == "-" then
        exon_trans_line.insert(index, " ")
        fill_line.insert(index, " ")
        cdna_trans_line.insert(index, " ")
        cdna_line[index] = cdna_line[index].chr.downcase if cdna_line[index]

        gaps << index
      end
    end




    exon_dna_line.strip!
    exon_trans_line.strip!
    fill_line.strip!
    cdna_trans_line.strip!
    cdna_line.strip!


    if frame == 1 then
#      exon_dna_line = " " + exon_dna_line
      exon_trans_line = "  " + exon_trans_line
      fill_line = "  " + fill_line
      cdna_trans_line = "  " + cdna_trans_line
#      cdna_line = " " + cdna_line
    elsif frame == 2 then
#      exon_dna_line = "  " + exon_dna_line
      exon_trans_line = "" + exon_trans_line
      fill_line = "" + fill_line
      cdna_trans_line = "" + cdna_trans_line
#      cdna_line = "  " + cdna_line
    else
#      exon_dna_line = "" + exon_dna_line
      exon_trans_line = " " + exon_trans_line
      fill_line = " " + fill_line
      cdna_trans_line = " " + cdna_trans_line
#      cdna_line = "" + cdna_line
    end


    index = -1
    fill_line.each_char do |a|

      index += 1

      if a == "|" then
        first_index = [0, index-1].max
        second_index = [index+1, exon_dna_line.length-1, cdna_line.length-1].min

        mismatch = (exon_dna_line[first_index..second_index].upcase != cdna_line[first_index..second_index].upcase) if second_index >= first_index
        silent_mutation = (exon_trans_line[index] == cdna_trans_line[index])

        if(mismatch)
          if(silent_mutation)
            fill_line[index] = "!"
          else
            fill_line[index] = "X"
          end
        end
      end

    end
    
    #changes to_array!!!
    exon_dna_line = cut_dna_lines(exon_dna_line,genomic_dna_offset)
    exon_trans_line = cut_aa_lines(exon_trans_line,nil)
    fill_line = cut_aa_lines(fill_line,nil)
    cdna_trans_line = cut_aa_lines(cdna_trans_line,nil)
    cdna_line = cut_dna_lines(cdna_line, cdna_offset)

    output = ""

    cdna_line.each_with_index do |line,i|
      output << "  " << (exon_dna_line[i] || ' ') << "\n"
      output << "  " << (exon_trans_line[i] || ' ') << "\n"
      output << "  " << (fill_line[i] || ' ') << "\n"
      output << "  " << (cdna_trans_line[i] || ' ') << "\n"
      output << "  " << cdna_line[i] << "\n"
      output << "\n"
    end


    return "#{output}"
  end

  
  private



  def cut_dna_lines(dna, dna_start = 0, cut_length = 72)

    splitted_dnas = []
    cutted_dna_start = dna_start

    (0..((dna.size.to_i-1) / cut_length.to_i)).each do |split_index|
      
      cutted_dna = dna[split_index*cut_length, cut_length]
      cutted_dna_end = cutted_dna_start + cutted_dna.scan(/[A-Z]/i).size
      fill = " " * (11 - cutted_dna_end.abs.to_s.size) + " " * (cut_length - cutted_dna.length)
      cutted_dna_start = cutted_dna_end 

      splitted_dnas << cutted_dna + fill + cutted_dna_end.abs.to_s
        
    end

    return splitted_dnas

  end

  # If aa_start is nil a "|" will be shown at the end instead of the aminoacid position
  def cut_aa_lines(aa, aa_start = nil, cut_length = 72)


    splitted_aas = []
    cutted_aa_start = aa_start

    (0..((aa.size.to_i-1) / cut_length.to_i)).each do |split_index|

      cutted_aa = aa[split_index*cut_length, cut_length]

      # If aa_start is nil a "|" will be shown at the end instead of the aminoacid position
      if aa_start then

        cutted_aa_end = cutted_aa_start + cutted_aa.scan(/[A-Z]/i).size if aa_start
        fill = " " * (11 - cutted_aa_end.abs.to_s.size) + " " * (cut_length - cutted_aa.length)
        cutted_aa_start = cutted_aa_end

        splitted_aas << cutted_aa + fill + cutted_aa_end.abs.to_s

      else

        fill = " " * 10 + " " * (cut_length - cutted_aa.length)
        

        splitted_aas << cutted_aa + fill + "|"
        
      end

    end

    return splitted_aas
  end

  
end

