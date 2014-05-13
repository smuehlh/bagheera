BASE_PATH = Rails.root.join('db', 'ref_data').to_s
BASE_PATH_PROTEIN = File.join(BASE_PATH, "protein")
BASE_PATH_TRNA = File.join(BASE_PATH, "trna")
REF_DATA = "alignment_gene_structure.json"
TRNA_REF_DATA = "tRNAserCAG_wo_intron.fas"
PATH_REF_WO_EXAMPLE = "without_ca_b/"
Tmp_path = File.join(Dir::tmpdir, "cug/")

Max_filesize = 26214400 # 25 MB

# FORMATDB = "/usr/bin/formatdb"
# BLASTALL = "/usr/bin/blastall"
# FASTACMD = "/usr/bin/fastacmd"

AUGUSTUS = "/usr/local/bin/augustus/src/augustus"
# DIALIGN2 = "/usr/local/bin/dialign_package/src/dialign2-2"
FastTree = "/usr/local/bin/FastTree"
MAFFT = "/usr/local/bin/mafft"

PAIR_ALIGN =  Rails.root.join('lib', 'pair_align').to_s #"/usr/local/bin/pair_align"

TRNASCAN = "tRNAscan-SE"
# Progress_file = Rails.root.join('public', 'status.html').to_s

#host depending settings
if ENV && ENV["PWD"] && (ENV["PWD"].include?("fab8") || ENV["PWD"].include?("jenkins")) then
        LUCULLUS_URL = "http://fab8:8080/tpl_os"
        GBLOCKS = "/usr/local/Gblocks_0.91b/Gblocks"
        New_blast_path = "/usr/local/ncbi_2.2.28+/bin"
else
        LUCULLUS_URL = "http://www.motorprotein.de/tpl_os"
        GBLOCKS = "/usr/local/bin/Gblocks_0.91b/Gblocks"
        New_blast_path = "/usr/local/ncbi-blast-2.2.28+/bin"
end
        