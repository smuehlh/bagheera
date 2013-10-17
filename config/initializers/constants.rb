# BASE_PATH = Dir::tmpdir + "/cug/" # all data files are stored in /tmp
BASE_PATH = Rails.root.join('db', 'ref_data').to_s
REF_DATA = "alignment_gene_structure.json"
PATH_REF_WO_EXAMPLE = "without_ca_b/"

# FORMATDB = "/usr/bin/formatdb"
# BLASTALL = "/usr/bin/blastall"
# FASTACMD = "/usr/bin/fastacmd"
FORMATDB = "/usr/local/ncbi_2.2.28+/bin/makeblastdb"
BLASTALL = "/usr/local/ncbi_2.2.28+/bin/tblastn"
FASTACMD = "/usr/local/ncbi_2.2.28+/bin/blastdbcmd"
Legacy_blast = "/usr/local/ncbi_2.2.28+/bin/legacy_blast.pl"
New_blast_path = "/usr/local/ncbi_2.2.28+/bin"

AUGUSTUS = "/usr/local/bin/augustus/src/augustus"
# DIALIGN2 = "/usr/local/bin/dialign_package/src/dialign2-2"
FastTree = "/usr/local/bin/FastTree"

PAIR_ALIGN =  Rails.root.join('lib', 'pair_align').to_s #"/usr/local/bin/pair_align"
Progress_file = Rails.root.join('public', 'status.html').to_s

#host depending settings
if ENV && ENV["PWD"] && (ENV["PWD"].include?("fab8") || ENV["PWD"].include?("jenkins")) then
        LUCULLUS_URL = "http://fab8:8080/tpl_os"
        GBLOCKS = "/usr/local/Gblocks_0.91b/Gblocks"
        MAFFT = "/usr/bin/mafft"
else
        LUCULLUS_URL = "http://www.motorprotein.de/tpl_os"
        GBLOCKS = "/usr/local/bin/Gblocks_0.91b/Gblocks"
        MAFFT = "/usr/local/bin/mafft"
end
