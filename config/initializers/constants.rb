BASE_PATH = Dir::tmpdir + "/cug/" # all data files are stored in /tmp
REF_DATA = "alignment_gene_structure.json"

FORMATDB = "/usr/bin/formatdb"
BLASTALL = "/usr/bin/blastall"
FASTACMD = "/usr/bin/fastacmd"
AUGUSTUS = "/usr/local/bin/augustus/src/augustus"
DIALIGN2 = "/usr/local/bin/dialign_package/src/dialign2-2"
MAFFT = "/usr/local/bin/mafft"

PAIR_ALIGN =  Rails.root.join('lib', 'pair_align').to_s #"/usr/local/bin/pair_align"

#host depending settings
if ENV && ENV["PWD"] && (ENV["PWD"].include?("fab8") || ENV["PWD"].include?("jenkins")) then
	LUCULLUS_URL = "http://fab8:8080/tpl_os"
else
	LUCULLUS_URL = "http://www.motorprotein.de/bagheera/tpl_os"
end
