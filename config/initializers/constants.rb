BASE_PATH = Dir::tmpdir + "/cug/" # all data files are stored in /tmp
REF_DATA = "alignment_gene_structure.json"
#host depending settings
if ENV && ENV["PWD"] && (ENV["PWD"].include?("fab8") || ENV["PWD"].include?("jenkins")) then
	LUCULLUS_URL = "http://fab8:8080/tpl_os"
	Verbose_error = true
else
	# FIXME
	# replace "bagheera" by correct URL & and add it to proxy (same as cymo, peakr)
	LUCULLUS_URL = "http://www.bagheera.org/tpl_os"
	Verbose_error = false
end
