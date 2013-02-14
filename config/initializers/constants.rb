#host depending settings
if ENV && ENV["PWD"] && (ENV["PWD"].include?("fab8") || ENV["PWD"].include?("jenkins")) then
	LUCULLUS_URL = "http://fab8:8080/tpl_os"
else
	# FIXME
	# replace "bagheera" by correct URL & and add it to proxy (same as cymo, peakr)
	LUCULLUS_URL = "http://www.bagheera.org/tpl_os"
end
