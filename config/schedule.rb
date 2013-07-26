every :sunday, :at => "4am" do

	command "rm -f /tmp/cug/cron.log && rm -f /tmp/cug/cron_ruby.log && rm -f /tmp/cug/err.log"

	# delete reference data file
	command "rm -f /tmp/cug/new_alignment_gene_structure.json"

	# download newest data from cymobase
	command "wget --spider http://fab8:2001/api_cug_alignment/all -a '/tmp/cug/cron.log'"

	# prepare the data:
	# => separate actin, myosin and kinesin by class
	# => make all alignments accectable by mafft (they are not if they are not of same length)
	# this task will wait for the wget to be finished
	runner "CronjobController.prepare_ref_data", :output => '/tmp/cug/cron_ruby.log'
	
	# => delete data files older than one day
	runner "CronjobController.delete_old_data"

end

every :sunday, :at => "5am" do
	# copy reference data to cymo
	runner "CronjobController.copy_refdata_to_cymo", :output => '/tmp/cug/cron_copytocymo.log'
end

