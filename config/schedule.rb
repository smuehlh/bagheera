every :sunday, :at => "4am" do

	command "rm /tmp/cug/cron.log && rm /tmp/cug/err.log"

	# delete reference data file
	command "rm /tmp/cug/new_alignment_gene_structure.json"

	# download newest data from cymobase, be not verbose, and not quiet
	command "wget --spider http://fab8:2001/api_cug_alignment/all -a '/tmp/cug/cron.log'"

	# prepare the data:
	# => separate actin, myosin and kinesin by class
	# => make all alignments accectable by mafft (they are not if they are not of same length)
	# this task will wait for the wget to be finished
	runner "CronjobController.prepare_ref_data", :output => '/tmp/cug/cron.log'
	
	# => delete data files older than one day
	runner "CronjobController.delete_old_data"

end