every :sunday, :at => "12pm" do
	# download newest data from cymobase
	command "wget http://fab8:2001/api_cug_alignment/all -o '/tmp/cug/cron.log'"

	# prepare the data:
	# => separate actin, myosin and kinesin by class
	# => make all alignments accectable by mafft (they are not if they are not of same length)
	runner "CronjobController.prepare_ref_data", :output => '/tmp/cug/cron.log'
	
	# => delete data files older than one day
	runner "CronjobController.delete_old_data"
end