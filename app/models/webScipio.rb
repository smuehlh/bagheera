module WebScipio

	extend self
	require "net/http"
	require "uri"
	require 'yaml'

	Scipio_default_params = {:blattile => 7, :minid => 90, :maxmis => 7, :min_score => 0.3, :exhaust_align_size => 15000}
	Scipio_less_stringent_params = {:blattile => 7, :minid => 30, :maxmis => 0, :min_score => 0.1, :exhaust_align_size => 15000}

	def call_webscipio(species, protein_seq, is_use_scipio_default_params)

		yaml = ""

		err_msg = catch(:problem) {
		# worked_or_throw_error

			# find all availabe genomes
			is_success, genomes = is_species_valid(species)
			Helper.worked_or_throw_error(is_success, "No genome file found for species #{species}")

			# find 'best' genome
			best_genome = select_best_genome( genomes )
			if best_genome.nil? || best_genome[:path].nil? || best_genome[:path].empty? then 
				Helper.worked_or_throw_error( false, "Cannot select a genome for species #{species}." )
			end

			# scipio-search in best genome
			status, scipio_result = run_scipio( best_genome, species, protein_seq, is_use_scipio_default_params )

			# process results
			is_success = is_scipio_run_successful(status)
			Helper.worked_or_throw_error( is_success, "Gene prediction failed. Nothing found.<br> Try it again with relaxed parameters." )

			yaml =  YAML.dump(scipio_result)
		}

		if ! yaml.blank? then 
			err_msg = ""
		end

		return yaml, err_msg

	rescue NoMethodError, TypeError, NameError, RuntimeError, Errno::ENOENT, Errno::ETIMEDOUT
		Helper.raise_runtime_error("An error occured with Webscipio")
	end

	def is_species_valid(species_query)
		url = URI.parse("http://#{Webscipio_server}:#{Webscipio_port}/api_searches")
		post_parameters = {'search_genomes' => 'true', 'query' => species_query}
		response = Net::HTTP.post_form(url, post_parameters)
		id = response.body

		url = URI.parse("http://#{Webscipio_server}:#{Webscipio_port}/api_searches/#{id}.yaml")
		response = Net::HTTP.get_response(url)
		yaml_string = response.body
		genomes = YAML::load(yaml_string)

		return genomes.any?, genomes
	end

	def select_best_genome(genomes)
		if genomes.size == 1 then
			# it's just one genome, simply use it! 
			return genomes.first
		else
			# select for: type, version, size
			# sort genomes by reverse order: worst genome is listed first, best last

			# in reverse order: from worst to best case
			genome_assembly_types = ["reads", "contigs", "supercontigs", "ultracontigs", "chromosome"]
			sorted_genomes = genomes.sort_by do |hash|
				ind_assembly_type =  genome_assembly_types.index(hash[:type]) || 0 # if its not found in type, assume worst type
				version = [ hash[:major_version] || 0,
					hash[:minor_version] || 0, 
					hash[:mini_version] || 0] # assume worst version if its not specified
				size = hash[:size] || 0 # assume worst size if its not speciified

				# sort by type, version, size
				[ind_assembly_type, version, size]
			end

			return sorted_genomes.last
		end
	end

	def run_scipio(genome, species, fasta, is_use_scipio_default_params )
		url = URI.parse("http://#{Webscipio_server}:#{Webscipio_port}/api_searches")
		post_parameters = {'scipio_run' => 'true', 'target_file_path' => genome[:path], 'query' => fasta}
		if is_use_scipio_default_params then 
			post_parameters.merge!(Scipio_default_params)
		else
			post_parameters.merge!(Scipio_less_stringent_params)
		end
		response = Net::HTTP.post_form(url, post_parameters)
		id = response.body

		url = URI.parse("http://#{Webscipio_server}:#{Webscipio_port}/api_searches/#{id}.yaml")
		scipio_status_result = ["running", ""] # scipio returns "running" until its done!
		while( scipio_status_result[0] == "running" ) do 
			response = Net::HTTP.get_response(url)
			yaml_string = response.body
			scipio_status_result = YAML::load(yaml_string)
			sleep(5)
		end

		scipio_status = scipio_status_result[0]
		scipio_result = YAML::load(scipio_status_result[1])

		return scipio_status, scipio_result

	end

	def is_scipio_run_successful(scipio_status)
		return scipio_status == "finished"
	end

end