class ReferenceData
	attr_reader :path, :data

	def initialize(file)
		@path = file
		@data = load_ref_data
	end

	# separate by class & delete unseparated data & store back in ref_data
	def split_prot_into_classes(prot)
		prot_fam_obj = ProteinFamily.new(prot,@data[prot])	

		# first: delete unseparated data
		del_prot(prot)

		# second: separate this obj into classes
		keys = prot_fam_obj.ref_alignment.keys & prot_fam_obj.ref_genes.keys
		keys.each do |name|

			prot_abbr = Split_prot_abbrs[prot]

			klass = ""
			new_prot_key = ""

			matches = name.match(/#{prot_abbr}([0-9]+|Mhc)/)

			# Myosin protein has two abbreviations: Mhc for class 2 and Myo for all other classes
			if ! matches && name.include?("Mhc") then
				matches = ["", "Mhc"]
			end

			if matches && matches[1] then
				klass = matches[1]
				new_prot_key = prot + " Class " + klass
			else
				new_prot_key = prot
			end

			if ! @data.has_key?(new_prot_key) then
				add_prot(new_prot_key)
			end

			# third: save genes and alignments
			if prot_fam_obj.ref_genes.has_key?(name)
				@data[new_prot_key]["genes"][name] = prot_fam_obj.ref_genes[name]
			end

			if prot_fam_obj.ref_alignment.has_key?(name)
				@data[new_prot_key]["alignment"] += Helper::Sequence.str2fasta(name,prot_fam_obj.ref_alignment[name])
				@data[new_prot_key]["alignment"] += "\n"
			end

		end
	end

	def add_prot(key)
		@data[key] = {}
		@data[key]["genes"] = {}
		@data[key]["alignment"] = ""
	end

	def del_prot(key)
		@data.delete(key)
	end

	def update_alignment(key, fasta_hash)
		@data[key]["alignment"] = fasta_hash.map{|k,v| ">#{k}\n#{v}"}.join("\n") + "\n"
	end

	def update_genes(key, genes_hash)
		@data[key]["genes"] = genes_hash
	end

	def create_protfam_obj(key)
		obj = ProteinFamily.new(key,@data[key])	
	end

	def load_ref_data
		Helper.file_exist_or_die(@path)
		return JSON.load(File.read(@path)) 
	end

	def save_ref_data(this_path=@path)
		fh = File.new(this_path, "w")
		fh.puts JSON.dump(data)
		fh.close
	end
end