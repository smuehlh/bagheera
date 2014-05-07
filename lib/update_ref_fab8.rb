#!/usr/local/bin/ruby

### Update Bagheeras reference data on fab8
### Called monthly in crontab

# 12 3 1 * * ruby update_ref_neu.rb > /tmp/cug/cron-ruby.log 2> /tmp/cug/err-ruby.log
require 'timeout'
require 'fileutils'
require 'json'
require 'open3'
# require 'ruby-debug'

# load helper classes and modules
Dir[File.join(File.dirname(__FILE__), 'update_ref_fab8_lib', '*.rb')].each {|file| require File.absolute_path(file) }

### files and pathes to the new refdata and their final location
Path = "/tmp/new_cug"
Json = "alignment_gene_structure.json"
tmp_path_new_data = ""
final_pathes_new_data = ["/fab8/server/bagheera/db/ref_data", "/fab8/smuehlh/rails_projects/bagheera/db/ref_data"]
subfolder = "without_ca_b" # should be expanded to tmp_path_new_data

### the final and temporary pathes for new data must exist and be writable
if File.writable?(Path) then
	tmp_path_new_data = Path
else
	# a fresh dir under /tmp should be writable!
	tmp_path_new_data = "/tmp/new_cug2" 
	Helper.mkdir_or_die(tmp_path_new_data)
end
Helper.mkdir_or_die( File.join(tmp_path_new_data, subfolder) )

final_pathes_new_data.delete_if do |path|
	# del if both are NOT true
	! ( Helper.does_dir_exist(path) && File.writable?(path) )
end

### path to new json file
path_to_native_json = File.join(Path, Json) # this comes directly from cymo-API
path_to_revised_json = File.join(tmp_path_new_data, Json) # this is the revised json (adapted to actual needs ... )

### details about how to handle protein families
Del_unusual_prots = ["Calcineurin", "Calmodulin", "Centrin", "Frequenin", 
	"Dynein Light Intermediate Chain",
	"Myosin heavy chain", "Myosin essential light chain", "Myosin regulatory light chain", "Myosin heavy chain Class 17",
	"Kinesin", "Kinesin Class 4", "Kinesin Class 16"]
Split_prot_families = ["Myosin heavy chain", "Actin related protein", "Kinesin", "Tubulin", "Capping Protein"]
Split_prot_abbrs = { "Myosin heavy chain" => "Myo", 
	"Actin related protein" => "Arp", 
	"Kinesin" => "Kinesin", 
	"Tubulin" => "Tub", 
	"Capping Protein" => "CAP" }


### define process housekeeping
at_exit {
	# clean up
	Helper.del_file_or_dir( tmp_path_new_data )
}
### define function
def callCymoAPI(temp)
	wget_path = `which wget`.chomp
	system(wget_path, "--spider", "http://fab8:2001/api_cug_alignment/all", "-o", File.join(temp, "cron.log"))
end

def validate_save_alignment_and_precalc_profile(ref_data_obj, prot, prot_obj, tmp_path)
	# 1) adapt alignments: ensure same lenght of all seqs and remove common gaps
	prot_obj.ensure_same_length
	prot_obj.remove_common_gaps

	# 2) save alignments to file and to reference data

	f_out = File.join(tmp_path, "#{prot_obj.prot_filesave_name}.fasta")	
	Helper::Sequence.save_alignment(f_out, prot_obj.ref_alignment)
	ref_data_obj.update_alignment(prot, prot_obj.ref_alignment)
	ref_data_obj.update_genes(prot, prot_obj.ref_genes)

	# 3) test if mafft can handle alignment files
	prot_obj.ensure_mafft_is_fine(f_out)

	# 4) precalculate protein profiles
	f_prfl = f_out.sub("fasta", "prfl")
	Helper.calc_protein_profile(f_out, f_prfl)

end

### "main script"
# prepare directory structure

Helper.del_file_or_dir(path_to_native_json)

# fetch reference data from cymobase
# max_secs = 60*10 # wait max. 10 minutes for cymoapi
puts "Start wget: #{Time.now}"
begin
	status = Timeout::timeout(max_secs) { callCymoAPI(tmp_path_new_data) }
rescue Timeout::Error => exc
	Helper.abort(exc)
end
Helper.worked_or_die(status && Helper.does_file_exist(path_to_native_json), "wget failed")

ref_data_obj = ReferenceData.new(path_to_native_json)

# split protein families into protein classes
Split_prot_families.each do |fam|
	ref_data_obj.split_prot_into_classes(fam)
end

# delete unusual proteins 
Del_unusual_prots.each do |prot| 
	ref_data_obj.del_prot(prot)
end
ref_data_obj.save_ref_data(path_to_revised_json)

# 'reload' ref_data to have new set of keys in ref_data_obj
ref_data_obj = ReferenceData.new(path_to_revised_json)
ref_data_wo_cab_obj = ReferenceData.new(path_to_revised_json)

prot_list = ref_data_obj.data.keys
prot_list.sort.each do |prot|
	puts prot

	# complete set of reference data

	prot_obj = ref_data_obj.create_protfam_obj(prot)
	tmp_path = tmp_path_new_data

	validate_save_alignment_and_precalc_profile(ref_data_obj, prot, prot_obj, tmp_path)

	# set of reference data without example species 'Ca_b'
	prot_obj = nil
	prot_obj = ref_data_wo_cab_obj.create_protfam_obj(prot)
	# remove all 'Ca_b' genes from genes and alignments
	prot_obj.ref_alignment.delete_if { |k,v| k =~ /Ca_b/ }
	prot_obj.ref_genes.delete_if { |k,v| k =~ /Ca_b/ }

	tmp_path = File.join(tmp_path_new_data, subfolder)

	validate_save_alignment_and_precalc_profile(ref_data_wo_cab_obj, prot, prot_obj, tmp_path)

end

puts "Save reference data to file"
ref_data_obj.save_ref_data(path_to_revised_json)
ref_data_wo_cab_obj.save_ref_data( File.join( tmp_path_new_data, subfolder, Json ) )

puts "Move everything to place"
final_pathes_new_data.each do |path|

	# delete old stuff in path, files and folders
	Helper.del_file_or_dir( Dir.glob(File.join(path, "*.*")) ) # delete old fasta, prfl, json and log files
	Helper.del_file_or_dir( File.join( path, subfolder ) )

	# move new data to place
	scr = tmp_path_new_data + '/.'
	FileUtils.cp_r scr, path
end
