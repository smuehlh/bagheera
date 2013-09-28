#!/usr/local/bin/ruby

### Clean up Bagheeras temporary files
### Called dayly in crontab

# 00 23 * * * ruby /fab8/server/bagheera/lib/del_old_data.rb

class Fixnum
	SECONDS_IN_DAY = 24 * 60 * 60
	def days
		self * SECONDS_IN_DAY
	end
	def ago
		Time.now - self
		end
end

max_days = 2

# delete all session directories (only numbers) older than 2 days
path = File.join("/tmp/cug","**[0-9]")
Dir.glob(path) do |dir|
	if File.mtime(dir) <= max_days.days.ago 
		begin
			FileUtils.rm_rf(dir)
		rescue => err
		end
	end
end

# delete also alignments for lucullus
path = File.join("/tmp","cymobase_alignment_cug*")
Dir.glob(path) do |file|
	if File.mtime(file) <= max_days.days.ago
		begin
			File.rm(file)
		rescue => err
		end
	end
end

# may not always work?
# find "/tmp/cug" -maxdepth 1 -type d -mtime +2 -regextype sed -regex ".*/[0-9]*" -exec rm -rf {} \;
