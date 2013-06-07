module Yaml2Log

  YAML2LOG_COMMAND = "#{Rails.root.to_s}/script/yaml2log.pl --width=85"

  @@pipes = Queue.new # start 5 scripts (perhaps more, depending on traffic)
  5.times {@@pipes << IO.popen(YAML2LOG_COMMAND, "w+")}

  def Yaml2Log.logfile(yamlstring)

    pipe = @@pipes.pop      # retrieve one instance of a script
    
    pipe.write(yamlstring)    # write the data
    pipe.write("\n#END\n")
    pipe.write("\n")
    pipe.flush
    
    result = ""

    while true              # read the result
      begin
        string = pipe.readpartial(4096)
        result += string
        break if /#END/.match result
      rescue EOFError
        pipe = IO.popen(YAML2LOG_COMMAND, "w+")
        break
      end
    end
  
    @@pipes << pipe   # push back the pipe instance

    return result.gsub("#END", "")
    
  end

end
