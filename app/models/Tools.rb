class Tools

end

class String
	# Extends class string by an naturalized sort method
	# Sorting like finder: by name and number in ascending order
	def naturalized
		scan(/[^\d\.]+|[\d\.]+/).collect { |f| f.match(/\d+(\.\d+)?/) ? f.to_f : f }
	end
end
