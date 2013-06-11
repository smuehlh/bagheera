class Tools 

end

class String
	# Extends class string by an naturalized sort method
	# Sorting like finder: by name and number in ascending order
	def naturalized
		scan(/[^\d\.]+|[\d\.]+/).collect { |f| f.match(/\d+(\.\d+)?/) ? f.to_f : f }
	end

end

class Array
  def find_each_index find
	found, index, q = -1, -1, []
	while found
	  found = self[index+1..-1].index(find)
	  if found
		index = index + found + 1
		q << index
	  end
	end
	q
  end
end
