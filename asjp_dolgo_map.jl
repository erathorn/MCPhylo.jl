mapping_classes = Dict{Int64,Vector{Char}}(

    1=> ['a', 'e', 'E', '3', 'u', 'o', 'i'],
    2=>['c', 'C', 'T', 'k', 'g', 'x', 'q', 'G'],
    3=>['p', 'b', 'f'],
    4=>['h'],
    5=>['j', 'y'],
    6=>['m'],
    7=>['n', '5', 'N'],
    8=>['s', 'z', 'S', 'Z'],
    9=>['l', 'L', 'r'],
    10=>['t', 'd', '8'],
    11=>['w', 'v']
)

mapping = Dict{Int64,Char}(0=> 'a', 1=> 'e', 2=> 'E', 3=> '3', 4=> 'u', 5=> 'o', 6=> 'i', 7=> 'c', 8=> 'C', 9=> 'T', 10=> 'k', 11=> 'g', 12=> 'x',
           13=> 'q', 14=> 'G', 15=> 'p', 16=> 'b', 17=> 'f', 18=> 'h', 19=> 'j', 20=> 'y', 21=> 'm', 22=> 'n', 23=> '5', 24=> 'N',
           25=> 's', 26=> 'z', 27=> 'S', 28=> 'Z', 29=> 'l', 30=> 'L', 31=> 'r', 32=> 't', 33=> 'd', 34=> '8', 35=> 'w', 36=> 'v')



function letter2class(mapdict)
	retdict = Dict{Char,Int64}()
	for (key, value) in pairs(mapdict)
		for letter in value
			retdict[letter] = key
		end
	end
	retdict
end


function get_classes(maping)
	rmap::Vector{Tuple{Int64,Int64}} = []
	for (i,j) in Iterators.product(keys(mapping_classes),keys(mapping_classes))
		if !((i,j) in rmap || (j,i) in rmap)
			push!(rmap, (i,j))
		end
	end
	Dict(itm=>ind for (ind,itm) in enumerate(rmap))
end


function set_map(i2l, l2c, rmap)
	mymap::Vector{Int64} = []
	arr=Array{Int64,2}(undef, 37,37)
	for (index, letter) in pairs(i2l)
		for (index2, letter2) in pairs(i2l)
			class1 = l2c[letter]
			class2 = l2c[letter2]
			if class1 == class2
				push!(mymap, rmap[(-1, -1)])
				v = rmap[(-1, -1)]
			elseif class1 == 1 || class2 == 1
				push!(mymap, rmap[(-2, -2)])
				v = rmap[(-2, -2)]
			else
				if (class1, class2) in keys(rmap)
					v = rmap[(class1, class2)]
				else
					v = rmap[(class2, class1)]
				end

				push!(mymap, v)
			end
			arr[index, index2] = v
			arr[index2, index] = v
		end
	end
	mymap, arr
end
