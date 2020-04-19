
mutable struct Alignment <: DiscreteMatrixDistribution
        emp::AbstractArray{Float64,2}
        px::AbstractArray{Float64}
        py::AbstractArray{Float64}
		a::Float64
		r::Float64
		timemat::AbstractArray{Float64}
        langs::Int64
        concs::Int64
		wsize::Int64
end

Base.minimum(d::Alignment) = -Inf
Base.maximum(d::Alignment) = +Inf
Base.size(d::Alignment) = (d.langs, d.concs, d.wsize)

function logpdf(d::Alignment, x::AbstractArray{T, 3})::Float64 where T<:Real
	calcl_ll(convert.(Int64,x), d.emp, d.px, d.py, d.a, d.r, d.timemat)
end

function vl(d::Alignment, x::AbstractArray{T, 3}) where T<:Real
	calc_viterbi(convert.(Int64,x), d.emp, d.px, d.py, d.a, d.r, d.timemat)
end

function cognate_statements(d::Alignment, x::AbstractArray{T,3}; quant::Float64=0.9) where T <: Real
	res = calc_viterbi(convert.(Int64,x), d.emp, d.px, d.py, d.a, d.r, d.timemat)
	out = Array{Bool, 3}(undef, d.langs, d.langs, d.concs)
	for l1=1:d.langs, l2=1:l1
		if l1 != l2
			mv = res[l1,l2,:,:]
			dv = mv[diagind(mv)]
			mv[diagind(mv)] .= -Inf
			nmv = mv[mv .!= -Inf]
			out[l1, l2, :] .= dv .> quantile(nmv, quant)
		end
	end
	out
end





function forward_log(s1::AbstractArray{Int64}, s2::AbstractArray{Int64}, m::Int64, n::Int64, em_p::AbstractArray{Float64,2},
					gx_p::AbstractArray{Float64}, gy_p::AbstractArray{Float64}, trans_p::AbstractArray{Float64})::Float64

	P = 0.0
	trellis = Array{Float64,3}(undef, m+1, n+1,3)
	trellis .= -Inf

    t00 = trans_p[1]
    t01 = trans_p[2]
    t02 = trans_p[3]
    t0E = trans_p[4]
    t10 = trans_p[5]
    t11 = trans_p[6]
    t12 = trans_p[7]
    t1E = trans_p[8]
    t20 = trans_p[9]
    t21 = trans_p[10]
    t22 = trans_p[11]
    t2E = trans_p[12]
    tS0 = trans_p[13]
    tS1 = trans_p[14]
    tS2 = trans_p[15]

    trellis[1,1,1] = tS0
	trellis[1,1,2] = tS1
	trellis[1,1,3] = tS2


	trellis[2,1,2] = gx_p[s1[1]] + log(exp(trellis[1,1,1]+t01) + exp(trellis[1,1,2]+t11)+ exp(trellis[1,1,2]+t21))


	for i = 3:m + 1
	    trellis[i,1,2] = gx_p[s1[i-1]] + trellis[i - 1,1,2] + t11
	end



	trellis[1,2,3] = gy_p[s2[1]] + log(exp(trellis[1,1,1]+t02) + exp(trellis[1,1,3]+t22)+ exp(trellis[1,1,2]+t12))
	for i = 3:n + 1
		trellis[1,i,3] = gy_p[s2[i-1]] + trellis[1,i - 1,3] + t22
	end

	for i = 2:m+1
	    x = s1[i - 1]
	    xv = gx_p[x]
		for j = 2:n+1
		    y = s2[j - 1]
		    yv = gy_p[y]

		    v = em_p[x,y]

			trellis[i,j,1] = v + log(exp(t00 + trellis[i - 1,j - 1,1]) + exp(t10+trellis[i - 1,j - 1,2]) + exp(t20 + trellis[i - 1,j - 1,3]))
			trellis[i,j,2] = xv + log(exp(trellis[i - 1,j,1] + t01) + exp(trellis[i - 1,j,2] + t11) + exp(trellis[i - 1,j,3] + t21))
			trellis[i,j,3] = yv + log(exp(trellis[i,j - 1,1] + t02) + exp(trellis[i,j - 1,3] + t22) + exp(trellis[i,j - 1,2] + t12))

		end
	end
	P = log(exp(t0E + trellis[m,n,1]) + exp(t1E + trellis[m,n,2]) + exp(t2E + trellis[m,n,3]))
	return P
end


function viterbi_log(s1::AbstractArray{Int64}, s2::AbstractArray{Int64}, m::Int64, n::Int64, em_p::AbstractArray{Float64,2},
					gx_p::AbstractArray{Float64}, gy_p::AbstractArray{Float64}, trans_p::AbstractArray{Float64})::Float64

	P = 0.0
	trellis = Array{Float64,3}(undef, m+1, n+1,3)
	trellis .= -Inf

    t00 = trans_p[1]
    t01 = trans_p[2]
    t02 = trans_p[3]
    t0E = trans_p[4]
    t10 = trans_p[5]
    t11 = trans_p[6]
    t12 = trans_p[7]
    t1E = trans_p[8]
    t20 = trans_p[9]
    t21 = trans_p[10]
    t22 = trans_p[11]
    t2E = trans_p[12]
    tS0 = trans_p[13]
    tS1 = trans_p[14]
    tS2 = trans_p[15]

    trellis[1,1,1] = tS0
	trellis[1,1,2] = tS1
	trellis[1,1,3] = tS2


	trellis[2,1,2] = gx_p[s1[1]] + max(trellis[1,1,1]+t01,trellis[1,1,2]+t11,trellis[1,1,2]+t21)


	for i = 3:m + 1
	    trellis[i,1,2] = gx_p[s1[i-1]] + trellis[i - 1,1,2] + t11
	end



	trellis[1,2,3] = gy_p[s2[1]] + max(trellis[1,1,1]+t02,trellis[1,1,3]+t22,trellis[1,1,2]+t12)
	for i = 3:n + 1
		trellis[1,i,3] = gy_p[s2[i-1]] + trellis[1,i - 1,3] + t22
	end

	for i = 2:m+1
	    x = s1[i - 1]
	    xv = gx_p[x]
		for j = 2:n+1
		    y = s2[j - 1]
		    yv = gy_p[y]

		    v = em_p[x,y]

			trellis[i,j,1] = v + max(t00 + trellis[i - 1,j - 1,1], t10+trellis[i - 1,j - 1,2],t20 + trellis[i - 1,j - 1,3])
			trellis[i,j,2] = xv + max(trellis[i - 1,j,1] + t01, trellis[i - 1,j,2] + t11, trellis[i - 1,j,3] + t21)
			trellis[i,j,3] = yv + max(trellis[i,j - 1,1] + t02, trellis[i,j - 1,3] + t22, trellis[i,j - 1,2] + t12)

		end
	end
	P = max(t0E + trellis[m,n,1],t1E + trellis[m,n,2],t2E + trellis[m,n,3])
	return P
end



function random_model(s1, s2, m, n, gx, gy)::Float64
	score = 0.0
	@inbounds for i=1:m
		score += gx[i]
	end
	@inbounds for i=1:n
		score += gy[i]
	end
	score
end

function read_cognate_data(filename)
	open(filename, "r") do file
		global content = readlines(file)
	end # do

	header = popfirst!(content)
	header = split(header, "\t")
	index = zeros(Bool,10)
	index[1] = 1
	index[2] = 1
	index[8] = 1
	header[index]
	data = []
	while length(content) != 0
		line = popfirst!(content)
		line = split(line, "\t")[index]
		push!(data, line)
	end
	words = [i[3] for i in data]
	langs = Dict(key=>value for (value, key) in enumerate(String.(unique([i[1] for i in data]))))
	concs = Dict(key=>value for (value, key) in enumerate(String.(unique([i[2] for i in data]))))
	maximum(length.(words))
	alphabet = Dict(key=> value for (value, key) in enumerate(unique([j for i in words for j in i])))
	alldatamat = Array{Int64,3}(undef, length(langs), length(concs), 16)
	alldatamat .= -1

	for entry in data
		lang, conc, word = entry
		indlang = langs[lang]
		indconc = concs[conc]
		l = 0
		for (indl,letter) in enumerate(word)
			alldatamat[indlang, indconc, indl]=alphabet[letter]
			l+=1
		end
		alldatamat[indlang, indconc, 16] = l
	end
	return alldatamat, alphabet
end



function calcl_ll(datamat::Array{Int64,3}, submat::AbstractArray{Float64,2}, gx::AbstractArray{Float64}, gy::AbstractArray{Float64},
					a::Float64, r::Float64, timemat::AbstractArray{Float64})::Float64
	l_size, c_size, w_size = size(datamat)
	score = Base.Threads.Atomic{Float64}(0)
	@inbounds @views Base.Threads.@threads for l1=1:l_size
		for l2=1:l1
				if l1 != l2
					t = timemat[l1, l2]
					emmat = log.(exp(submat*t))
					tra = tr_td(r, t, a)
					for c=1:c_size
						w1 = datamat[l1, c, :]
						w2 = datamat[l2, c, :]
						if w1[1] != -1 && w2[1] != -1
							m = w1[end]
							n = w2[end]
							Base.Threads.atomic_add!(score, forward_log(w1, w2, m, n, emmat, gx, gy, tra))
						end
				end
			end
		end
	end
	score[]
end


function calc_viterbi(datamat::Array{Int64,3}, submat::AbstractArray{Float64,2}, gx::AbstractArray{Float64}, gy::AbstractArray{Float64},
					a::Float64, r::Float64, timemat::AbstractArray{Float64})
	l_size, c_size, w_size = size(datamat)
	scores = zeros(l_size, l_size, c_size, c_size)
	scores .= -Inf
	@inbounds @views Base.Threads.@threads for l1=1:l_size
		for l2=1:l_size
				if l1 != l2
					t = timemat[l1, l2]
					emmat = log.(exp(submat*t))
					tra = tr_td(r, t, a)
					for c=1:c_size
						w1 = datamat[l1, c, :]
						for c2=1:c_size
							w2 = datamat[l2, c2, :]
							if w1[1] != -1 && w2[1] != -1
								m = w1[end]
								n = w2[end]
								scores[l1, l2, c, c2] = viterbi_log(w1, w2, m, n, emmat, gx, gy, tra) - random_model(w1, w2, m, n, gx, gy)
							end
					end
				end
			end
		end
	end
	scores
end



function calcl_ll(datamat::Array{Int64,3}, submat::AbstractArray{Float64,2}, gx::AbstractArray{Float64}, gy::AbstractArray{Float64}, tra::AbstractArray{Float64})::Float64
	l_size, c_size, w_size = size(datamat)
	score = Base.Threads.Atomic{Float64}(0)
	@inbounds @views Base.Threads.@threads for c=1:c_size
		for l1=1:l_size
			for l2=1:l1
				if l1 != l2

					w1 = datamat[l1, c, :]
					w2 = datamat[l2, c, :]
					m = w1[end]
					n = w2[end]
					if m != -1 && n != -1

						 Base.Threads.atomic_add!(score, forward_log(w1, w2, m, n, submat, gx, gy, tra))
						#score += forward_log(w1, w2, m, n, submat, gx, gy, tra)#-random_model(w1, w2, m, n, gx, gy)
					end
				end
			end
		end
	end
	score[]
end


function tr_td(r::Float64, t::Float64, a::Float64)
    indel::Float64 = 1.0 - exp(-2.0*r*t)
    v1::Float64 = log(1.0-indel)
    v2::Float64 = log(0.5*indel)
    v3::Float64 = log((1.0-a)*(1-indel))
    v4::Float64 = log((1.0-a)*indel)

	out = zeros(15)
    out[1] = v1 #t00
    out[2] = v2 #t01
    out[3] = v2#t02
    out[4] = v1#t0E
    out[5] = v3#t10
    out[6] = log(a)#t11
    out[7] = v4#t12
    out[8] = v3#t1E
    out[9] = v3#t20
    out[10] = v4#t21
    out[11] =log(a)#t22
    out[12] =v3#t2E
    out[13] =v1#tS0
    out[14] =v2#tS1
    out[15] =v2#tS2
	return out
end

function maptrprobs(trpobs1, trprobs2)
	δ, τm, rm1 = trpobs1
	ϵ, τxy, λ , rm2 = trprobs2
	tr = zeros(15)
	t00 = 1-2*δ-τm
	tr[1] = rm1
	tr[2] = δ
	tr[3] = δ
	tr[4] = τm
	tr[5] = rm2
	tr[6] = ϵ
	tr[7] = λ
	tr[8] = τxy
	tr[9] = rm2
	tr[10] = ϵ
	tr[11] = λ
	tr[12] = τxy
	tr[13] = rm1
	tr[14] = δ
	tr[15] = δ
	return tr

end


function indices2values(maparr::Array{Int64,2}, vec::AbstractArray{Float64}, freqs::AbstractArray{Float64}, nchar::Int64)::Array{Float64}
	x = zeros(Float64, size(maparr))
	@inbounds @simd for i=1:nchar
		for j=1:i
			if i!=j
				x[j,i] = vec[maparr[j,i]]*freqs[i]
				x[i,j] = vec[maparr[i,j]]*freqs[j]
			end
		end
	end
	di = -sum(x, dims = 2)
	x[diagind(nchar,nchar)] .= di[:]
	x * -1.0/dot(freqs, di)
end



function project(mtime::AbstractArray{Float64}, sl::Int64)
	tm = zeros(sl, sl)
	k = 1
	@inbounds for i = 1:sl
		for j = 1:i
			if i != j
				tm[i,j] = mtime[k]
				tm[j,i] = mtime[k]
				k+=1
			end
		end
	end
	tm
end


"""
    to_distance_matrix(tree::T)::Array{Float64,2} where T <:AbstractNode

Calculate the distance matrix over the set of leaves.
"""
function mto_distance_matrix_i(tree::T)::Array{Float64,2} where T <:MCPhylo.AbstractNode
    leaves::Vector{T} = MCPhylo.get_leaves(tree)
	sort!(leaves, by=y->parse(Int64, y.name))
    ll = length(leaves)
    distance_mat = zeros(Float64, ll, ll)
    for i in 1:ll
        for j in 1:i
            if i!=j
                d = MCPhylo.node_distance(tree, leaves[i], leaves[j])
                distance_mat[i,j] = d
                distance_mat[j,i] = d
            end # if
        end # for
    end #for
    distance_mat
end # function to_distance_matrix

function mto_distance_matrix(tree::T)::Array{Float64,2} where T <:MCPhylo.TreeVariate
    mto_distance_matrix_i(tree.value)
end # function to_distance_matrix
