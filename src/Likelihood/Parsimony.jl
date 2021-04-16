
"""
      parsimony(tree::N, char::Dict{String,String}, gap::String="-")::Float64 where N<: GeneralNode

Do parsimony reconstruction for a tree and a set of characters. The characters are
supplied as the char dictionary with the key being the name of a leave and the value
string representation of the character.
"""
function parsimony(tree<:GeneralNode, char::Dict{String,String}, gap::String="-")::Float64
      po = post_order(tree)
      states = filter(x -> x != "-", unique(collect(values(char))))
      nStates = length(states)
      if nStates == 0
            return 0.0
      end
      nNodes = length(po)
      costMatrix = zeros(nNodes, nStates) .- 1
      for node in po
            if node.nchild == 0
                  # i am a leave node
                  s = char[node.name]
                  costMatrix[node.num, :] .= s == gap ? 0 : abs.(log.(states .== s))
            else
                  # i am an internal node
                  daughters = node.children
                  localCosts = zeros(nStates, nStates, length(daughters))
                  for (j, s1) in enumerate(states),
                        (k, s2) in enumerate(states),
                        (l, d) in enumerate(daughters)

                        localCosts[j, k, l] = costMatrix[d.num, k] + Int(s1 != s2)
                  end
                  costMatrix[i, :] = vec(
                        mapslices(
                              sum,
                              mapslices(minimum, localCosts, dims = 2),
                              dims = 3,
                        ),
                  )
            end
      end
      minimum(costMatrix[tree.num, :])
end
