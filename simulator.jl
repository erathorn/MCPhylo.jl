
function sample_data(tree, eq_f, rates, subst_model, nsites)
    pro = pre_order(tree)
    nnodes = length(pro)
    sim_data = zeros(Int,nnodes, nsites)
    U, D, Uinv, mu = subst_model(eq_f, rates)
    out_data = zeros(Float64, length(eq_f), nsites, nnodes)
    for node in pro
        if node.root
            sim_data[node.num, : ] .= rand(Categorical(eq_f), nsites)
        else

            t = node.inc_length
            trans_mat = BLAS.gemm('N', 'N', 1.0, BLAS.symm('R', 'L', diagm(exp.(D .* (mu*t))), U), Uinv)
            for site in 1:nsites
                sampled_val = rand(Categorical(trans_mat[sim_data[node.mother.num, site],:]))
                sim_data[node.num, site] = sampled_val
                if node.nchild == 0
                    out_data[sampled_val, site, node.num] = 1
                end
            end
            
        end
    end
    out_data
end


function max_psrf(sim::ModelChains)::Float64
    
    gd = gelmandiag(sim)
    gd_values = gd.value
    indices = isnan.(gd_values)
    gd_values[indices] .= -1
    psrf = maximum(gd_values[:,1])
    psrf
end
