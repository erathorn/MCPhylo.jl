

function calc_b_cubed(true_labels, labels)
    """
    Calculate the B-cubed (precision, recall, F-score) of a list of cognate set
    labels against the gold-standard version of the same list.
    This function is a (just slightly) modified version of the b_cubed function
    of PhyloStar's CogDetect library.
    """
    precision = zeros(length(true_labels))
    recall = zeros(length(true_labels))

    for (i, l) in enumerate(labels)
        match = 0.0
        prec_denom = 0.0
        recall_denom = 0.0
        for (j, m) in enumerate(labels)
            if l == m
                prec_denom += 1.0
                if true_labels[i] == true_labels[j]
                    match += 1.0
                    recall_denom += 1.0
                end
            elseif l != m
                if true_labels[i] == true_labels[j]
                    recall_denom += 1.0
                end
            end
        end
        precision[i] = match / prec_denom
        recall[i] = match / recall_denom
    end
    avg_precision = mean(precision)
    avg_recall = mean(recall)
    avg_f_score = 2.0 * (avg_precision * avg_recall) / (avg_precision + avg_recall)

    return avg_precision, avg_recall, avg_f_score
end

function calc_f_score(true_clusters, pred_clusters)
    """
    Calculate the B-cubed F-score of a dataset's cognate sets against their
    gold-standard. This is the metric used to evaluate the performance of the
    cognacy identification algorithms.
    Both args should be dicts mapping concepts to frozen sets of frozen sets of
    Word named tuples. The first comprises the gold-standard clustering and the
    second comprises the inferred one.
    It is assumed that both clusterings comprise the same data and that there
    is at most one word per concept per doculect. An AssertionError is raised
    if these assumptions do not hold true.
    """
    f_scores::Array{Float64} = []
    prec::Array{Float64} = []
    rec::Array{Float64} = []

    for concept in keys(true_clusters)

        #@assert concept in keys(pred_clusters)

        true_labels = Dict()
        pred2 = Dict()

        for (index, cog_set) in enumerate(true_clusters[concept])
            label = "true:"*string(index)

            for lang in true_clusters[concept][cog_set[1]]
                true_labels[lang] = label
            end
        end

        for (index, cog_set) in enumerate(pred_clusters[concept])
            label = "pred:"*string(index)
            for lang in pred_clusters[concept][cog_set[1]]
                pred2[lang] = label
            end
        end

        pred_labels = Dict()
        for (k, v) in pred2
            if k in keys(true_labels)
                pred_labels[k] = v
            end
        end
        @assert Set(keys(true_labels)) == Set(keys(pred_labels))
        sorted_keys = sort([i for i in keys(true_labels)])
        true_labels_list = [true_labels[label] for label in sorted_keys]
        pred_labels_list = [pred_labels[label] for label in sorted_keys]
        p, r, f = calc_b_cubed(true_labels_list, pred_labels_list)
        push!(f_scores, f)
        push!(prec, p)
        push!(rec, r)

    end

    return prec, rec, f_scores
end
