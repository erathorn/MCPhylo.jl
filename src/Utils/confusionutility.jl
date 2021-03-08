"""
    generate_confusion_matrix(pairwise_results)::Array{Int64,2}

Generates a simple confusion matrix.

Confusion matrix is returned as an Array of Int values.

* `pairwise_results` : Array containing arrays of float values; for each internal array, the float value at index 1 should represent the expected value, and index 2, the returned value. Both float values should be EITHER 1.0 OR 0.0.
"""
function generate_confusion_matrix(pairwise_results)::Array{Int64,2}

    confusion = zeros(Int64,2,2)
    for pair in pairwise_results
      pair[1] == 1.0 ? xaxis = 1 : xaxis = 2
      if xaxis == 1
        pair[2] == 1.0 ? yaxis = 1 : yaxis = 2
      else
        pair[2] == 0.0 ? yaxis = 2 : yaxis = 1
      end #ifelse
      confusion[xaxis,yaxis]+= 1
    end #for
    return confusion
end #generate_confusion_matrix
"""
    get_summary_statistics(confusion::Array{Int64,2})::Tuple{Float64,Float64,Float64}

Calculates accuracy, precision, and fscore for a confusion matrix.

Returns those values, in that order, as a tuple of floats.

* `confusion` : the matrix on which to perform the calculations; ensure the matrix is formatted such that true positives occur at index (1,1), true negatives at (2,2), false negatives at (1,2), and false positives at (2,1).
"""
function get_summary_statistics(confusion::Array{Int64,2})::Tuple{Float64,Float64,Float64}
  #confusion = generate_confusion_matrix(pairwise_results)

  truepos = confusion[1,1]
  trueneg = confusion[2,2]
  falseneg = confusion[1,2]
  falsepos = confusion[2,1]

  accuracy = (truepos + trueneg)/length(pairwise_results)
  precision = truepos/(confusion[1,1]+confusion[1,2])
  fscore = truepos/(truepos + (0.5*(falsepos + falseneg)))

  accuracy,precision,fscore
end #get_summary_statistics
