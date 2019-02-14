### Saving results

function saveresults(number_zones, state_matrix, path_to_results, data_set, C)

    # prepare the resulting data as share of total cars in each zone
    results = zeros(number_zones,T);
    for i=1:C
        for t=1:T
            index = state_matrix[i,t];
            results[index,t] = results[index,t] + 1;
        end
    end
    results = results/C;
    results = convert(DataFrame, results);

    # create the header for the CSV file and save the results
    header_vector = fill("", T); # creates an empty string vector in the desired dimension
    for t = 1:T
        header_vector[t] = string("t = ", t, "h");
    end
    results_string = string(path_to_results,"/results_",data_set);
    CSV.write(results_string, results, header=header_vector);

end