### Saving results

function saveresults(number_zones, state_matrix, transition_matrix, path_to_results, data_set, C)

    # prepare the resulting location of cars, the share of parking and driving cars in each zone
    driving_cars = zeros(number_zones, T);
    parking_cars = zeros(number_zones, T);
    for i=1:C
        for t=1:T
            index = state_matrix[i,t];
            driving_activity = transition_matrix[i,t,1];
            parking_cars[index,t] += 1;
            if(driving_activity == 1)
                driving_cars[index,t] += 1;
            end
        end
    end

    # normalize with the total number of cars
    parking_cars = parking_cars/C;

    # calculating circadian rythm of driving activity
    traffic_resultsmatrix = sum(driving_cars,dims=1);
    min_sampled = minimum(traffic_resultsmatrix); 
    max_sampled = maximum(traffic_resultsmatrix); 
    for i=1:24
        traffic_resultsmatrix[i] = (traffic_resultsmatrix[i]-min_sampled)/(max_sampled - min_sampled);
    end

    # convert to DataFrame for saving as CSV files
    parking_cars = convert(DataFrame, parking_cars);
    traffic_resultsmatrix = convert(DataFrame, traffic_resultsmatrix);
    
    # create the header for the CSV file and save the results
    header_vector = fill("", T); # creates an empty string vector in the desired dimension
    for t = 1:T
        header_vector[t] = string("t = ", t, "h");
    end
    results_string_parking = string(path_to_results,"/results_parkingdensities_",data_set);
    results_string_driving = string(path_to_results,"/results_trafficactivity_",data_set);
    CSV.write(results_string_parking, parking_cars, header=header_vector);
    CSV.write(results_string_driving, traffic_resultsmatrix, header=header_vector);
    

end