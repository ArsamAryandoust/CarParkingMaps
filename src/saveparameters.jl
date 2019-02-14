### saving parameters

function saveparameters(N, path_to_results, T, number_zones, cars_per_zone, C, e, p_min, p_max, A_set, A_drive, △A)

    # preparing the parameter values for saving
    parameters = zeros(N, 10);
    for i = 1:N
        parameters[i,1] = T;
        parameters[i,2] = number_zones;
        parameters[i,3] = cars_per_zone;
        parameters[i,4] = C;
        parameters[i,5] = e;
        parameters[i,6] = p_min[i];
        parameters[i,7] = p_max[i];
        parameters[i,8] = A_set;
        parameters[i,9] = A_drive[i];
        parameters[i,10] = △A[i];
    end
    parameters = convert(DataFrame, parameters);

    # create the header for the parameter file and save the results
    header_vector = ["T (time steps)", "number_zones", "cars_per_zone", "C (number of cars)", "e (probability exponent)", "p_min", "p_max", "A_set", "A_drive", "△A" ];
    results_string = string(path_to_results,"/results_parameters.csv");
    CSV.write(results_string, parameters, header=header_vector);

end