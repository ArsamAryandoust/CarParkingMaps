### saving parameters

function saveparameters(path_to_results, T, number_zones, cars_per_zone, C, e_drive, p_min, p_max, e_dest, A_drive)

    # preparing the parameter values for saving
    parameters = zeros(1, 9);
    
    parameters[1,1] = T;
    parameters[1,2] = number_zones;
    parameters[1,3] = cars_per_zone;
    parameters[1,4] = C;
    parameters[1,5] = e_drive;
    parameters[1,6] = p_min;
    parameters[1,7] = p_max;
    parameters[1,8] = e_dest;
    parameters[1,9] = A_drive;

    parameters = convert(DataFrame, parameters);

    # create the header for the parameter file and save the results
    header_vector = ["T (time steps)", "number_zones", "cars_per_zone", "C (number of cars)", "e_drive (model parameter 1)", "p_min (model parameter 2)", "p_max (model parameter 3)", "e_dest (model parameter 4)", "A_drive" ];
    results_string = string(path_to_results,"/sampling_parameters.csv");
    CSV.write(results_string, parameters, header=header_vector);

end