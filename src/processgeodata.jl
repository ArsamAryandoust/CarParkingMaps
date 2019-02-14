### Processing GEOJSON file: INPUT is path_to_json_data and OUTPUT is distance_matrix_km, number_zones

function processgeodata(path_to_json_data, path_to_data, csv_dataset_list, path_to_results)

    open(path_to_json_data, "r") do json_file 
        global julia_file
        julia_file = JSON.parse(json_file);
    end

    a = julia_file["features"]
    first_zone_id = parse(Int,a[1]["properties"]["MOVEMENT_ID"]);
    if (first_zone_id == 0) # catches the case that the first zone ID is a zero
        number_zones = parse(Int,a[length(a)]["properties"]["MOVEMENT_ID"]) + 1;
    else
        number_zones =  parse(Int,a[length(a)]["properties"]["MOVEMENT_ID"]);
    end
    longitude_matrix = zeros(number_zones, 100000);
    latitude_matrix = zeros(number_zones, 100000);

    for i=1:length(a)
        entry = a[i];
        id = parse(Int,entry["properties"]["MOVEMENT_ID"]);
        if (id == 0) # catches the case that the first zone ID is a zero
            id = number_zones;
        end
        coordinates = entry["geometry"]["coordinates"];
        counter = 0;
        
        if(typeof(coordinates) == Float64)
            counter = counter + 1;
            longitude_matrix[id, counter] = coordinates[1];
            latitude_matrix[id, counter] = coordinates[2];
        else
            for j=1:length(coordinates)
                nested_coordinates = coordinates[j]
                if(typeof(nested_coordinates) == Float64)
                    counter = counter + 1;
                    longitude_matrix[id, counter] = coordinates[1];
                    latitude_matrix[id, counter] = coordinates[2];
                else
                    for k=1:length(nested_coordinates)
                        nested_coordinates = coordinates[j][k];
                        if(typeof(nested_coordinates) == Float64)
                            counter = counter + 1;
                            longitude_matrix[id, counter] = coordinates[j][1];
                            latitude_matrix[id, counter] = coordinates[j][2];
                        else
                            for l=1:length(nested_coordinates)
                                nested_coordinates = coordinates[j][k][l];
                                if(typeof(nested_coordinates) == Float64)
                                    counter = counter + 1;
                                    longitude_matrix[id, counter] = coordinates[j][k][1];
                                    latitude_matrix[id, counter] = coordinates[j][k][2];
                                else
                                    for m=1:length(nested_coordinates)
                                        nested_coordinates = coordinates[j][k][l][m];
                                        if(typeof(nested_coordinates) == Float64)
                                            counter = counter + 1;
                                            longitude_matrix[id, counter] = coordinates[j][k][l][1];
                                            latitude_matrix[id, counter] = coordinates[j][k][l][2];
                                        else
                                            for n=1:length(nested_coordinates)
                                                nested_coordinates = coordinates[j][k][l][m][n];
                                                if(typeof(nested_coordinates) == Float64)
                                                    counter = counter + 1;
                                                    longitude_matrix[id, counter] = coordinates[j][k][l][m][1];
                                                    latitude_matrix[id, counter] = coordinates[j][k][l][m][2];
                                                else
                                                    for o=1:length(nested_coordinates)
                                                        nested_coordinates = coordinates[j][k][l][m][n][o];
                                                        if(typeof(coordinates) == Float64)
                                                            counter = counter + 1;
                                                            longitude_matrix[id, counter] = coordinates[j][k][l][m][n][1];
                                                            latitude_matrix[id, counter] = coordinates[j][k][l][m][n][2];
                                                        else
                                                            println("Double-check how often this json file is nested at:")
                                                            println("i = ", i)
                                                            println("j = ", j)
                                                            println("k = ", k)
                                                            println("l = ", l)
                                                            println("m = ", m)
                                                            println("n = ", n)
                                                            println("o = ", o)
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    ## averaging coordinates for each zone 
    average_lat = zeros(number_zones,1);
    average_long = zeros(number_zones,1);
    summing_lat = 0;
    summing_long = 0;
    for i=1:number_zones
        counter = 0;
        summing_lat = 0;
        summing_long = 0;
        for j=1:length(longitude_matrix[1,:])
            if(longitude_matrix[i,j] == 0)
                break;
            else
                counter = counter + 1;
                summing_lat = summing_lat + latitude_matrix[i,j];
                summing_long = summing_long + longitude_matrix[i,j];
            end
        end
        average_lat[i, 1] = summing_lat / counter;
        average_long[i, 1] = summing_long / counter;
    end

    ## calculating geographic distances between zones 
    c_lat_long = 111.3; # this is a natural constant 
    conv_deg_rad = 0.01745; # this is a constant for the conversion from degrees into radians. It serves for input into trig. functions. 
    distance_matrix_km = zeros(number_zones, number_zones);
    for i=1:(number_zones-1)
        long1 = average_long[i];
        lat1 = average_lat[i];
        for j=i:number_zones
            long2 = average_long[j];
            lat2 = average_lat[j];
            distance_km = c_lat_long * sqrt((cos((lat1+lat2)/2 * conv_deg_rad))^2 * (long1-long2)^2 + (lat1-lat2)^2);
            distance_matrix_km[i,j] = distance_km;
            distance_matrix_km[j,i] = distance_km;
        end 
    end

    ## save the results
    results = zeros(number_zones, 2)
    for i = 1:number_zones
        results[i,1] = average_lat[i];
        results[i,2] = average_long[i];
    end
    results = convert(DataFrame, results);
    header_vector = ["latitude", "longitude"]; # creates an empty string vector in the desired dimension
    results_string = string(path_to_results,"/zoneID_coordinates.csv");
    CSV.write(results_string, results, header=header_vector);

    return distance_matrix_km, number_zones

end