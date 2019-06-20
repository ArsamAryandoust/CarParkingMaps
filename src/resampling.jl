### Resampling: INPUT state_matrix, transition_matrix, C, T, number_zones, p_drive, p_dest, datamatrix, distance_matrix_km and OUTPUT state_matrix, transition_matrix

function resampling(state_matrix, transition_matrix, C, number_zones, p_drive, p_dest, datamatrix, distance_matrix_km)


        ## Sampling an entire day with imporved intitial value.
        for t=1:T
            

            # Sampling if a car drives
            for i=1:C
                origin = state_matrix[i,t];
                RndVar = rand();
                driving_probability = p_drive[origin,t];
                if(RndVar <= driving_probability)
                    drive = 1;
                else
                    drive = 0;
                    transition_matrix[i,t,2]= origin;
                end
                transition_matrix[i,t,1] = drive;
            end
            
            
            # Sampling travel destinations
            for i=1:C
                drive = transition_matrix[i,t,1];
                if (drive == 1)
                    RndVar = rand();
                    range_up = 0;
                    range_low = 0;
                    destination = 0;
                    origin = state_matrix[i,t];
                    distribution = p_dest[origin, :, t];
                    if(sum(distribution) == 0)
                        destination = origin;
                    else
                        for j=1:number_zones  
                            range_up = range_up + distribution[j];
                            if(range_low < RndVar <= range_up)
                                destination = j;
                                break;
                            end
                            range_low = range_up;
                        end
                    end
                    transition_matrix[i,t,2]= destination;
                end
            end
            
            
            # sampling travel times and distances 
            for i=1:C
                drive = transition_matrix[i,t,1];
                if(drive == 1)
                    origin = state_matrix[i,t,1];
                    destination = round.(Int, transition_matrix[i,t,2]);
                    if(origin == destination) # this is the case that a car drives to a destination within the origin zone
                        transition_matrix[i,t,3] = 5 * 60; # assuming arbitrarily a driving time of 5 minutes given in seconds
                        transition_matrix[i,t,4] = 1; # assuming arbitraliy a driving distance of 1 km given in km. 
                    else
                        # sampling travel times
                        mean = datamatrix[origin, destination, t, 1];
                        std_deviation = datamatrix[origin, destination, t, 2];
                        if(std_deviation == 0) # catches the case of missing std deviation data where it is set to 0 in the original datasets
                            std_deviation = 0.1 * mean;
                        end
                        travel_time = rand(Truncated(Normal(mean, std_deviation),0.9 * mean, 1.1* mean)); # travel time chosen from truncated normal distribution 
                        transition_matrix[i,t,3] = travel_time;

                        # sampling travel distances
                        mean = distance_matrix_km[origin, destination]; # constant zone distances as mean
                        std_deviation = 0.1 * mean; # arbitrary chosen standard deviation as 10% of the mean value, to create some variation.
                        travel_distance = rand(Truncated(Normal(mean, std_deviation),0.9 * mean, 1.1 * mean)); # travel distance chosen from truncated normal distribution.
                        transition_matrix[i,t,4] = travel_distance;
                    end
                end
            end

            # state update
            if(t<T)
                state_matrix[:,t+1] = transition_matrix[:,t,2];
            end
            println(trunc(Int,t/0.23),"% of the dataset is simulated")
        end

        return state_matrix, transition_matrix 

end