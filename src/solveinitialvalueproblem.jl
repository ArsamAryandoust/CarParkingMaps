### Initial Value Problem: INPUT is state_matrix, p_drive, p_dest, C, T and OUTPUT initial_state


function solveinitialvalueproblem(state_matrix, transition_matrix, p_drive, p_dest, C, number_zones)

    initial_state = zeros(C,1);

    for t=1:(T-1)
        
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
                transition_matrix[i,t,2] = destination;
            end
        end

        
        # state update
        state_matrix[:,t+1] = transition_matrix[:,t,2];
        println(trunc(Int,t/0.23),"% of the Initial Value Problem is solved")
    end

    initial_state = state_matrix[:,24];
    initial_state = round.(Int, initial_state)

    return initial_state 

end