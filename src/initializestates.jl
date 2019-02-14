### Initializing state and state transitions: INPUT is C, cars_per_zone, T and OUTPUT state_matrix, transition_matrix


function initializestates(C)

    state_matrix = zeros(C, T); # C X T  matrix for saving states of the system, i.e. the location of cars in each timestep. The rows stand for car ID and the columns for the timestep.
    transition_matrix = zeros(C, T, 4); # C X T X 4 matrix for saving state transitions of the system in each timestep. The entries for each car ID (C) and each time step (T) are the driving activity, destination, travel time and travel distance.
    zone = 0;
        
    ## assigning a zone uniformely to each car
    for i = 0:cars_per_zone:(C-cars_per_zone)
        zone = zone + 1;
        for j=1:cars_per_zone
            state_matrix[i+j, 1] = zone; # each car receives a location zone in time step 1, i.e. at 1h o'clock 
        end  
    end
    state_matrix = round.(Int, state_matrix);

    return state_matrix, transition_matrix


end