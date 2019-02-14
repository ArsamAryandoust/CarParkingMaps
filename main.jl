#Pkg.update()
#Pkg.add("CSV")
#Pkg.add("Distributions")
#Pkg.add("JSON")
#Pkg.add("FileIO")
#Pkg.add("DataFrames")
using CSV
using Distributions
using JSON
using FileIO
using DataFrames

### Provide the path to the source files
include("src/processgeodata.jl");
include("src/createdatamatrix.jl");
include("src/createpdrive.jl");
include("src/createpdestin.jl")
include("src/initializestates.jl")
include("src/solveinitialvalueproblem.jl")
include("src/resampling.jl")
include("src/averagedrivingtime.jl")
include("src/correctparameters.jl")
include("src/saveresults.jl")
include("src/saveparameters.jl")
include("src/createresultsdirectory.jl")

### Provide the path to the Uber datasets and the path to the folder where you want to save your results
path_to_cities = "/data/Uber Movement/";
path_to_results_folder = "/results/";
city_list = readdir(path_to_cities); #getting the names of directories

### Define these parameters manually
normalize = false; # a boolean expression that adapts the code accordingly, if we choose to normalize the traffic system
e = 2; # An exponent that determines the dynamics of driving activity. Can be changed.
cars_per_zone = 10; # the number of cars to be sampled per each city zone
T = 24; # The datasets we analyze have an hourly time resolution
A_set = 0.05; # the set point value for normalizing average driving time
max_iteration = 6; # maximum number of iteration for convergence of normalization algorithm


for city_index = 1:1#length(city_list)
    
    ### Set the path to city files that we want to simulate
    city = city_list[city_index];
    path_to_data = string(path_to_cities, city);
    dataset_list = readdir(path_to_data);
    path_to_json_data = string(path_to_data,"/",dataset_list[length(dataset_list)])
    csv_dataset_list = dataset_list[1:(length(dataset_list)-1)]; 
   
    ### create a directory for saving the results of the city
    path_to_results = createresultsdirectory(path_to_results_folder, city)

    ### process the GEOJSON data: INPUT is julia_file and OUTPUT is distance_matrix_km , number_zones
    distance_matrix_km, number_zones = processgeodata(path_to_json_data, path_to_data, csv_dataset_list, path_to_results);
    
    ### determine the number of cars that are simulated for the current city. 
    C = number_zones * cars_per_zone;
    C = trunc(Int, C); # transforms C from type Float64 to Integer. Necessary for subsequent matrix definitions. 
    println("Simulation starts for the city of ", city, " with a vehicle fleet of C = ", C);
    
    ### Initializing the values for Normalization Algorithm
    N = 1; # initializing count of iterations for normalization algorithm
    p_min = zeros(max_iteration,1); # saves the tuned p_min parameter in each iteration
    p_max = zeros(max_iteration,1); # saves the tuned p_max parameter in each iteration
    A_drive = zeros(max_iteration,1); # saves the computed travel time of sample in each iteration
    △A = zeros(max_iteration,1); # the deviation between A_set and average driving time of sample 
    p_min[1] = 0.1; # the initial parameters for unnormalized traffic sampling
    p_max[1] = 0.9; # the initial parameters for unnormalized traffic sampling
    

    while N  < max_iteration
    
        count_dataset = 0; # count the number of datasets 

        for dataset_index = 1:length(csv_dataset_list)
            
            ### Set the path to dataset that we want to simulate
            data_set = csv_dataset_list[dataset_index];
            path_to_csv_data = string(path_to_data,"/",data_set); 
            println("Simulating ", data_set)
            
            ### Processing dataset: INPUT is path_to_csv_data, number_zones and OUTPUT is datamatrix 
            datamatrix = createdatamatrix(path_to_csv_data, number_zones);

            ### Creating driving activity: INPUT is datamatrix, p_min, p_max, number_zones and OUTPUT is p_drive
            p_drive = createpdrive(datamatrix, p_min[N], p_max[N], number_zones);

            ### Creating destination popularity: INPUT is datamatrix,  number_zones and OUTPUT is p_dest
            p_dest = createpdestin(datamatrix, number_zones);

            ### Initializing state and state transitions: INPUT is C and OUTPUT state_matrix, transition_matrix
            state_matrix, transition_matrix = initializestates(C);

            ### Initial Value Problem: INPUT is state_matrix, p_drive, p_dest, C and OUTPUT initial_state
            initial_state = solveinitialvalueproblem(state_matrix, transition_matrix, p_drive, p_dest, C, number_zones)
            state_matrix[:,1] = initial_state; 

            ### Resampling: INPUT state_matrix, transition_matrix, C, number_zones, p_drive, p_dest, datamatrix, distance_matrix_km and OUTPUT state_matrix, transition_matrix
            state_matrix, transition_matrix = resampling(state_matrix, transition_matrix, C, number_zones, p_drive, p_dest, datamatrix, distance_matrix_km);

            ### Calculating average driving time: INPUT C, A_drive[N], transition_matrix and OUTPUT A_drive 
            A_drive[N] = averagedrivingtime(C, A_drive[N], transition_matrix)
            count_dataset = count_dataset + 1; # increment the number of datasets contained in this average 

            ### saving results
            saveresults(number_zones, state_matrix, path_to_results, data_set, C);

        end

        ### averages driving time by the number of datasets contained in A_drive and calculate deviation from set point value 
        A_drive[N] = A_drive[N]/ count_dataset; 
        △A[N] = A_set - A_drive[N];

        ### save parameters
        saveparameters(N, path_to_results, T, number_zones, cars_per_zone, C, e, p_min, p_max, A_set, A_drive, △A);


        ### Case conrols for normalization algorithm
        if (normalize==true)
            ### Case control: check if solution is improved and reset parameter values if this is not the case.
            if(trunc(Int, △A[N] * 100) == 0)
                println("A good solution is found for p_min = ",p_min[N]," and p_max = ",p_max[N], ". It leads to a △A = ", △A[N]);
                break
            else
                println("The solution is not good enough for p_min = ",p_min[N]," and p_max = ",p_max[N], ". It leads to a △A = ", △A[N])
            end

            if(N >= 2 && abs(△A[N]) > abs(A_set - A_drive[N]))
                println("The parameters p_min = ",p_min[N]," and p_max = ",p_max[N]," lead to |△A| = ",abs(△A[N]),". This is worse than the last solution so we reset the parameters to p_min = ",p_min[N-1]," and p_max = ",p_max[N-1]);
                p_min[N] = p_min[N-1];
                p_max[N] = p_max[N-1];
                △A[N]= A_set - A_drive[N-1];
                A_drive[N] = A_drive[N-1];
            end

            ### Case control: choosing which equation to solve, depending on whether parameters reached bounds or not.
            if(p_max[N]==1)
                if(p_min[N]==p_max[N])
                    println("something went wrong! p_min is equal to p_max!");
                    break
                else
                    p_min[N+1] = p_min[N] + (p_max[N] + p_min[N]) * △A[N] / (N * A_drive[N]);
                    p_max[N+1] = p_max[N];
                end
            else
                if(p_min[N]==0)
                    if(p_min[N]==p_max[N])
                        println("something went wrong! p_min is equal to p_max!");
                        break
                    else
                        p_max[N+1] = p_max[N] + (p_max[N] + p_min[N]) * △A[N] / (N * A_drive[N]);
                        p_min[N+1] = p_min[N];
                    end
                else
                    if(p_min[N]==p_max[N])
                        p_min[N+1] = p_min[N] + (p_max[N] + p_min[N]) * △A[N] / (N * A_drive[N]);
                        p_max[N+1] = p_max[N];
                    else
                        p_max[N+1] = p_max[N] + (p_max[N] + p_min[N]) * △A[N] / (N * A_drive[N]);
                        p_min[N+1] = p_min[N];
                    end
                end
            end

            ### Case control: correcting new parameters in case they cross any feasible bounds: INPUT p_min[N+1], p_max[N+1], p_min[N], p_max[N] OUTPUT p_min_next, p_max_next
            p_min[N+1], p_max[N+1]= correctparameters(p_min[N+1], p_max[N+1], p_min[N], p_max[N])

            ### increment iteration counter for normalization algorithm
            N = N + 1;

            ### indicate the case that the maximum number of iterations is reached
            if(N == max_iteration)
                println("Maximum number of iterations is reached. Best found solution is p_min = ",p_min[N]," and p_max = ",p_max[N])
            end

        else
            ### set N as max_iteration to exit the while loop
            N = max_iteration;
            println("Traffic has been sampled for ", city, " without normalization.")
        end

    end



end
