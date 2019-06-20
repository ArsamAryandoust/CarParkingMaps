#Pkg.update()
#Pkg.add("CSV")
#Pkg.add("Distributions")
#Pkg.add("JSON")
#Pkg.add("FileIO")
#Pkg.add("DataFrames")
using CSV, Distributions, JSON, FileIO, DataFrames


### Import source files
include("src/processgeodata.jl");
include("src/createdatamatrix.jl");
include("src/createpdrive.jl");
include("src/createpdestin.jl");
include("src/initializestates.jl");
include("src/solveinitialvalueproblem.jl");
include("src/resampling.jl");
include("src/averagedrivingtime.jl");
include("src/correctparameters.jl");
include("src/saveresults.jl");
include("src/saveparameters.jl");
include("src/createresultsdirectory.jl");


### Provide the path to the Uber datasets and the path to the folder where you want to save your results
path_to_cities = "/Users/undisputed/Documents/data/Uber Movement/";
path_to_results_folder = "/Users/undisputed/Documents/data/Results/";
city_list = readdir(path_to_cities); #getting the names of directories


############### Insert a list of particular cities that you want to sample here ##############
#city_list = [" ***Insert a city here*** "]
city_list = ["Sydney"]
##############################################################################################


############################### Define these parameters manually #############################
e_drive = 0.5; # An exponent that determines the dynamics of driving activity. Can be changed.
e_dest = 2; # An exponent that determines the dynamics of destination popularity. Can be changed.
p_min = 0.1; # A parameter for p_drive. Can be changed.
p_max = 0.9; # A parameter for p_drive. Can be changed.
cars_per_zone = 1000; # the number of cars to be sampled per each city zone
T = 24; # The datasets we analyze have an hourly time resolution
##############################################################################################


for city_index = 1:length(city_list)
    
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
    
    ### initialize values for counting
    A_drive = 0;
    count_dataset = 0;

    ### Starts the simulation for the data sets of the current city
    for dataset_index = 1:length(csv_dataset_list)
            
        ### Set the path to dataset that we want to simulate
        data_set = csv_dataset_list[dataset_index];
        path_to_csv_data = string(path_to_data,"/",data_set); 
        println("Simulating ", data_set)
            
        ### Processing dataset: INPUT is path_to_csv_data, number_zones and OUTPUT is datamatrix 
        datamatrix = createdatamatrix(path_to_csv_data, number_zones);

        ### Creating driving activity: INPUT is datamatrix, p_min, p_max, number_zones and OUTPUT is p_drive
        p_drive = createpdrive(datamatrix, distance_matrix_km, number_zones);

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
        A_drive = averagedrivingtime(C, A_drive, transition_matrix)
        count_dataset += 1; # increment the number of datasets contained in this average 

        ### saving results
        saveresults(number_zones, state_matrix, transition_matrix, path_to_results, data_set, C);

    end

    ### Average A_drive for all sampled datasets of the current city
    A_drive  = A_drive/count_dataset;

    ### save parameters
    saveparameters(path_to_results, T, number_zones, cars_per_zone, C, e_drive, p_min, p_max, e_dest, A_drive);



end

