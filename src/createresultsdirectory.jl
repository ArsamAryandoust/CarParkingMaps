### create a directory for saving the results of the city

function createresultsdirectory(path_to_results_folder, city)
    
    path_to_results = string(path_to_results_folder,city);

    folder_list = readdir(path_to_results_folder);
    for i = 1:length(folder_list)
        folder_name = folder_list[i];
        if (folder_name == city) # catches the case that a results folder already exists for this city
            return path_to_results
        end
    end

    mkdir(path_to_results);
    return path_to_results

end