### Processing dataset: INPUT is path_to_csv_data and OUTPUT is datamatrix 

function createdatamatrix(path_to_csv_data, number_zones)
    rawdata = CSV.read(path_to_csv_data);
    rawdata = convert(Matrix, rawdata[:,1:5]);
    size_rawdata, placeholder_columns = size(rawdata);
    datamatrix= zeros(number_zones, number_zones, T, 2);
    for i=1:size_rawdata
        if(rawdata[i,1] == 0)
            rawdata[i,1] = number_zones; # in case a zone has an ID of "0", we give it a new zone ID. 
        end
        if(rawdata[i,2] == 0)
            rawdata[i,2] = number_zones; # in case a zone has an ID of "0", we give it a new zone ID. 
        end
        if(rawdata[i,3] == 0)
            rawdata[i,3] = 24; # Midnight is denoted by a value of zero. We replace it with 24 to get matrix indexes right.
        end
        index1 = convert(Int64, rawdata[i,1]);
        index2 = convert(Int64, rawdata[i,2]);
        index3 = convert(Int64, rawdata[i,3]);
        datamatrix[index1,index2,index3,1] = rawdata[i,4]; # mean_travel_time goes here
        datamatrix[index1,index2,index3,2] = rawdata[i,5]; # standard_deviation_travel_time goes here
    end

    return datamatrix
    
end