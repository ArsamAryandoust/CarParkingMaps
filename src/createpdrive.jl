### Creating driving activity: INPUT is datamatrix, p_min, p_max, e, T, number_zones and OUTPUT is p_drive

function createpdrive(datamatrix, distance_matrix_km, number_zones)
    p_drive = zeros(number_zones, T);
    mean_sum = zeros(number_zones, T);
    min_t = zeros(number_zones,1);
    max_t = zeros(number_zones,1);

    # Calculating the sums, minimum sums and maximum sums of mean travel time out of a zone
    for i=1:number_zones    
        for t=1:T
            mean_sum[i, t] = 0;
            counter = 0;
            for j=1:number_zones
                if(datamatrix[i,j,t,1] != 0)
                    mean_sum[i, t] = mean_sum[i, t] + datamatrix[i,j,t,1]/distance_matrix_km[i,j];
                    counter += 1;
                end
            end
            mean_sum[i, t] = mean_sum[i, t]/counter;
        end
        max_t[i] = maximum(mean_sum[i,:]);
        min_t[i] = minimum(mean_sum[i,:]);
    end

    ## Calculating the Bernoulli distributions p_drive
    for i=1:number_zones
        if(max_t[i] > 0) # A zero would lead to infinity given the formula below and create NaN entries.
            for t=1:T
                # A linear function for assigning a value of 10% to p_drive if min_t is given, and a value of 90% if max_t is given.
                p_drive[i,t] = p_min + (p_max - p_min) * ((mean_sum[i,t] - min_t[i]) / (max_t[i] - min_t[i]))^e_drive;  
            end
        end
    end
    
    return p_drive;

end