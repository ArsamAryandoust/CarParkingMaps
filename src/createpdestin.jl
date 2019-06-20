### Creating destination popularity: INPUT is datamatrix, e, T, number_zones and OUTPUT is p_dest

function createpdestin(datamatrix, number_zones)
    p_dest = zeros(number_zones, number_zones, T);
    max_x = zeros(number_zones, number_zones);
    min_x = zeros(number_zones, number_zones);
    normalization_factor = zeros(number_zones, T);
        
    ## Calculating the minima and maxima of mean travel time 
    for i=1:number_zones 
        # Calculating the minima and maxima
        for j=1:number_zones # destination zone ID
            max_x[i,j] = maximum(datamatrix[i, j, :, 1])
            min_x[i,j] = minimum(datamatrix[i, j, :, 1])
        end
    end
        
    ## Calculating the unnormalized p_dest
    for i=1:number_zones 
        for j=1:number_zones
            for t=1:T
                if (max_x[i, j] > 0)
                    mean = datamatrix[i, j, t, 1];
                    p_dest[i, j, t] = ((mean - min_x[i, j])/ (max_x[i, j] - min_x[i, j]))^e_dest;
                end
            end
        end
    end
        
    ## Calculating normalization factors
    for i=1:number_zones
        for t=1:T
            normalization_factor[i, t] = sum(p_dest[i, :, t]);
        end
    end
        
    ## Normalizing the probability distributions of p_dest
    for i=1:number_zones 
        for t=1:T
            if(normalization_factor[i, t] > 0) # A zero would lead to infinity given the formula below and create NaN entries.
                for j=1:number_zones
                    p_dest[i, j, t] = p_dest[i, j, t] / normalization_factor[i, t];
                end
            end
        end
    end

    return p_dest

end