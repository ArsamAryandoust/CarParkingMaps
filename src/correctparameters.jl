### Case control: correcting new parameters in case they cross any feasible bounds: INPUT p_min[N+1], p_max[N+1], p_min[N], p_max[N] OUTPUT p_min_next, p_max_next

function correctparameters(p_min_next, p_max_next, p_min, p_max)

    if(p_min_next < 0)
        p_min_next = 0;
    else
        if(p_min_next > p_max)
            p_min_next = p_max;
        end
    end

    if(p_max_next > 1)
        p_max_next = 1;
    else
        if(p_max_next < p_min)
            p_max_next = p_min;
        end
    end    

    return p_min_next, p_max_next

end