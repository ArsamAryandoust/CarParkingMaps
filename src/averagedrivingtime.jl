### Calculating average driving time: INPUT C, A_drive[N], transition_matrix and OUTPUT A_drive 

function averagedrivingtime(C, A_drive, transition_matrix)
    sum_driving_time = 0;
    driving_time_entries = zeros(C,1);
    for t=1:T
        driving_time_entries = transition_matrix[:,t,3];
        sum_driving_time = sum_driving_time + sum(driving_time_entries);
    end
    A_drive = A_drive  + sum_driving_time/(C * T * 60 * 60);
    return A_drive
end

