function burst_trains = get_burst_trains(burst_matrix)
%take the burst_matrix (1's if burst, and 0 if no burst)

sampling_rate = 720*10;
burst_order = find (burst_matrix == 1);
spike_times = burst_order/sampling_rate;


    
tot_bursts = sum(sum(burst_matrix==1));

for number in burst_times
    data_points = number*

end 