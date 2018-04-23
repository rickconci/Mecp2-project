function [table] = get_bursts_forRicardo_v230117MS(spike_train,plots)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: spike train structure
% Output: burst train, burst duration, number of bursts
% Dependencies: findPattern2.m
% Version: 19 March 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set some parameters
min_nr_spk = 50; % minimum number of spikes/channel
min_nr_ch = 6; % minimum number of channles
min_distance = 20; % = 200ms

%elec = size(spike_train,2);
%points = size(spike_train,1);

% Load spike times
SPK = zeros(size(spike_train,2) ,size(spike_train,1)); % only 9 mins, because I removed the 1st min
for ii = 1:size(spike_train,2)
    if sum(spike_train.locs{ii}) > 0 % spike_train.locs is a vector with all spike times of one channel
        if isnan(spike_train.locs{ii}) == 0
            spk_tms = spike_train.locs{ii};
            SPK(ii,spk_tms') = 1;
        end
    end
end
clear spike_train

% Resample spike trains to ms resolution
for ch = 1:size(spike_train, 2)
    spike_train(ch,:) = (sum(reshape(SPK(ch,1:size(spike_train,1)),25,[])))>0;
end

% Number of spikes in 10ms windows
mat = reshape(spike_train,size(spike_train,2),10,72000);
Burst_train_spks = (squeeze(sum(sum(mat,2)))); % number of spikes
Burst_train_active_channel = (sum(squeeze((sum(mat,2)))>0)); % number of channels
Burst_train_thresh = Burst_train_spks>0;
min_spks = zeros(1,length(Burst_train_thresh));
sum_spk_pos = find(Burst_train_thresh);

% Find spike cluster
for j = 1:length(sum_spk_pos)-1
    dif = sum_spk_pos(j+1)-sum_spk_pos(j); % get inter-spike cluster interval
    if dif < 6  % minimal inter spike cluster interval: 60ms
        min_spks(sum_spk_pos(j):sum_spk_pos(j+1))=1; % insert ones if distance smaller than thresh
    end
end
min_spks(1)=0;
min_spks(end)=0;

% Find bursts
onset = findPattern2(min_spks,[0 1])+1;
offset = findPattern2(min_spks,[1 0]);
burst_vec_1 = zeros(1,length(Burst_train_thresh));
burst_vec_final = zeros(1,length(Burst_train_thresh));
burst_1sec_vec = zeros(1,length(Burst_train_thresh));

for burst = 1:length(offset)
    nubo_spks(burst) = sum(Burst_train_spks(onset(burst):offset(burst)));
    max_ch(burst) = max(Burst_train_active_channel(onset(burst):offset(burst)));
    bb = onset(burst):offset(burst);
    burst_times(burst) = {bb};
    burst_vec_1(onset(burst):offset(burst)) = 1;
end

% Exclude events below thresholds (min spks / min ch)
nubo_spks_thrsh = find(nubo_spks>min_nr_spk);
max_ch_thrsh = find(max_ch>min_nr_ch);
burst_over_thresh  = intersect(nubo_spks_thrsh,max_ch_thrsh);
burst_times_thresh = burst_times(burst_over_thresh);
bursts_over_thresh_1 = length(burst_over_thresh);
zts = zeros(1,72000);

for burst = 1:bursts_over_thresh_1;%-1
    bt = cell2mat(burst_times_thresh(burst));
    zts(bt)=1; % get new burst vector
    clear bt
end

% Combine burst events closer 200ms
zts_da = find(zts);
bursts_over_thresh_2 = zeros(1,72000);

for j = 1:length(zts_da)-1
    dif = zts_da(j+1)-zts_da(j);  % get inter-event interval
    if dif < min_distance 
        bursts_over_thresh_2(zts_da(j):zts_da(j+1))=1; % insert ones if distance smaller than thresh
    end
end
onset2 = findPattern2(bursts_over_thresh_2,[0 1])+1;
offset2 = findPattern2(bursts_over_thresh_2,[1 0]);

% Burst onset/offset times
for burst = 1:length(onset2)
    bb2 = onset2(burst):offset2(burst);
    burst_times2(burst) = {bb2};
    burst_length_ms(burst,:) = length(bb2)*10;
end

% Summarize results in table
table.burst_train = Burst_train_spks;
table.number_of_burst = length(onset2);
table.bursts_over_thresh_1 = zts;
table.bursts_over_thresh_2 = bursts_over_thresh_2;
table.burst_windows = burst_times2;
table.burst_duration = burst_length_ms;

% Plot bursts (optional)
if plots == 1
    figure,plot(table.burst_train);
    hold on; plot((table.bursts_over_thresh_2).*table.burst_train','r--');
    [da,wo] = find(table.bursts_over_thresh_2); hold on; plot(wo,da,'k.');
end
legend('final burst','combined bursts')
end