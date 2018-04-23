function batch_ave_firing_rate

% Goes through the .mat files, create and save raster plots for each file
files = dir('*.mat');  % where your .mat files are 

for file = 1:length(files)
    data = load(files(file).name,'tSpikes'); %cNeg01254_Spikes
    binary_spikes = data.tSpikes; 
    name = files(file).name;
    
    %ask if any electrodes need to be removed from recording
    name
    pause
    prompt0 = 'any faulty electrodes? 0 or [electrode1 electrode2]';
    electrodes = input(prompt0);
    if electrodes ~= 0
        for e = 1:length(electrodes)
            electrodes = sort(electrodes, 'descend');
            binary_spikes(:,electrodes(e)) = 0;
        end
    end
 
    newspikes = binary_spikes;
    old_sampling_rate = size(newspikes, 1)/720; %readings per second
    samplesize = (size(newspikes,1)/25); %get resolution of 0.001
    new_sampling_rate = 1000;
    %find spike times
    spikeTimesCell = findSpikeTimes(newspikes, "seconds", old_sampling_rate);

    %turn cell into sorted vector of all spike times
    spikeTimesVec = cell2mat(spikeTimesCell');
    sortedSpikeTimesVec = sort(spikeTimesVec)';
    ave_firing_rate = length(sortedSpikeTimesVec)/(12*720);
    
     featMatrix(file).ID = files(file).name(4:14); 
     featMatrix(file).DIV = str2num(files(file).name(19:20)); 
     featMatrix(file).genotype = files(file).name(1:2);
     %spike features
     featMatrix(file).ave_firing_rate = ave_firing_rate;

     
     writetable(struct2table(featMatrix), fullfile(pwd, ...
          'ave_firing_rates_wavelet.xlsx'));
end
end

