function batch_num_bursts

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
    
    %input into histogram function
    N = 2:15; % Range of N for ISI_N histograms 
    Steps = 10.^(-5:0.05:1.5); % Create uniform steps for log plot 
    HistogramISInTim(sortedSpikeTimesVec, N, Steps);
    title(name)
    
    %input N and ISI_N after seeing the histogram
    prompt1 = 'What is the value of N? N or 0 if no bursts';
    new_N = input(prompt1);
    
    if new_N == 0
        number_of_bursts = 0;
    elseif new_N ~= 0
   
    prompt2 = 'What is the value of ISI_N? ';
    new_ISI_N = input(prompt2);
    
    Spike.T = sortedSpikeTimesVec;
    Spike.C = size(newspikes, 2);
    %Initiate burst detection algorithm
    [Burst,SpikeBurstNumber] = BurstDetectISInTim( Spike, new_N, new_ISI_N ); 

    number_of_bursts = max(SpikeBurstNumber);
    end
     featMatrix(file).ID = files(file).name(4:14); 
     featMatrix(file).DIV = str2num(files(file).name(19:20)); 
     featMatrix(file).genotype = files(file).name(1:2);
     %spike features
     featMatrix(file).number_of_bursts = number_of_bursts;

     
     writetable(struct2table(featMatrix), fullfile(pwd, ...
          'number_of_bursts_NEO.xlsx'));

end
end
