function batchHistIsi

% Goes through the .mat files, create and save raster plots for each file
files = dir('*.mat');  % where your .mat files are 

for file = 1:length(files)
data = load(files(file).name,'tSpikes'); 
spikes = data.tSpikes; 

%find spike times
spikeTimesCell = findSpikeTimes(spikes, "seconds", 25000);

%turn cell into sorted vector of all spike times
spikeTimesVec = cell2mat(spikeTimesCell');
sortedSpikeTimesVec = sort(spikeTimesVec)';

%input into histogram function
N = [2:15]; % Range of N for ISI_N histograms 
Steps = 10.^[-5:0.05:1.5]; % Create uniform steps for log plot 
HistogramISInTim(sortedSpikeTimesVec, N, Steps);

%save as png
printformat = "histogram ISI %s.png";
name = files(file).name;

set(gcf, 'Name', sprintf("%s", name));
saveas(gcf,sprintf(printformat, name));

close all

end
end

