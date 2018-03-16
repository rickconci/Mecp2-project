function batchHistIsi

% Goes through the .mat files, create and save raster plots for each file
files = dir('*.mat');  % where your .mat files are 

for file = 1:length(files)
data = load(files(file).name,'tSpikes'); 
spikes = data.tSpikes; 

spikeTimes = findSpikeTimes(spikes, "seconds");



end

