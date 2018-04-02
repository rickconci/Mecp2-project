function batch_raster_plot

% Goes through the .mat files, create and save raster plots for each file

files = dir('*.mat');  % where your .mat files are 

for file = 1:length(files)
data = load(files(file).name,'tSpikes'); 
spikes = data.tSpikes; 

samplesize = size(spikes,1)/1000;
%downsample so takes shorter and can actually see spikes in raster plot
spikes_down = sparseDownSample(spikes, samplesize, 'sum');
%create raster plot
rastPlot(full(spikes_down));

%save raster plot as PNG
printformat = "rasterPlot %s.png";
name = files(file).name;

set(gcf, 'Name', sprintf("%s", name));
saveas(gcf,sprintf(printformat, name));

close all
end
    
end

