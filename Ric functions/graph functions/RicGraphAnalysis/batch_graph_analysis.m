function batch_graph_analysis

% Goes through the .mat files, create and save raster plots for each file
files = dir('*.mat');  % where your .mat files are 

for file = 1:length(files)
    data = load(files(file).name,'tSpikes'); 
    binary_spikes = data.tSpikes; 
    name = files(file).name;
    
    old_sampling_rate = size(binary_spikes, 1)/720; %readings per second
    samplesize = (size(binary_spikes,1)/25); %get resolution of 0.001
    new_sampling_rate = 1000;
    %find spike times
    spikeTimesCell = findSpikeTimes(binary_spikes, "seconds", old_sampling_rate);

    %turn cell into sorted vector of all spike times
    spikeTimesVec = cell2mat(spikeTimesCell');
    sortedSpikeTimesVec = sort(spikeTimesVec)';

    %input into histogram function
    N = [2:15]; % Range of N for ISI_N histograms 
    Steps = 10.^[-5:0.05:1.5]; % Create uniform steps for log plot 
    HistogramISInTim(sortedSpikeTimesVec, N, Steps);
    title(name)
    
    %input N and ISI_N after seeing the histogram
    prompt1 = 'What is the value of N? ';
    new_N = input(prompt1);
    prompt2 = 'What is the value of ISI_N? ';
    new_ISI_N = input(prompt2);
    
    Spike.T = sortedSpikeTimesVec;
    Spike.C = size(binary_spikes, 2);
    %Initiate burst detection algorithm
    [Burst,SpikeBurstNumber] = BurstDetectISInTim( Spike, new_N, new_ISI_N ); 

    number_of_bursts = max(SpikeBurstNumber)
    
    pause
    
    %relate the burst times back to binary spike matrix
    for ii = 1:max(SpikeBurstNumber)
      spikes_down = sparseDownSample(binary_spikes, samplesize, "sum");
      
      %spikes_down is a matrix with 1000 binary per second.. 
      spikes_down_full = full(spikes_down)';
      burst_matrix = spikes_down_full(:,round(Burst.T_start(ii)*new_sampling_rate, 3):round(Burst.T_end(ii)*new_sampling_rate,3));
      [row, col] = find(burst_matrix);
      for j = 1:length(col)-1
        dif = col(j+1)-col(j);  % get inter-event interval 
        if dif < new_ISI_N*1000
            burst_matrix(row, col:col+new_ISI_N*1000) = 1; %add 1s to increase signal of network burst
        end
      end 
      %rastPlot(burst_matrix');
      %find and save association matrix for each burst 
      ass_matrix = corrcoef(burst_matrix');
      mean_of_burst = mean(ass_matrix(:));
      ass_matrix(isnan(ass_matrix) == 1) = mean_of_burst;
      matFileName = sprintf('%s association matrix burst %d.mat', name, ii);
      save(matFileName, "ass_matrix")
      
      %find the real adjacency matrix of each association matrix
      burstAdjMatrix = get_real_adj_matrix(ass_matrix);
      %matFileName2 = sprintf('%s adjacency matrix burst %d_.mat', name, ii);
      %save(matFileName2, "burstAdjMatrix")
      
%Burst features 
      %burt duration 
      burst_duration = (Burst.T_start(ii)*new_sampling_rate - Burst.T_end(ii)*new_sampling_rate);
      
      %inburst firing rate (spikes/ms)
      inburst_firing_rate = Burst.S(ii)/burst_duration;
      
%Graph theory features
      %density
      %kden = density
      %N = number of vertices
      %K = number of edges
      [kden,N,K] = density_und(burstAdjMatrix);

      %average, max degree
      [deg] = degrees_und(burstAdjMatrix);
      ave_deg = mean(deg);
      max_deg = max(deg);
       
      %clustering coefficient is a vector
      clustering_coefficient =clustering_coef_bu(burstAdjMatrix);
      clustered_nodes = length(find(clustering_coefficient));

      %average path length
      %G = graph(burstAdjMatrix);
      %nodedistances = distances(G);
      %ave_path_length = mean(nodedistances);
      
      %Rich Club Topology
      %R = vector of rich-club coefficients for levels 1 to klevel.
      %Nk = number of nodes with degree>k
      %Ek = number of edges remaining in subgraph with degree>k
      [R,Nk,Ek] = rich_club_bu(burstAdjMatrix);
    
      %number of hubs
      hub_nodes = findHubs(burstAdjMatrix);
      num_hubs = length(hub_nodes);
      
      
      featMatrix(ii).batch = files(file).name(4:11); 
      featMatrix(ii).ID = files(file).name(4:14); 
      featMatrix(ii).DIV = str2num(files(file).name(19:20)); 
      featMatrix(ii).genotype = files(file).name(1:2);
      %burst features
      featMatrix(ii).burst_number = ii;
      featMatrix(ii).burst_duration = burst_duration;
      featMatrix(ii).inburst_firing_rate = inburst_firing_rate;
      %graph theory features 
      featMatrix(ii).density = kden;
      featMatrix(ii).ave_deg = ave_deg;
      featMatrix(ii).max_deg = max_deg;
      %featMatrix(ii).ave_path_length = ave_path_length;
      %featMatrix(ii).clustering_coefficient = clustering_coefficient;
      featMatrix(ii).clustered_nodes = clustered_nodes;
      featMatrix(ii).richclub = R;
      featMatrix(ii).num_hubs = num_hubs;
     
      
      %save all plots as in folder 
      figHandles = findobj('Type', 'figure');
      for k = 1:length(figHandles)
         printformat = "%s figure %d of burst %d.png";
         str = sprintf(printformat, name, k, ii);
         saveas(figHandles(k),str);
      end
      close all
    save('featMatrixM5', 'featMatrix'); 
    end 
     
end 

