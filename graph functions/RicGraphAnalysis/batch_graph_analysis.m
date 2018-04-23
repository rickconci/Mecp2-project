function batch_graph_analysis

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
    if size(binary_spikes, 2) == 58
        newspikes = zeros(size(binary_spikes,1),60);
        newspikes(:,1:3) = binary_spikes(:,1:3);
        newspikes(:,4) = zeros(size(binary_spikes,1), 1);
        newspikes(:,5:13) = binary_spikes(:,4:12);
        newspikes(:,14) = zeros(size(binary_spikes,1), 1);
        newspikes(:,15:60)= binary_spikes(:,13:58);
    elseif size(binary_spikes, 2) == 59
        newspikes = zeros(size(binary_spikes,1),60);
        newspikes(:,1:3) = binary_spikes(:,1:3);
        newspikes(:,4) = zeros(size(binary_spikes,1), 1);
        newspikes(:,5:60) = binary_spikes(:,4:59);
    elseif size(binary_spikes, 2) == 60
        newspikes = binary_spikes;
    end
    
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
    prompt1 = 'What is the value of N? ';
    new_N = input(prompt1);
    prompt2 = 'What is the value of ISI_N? ';
    new_ISI_N = input(prompt2);
    
    Spike.T = sortedSpikeTimesVec;
    Spike.C = size(newspikes, 2);
    %Initiate burst detection algorithm
    [Burst,SpikeBurstNumber] = BurstDetectISInTim( Spike, new_N, new_ISI_N ); 

    number_of_bursts = max(SpikeBurstNumber)
    
    pause
    
    %relate the burst times back to binary spike matrix
    for ii = 1:max(SpikeBurstNumber)
      spikes_down = sparseDownSample(newspikes, samplesize, "sum");
      
      %spikes_down is a matrix with 1000 binary per second.. 
      spikes_down_full = full(spikes_down)';
      burst_matrix = spikes_down_full(:,round(Burst.T_start(ii)*...
          new_sampling_rate, 3):round(Burst.T_end(ii)*new_sampling_rate,3));
      %show raster plot of spikes in burst
      %figure ('visible','off')
      %rastPlot(burst_matrix');
      
      if size(burst_matrix,1) == size(newspikes,2)
      %word = num2words(ii);
      concatenated_burst_struct(ii).burst_spikes = burst_matrix;
      
      %amplify signal during burst
      [row, col] = find(burst_matrix);
      for jk = 1:length(col)
         if size(burst_matrix,2)-col(end)>5
             %add 1s to increase signal of network burst
            burst_matrix(row(jk), col(jk):col(jk)+5) = 1; 
         elseif size(burst_matrix,2)-col(end)<5
             burst_matrix(row(jk), col(jk):col(jk)+5) = 1;
             burst_matrix(row(end), col(end):end-1) = 1;
         end
      end
      
      %show raster plot of spikes once amplified
      rasterburstprintformat = 'Raster Plot Burst %d';
      RasterBurstTitle = sprintf(rasterburstprintformat, ii);
      figure ('visible','off')
      rastPlot(burst_matrix');
      xlabel('Time (ms)');
      title(RasterBurstTitle);
      
      %find and save association matrix for each burst 
      ass_matrix = corrcoef(burst_matrix');
      %matFileName = sprintf('%s association matrix burst %d.mat', name, ii);
      %save(matFileName, "ass_matrix")
      
      %compare thresholded and binirised matrices
      ass_matrix(isnan(ass_matrix) == 1) = 0; 
      
      binary_matrix = binirise(ass_matrix);
      binary_matrix(logical(eye(size(binary_matrix)))) = 1;
      
      %burstAdjMatrix = get_real_adj_matrix(ass_matrix);
      %burstAdjMatrix(logical(eye(size(burstAdjMatrix)))) = 1;
      
      %plot heatmap of binary matrix
      figure ('visible','off')
      binaryprintformat = 'Binary Adj Matrix Burst %d';
      BinaryTitle = sprintf(binaryprintformat, ii);
      binary_matrix_heatmap = heatmap(binary_matrix,'Colormap', jet, ...
          'Title', BinaryTitle);
      
      
%Burst features 
      %burt duration 
      burst_duration = abs((Burst.T_start(ii)*new_sampling_rate -...
          Burst.T_end(ii)*new_sampling_rate));
      
      %inburst firing rate (spikes/ms)
      inburst_firing_rate = abs(Burst.S(ii)/burst_duration);
      
%graph features
       %kden = density %N = number of vertices %K = number of edges
      %[kdenburst,Nburst,Kburst] = density_und(burstAdjMatrix);
      [kdenbinary,Nbinary,Kbinary] = density_und(binary_matrix);
      
     
      %plot to MEA
      %threholdedprintformat2 = 'Thresholded Burst Network %d';
      %ThresholdedTitle2 = sprintf(threholdedprintformat2, ii);
      %GraphBurstAdj = plottoMEA(burstAdjMatrix, ThresholdedTitle2);
      
      binaryprintformat2 = 'Binary Burst Network %d';
      BinaryTitle2 = sprintf(binaryprintformat2, ii);
      GraphBinaryAdj = plottoMEA(binary_matrix, BinaryTitle2);
      
      
      
      
%similarity features
      
    %node similarity
      [sOut,tOut] = findedge(GraphBinaryAdj);
      %nodes in network burst = nnb
      nnb = unique(vertcat(unique(sOut), unique(tOut))); 
      num_nnb = length(nnb);
      nodes_in_network(ii).nodes = unique(vertcat(unique(sOut), unique(tOut)));        
      
      %word = num2words(ii);
      if length(binary_matrix) == size(newspikes,2)
        burst_matrix_struct(ii).burst_matrix = binary_matrix;
      end
      
      %density of connected nodes
      density = (Kbinary - 60)/((num_nnb*(num_nnb-1))/2);
      
      featMatrixBurst(ii).ID = files(file).name(4:14); 
      featMatrixBurst(ii).DIV = str2num(files(file).name(19:20)); 
      featMatrixBurst(ii).genotype = files(file).name(1:2);
      %spike features
      featMatrixBurst(ii).ave_firing_rate = ave_firing_rate;
      %burst features
      featMatrixBurst(ii).burst_number = ii;
      featMatrixBurst(ii).burst_N = new_N;
      featMatrixBurst(ii).burst_ISIn = new_ISI_N;
      featMatrixBurst(ii).burst_duration = burst_duration;
      featMatrixBurst(ii).inburst_firing_rate = inburst_firing_rate;
      %featMatrixBurst(ii).elec_in_burst = ;
      featMatrixBurst(ii).num_nodes_unthreholded = length(nodes_in_network(ii).nodes);
      featMatrixBurst(ii).num_edges_unthresholded = Kbinary - 60;
      featMatrixBurst(ii).density = density;
       
      if  density<1 && num_nnb > 6
          
       [realvsRN, eigenvector_mean_ratio, eigenvector_max_ratio,...
        num_rich_club_nodes, num_hubs ,num_hubs_RN_ratio, assortativity_ratio,...
        modularity_ratio] = nullstats(binary_matrix);
      
      featMatrixBurst2(ii).ID = files(file).name(4:14); 
      featMatrixBurst2(ii).DIV = str2num(files(file).name(19:20)); 
      featMatrixBurst2(ii).genotype = files(file).name(1:2);
      %spike features
      featMatrixBurst2(ii).ave_firing_rate = ave_firing_rate;
      %burst features
      featMatrixBurst2(ii).burst_number = ii;
      featMatrixBurst2(ii).burst_N = new_N;
      featMatrixBurst2(ii).burst_ISIn = new_ISI_N;
      featMatrixBurst2(ii).burst_duration = burst_duration;
      featMatrixBurst2(ii).inburst_firing_rate = inburst_firing_rate;
      %featMatrixBurs2t(ii).elec_in_burst = ;
      featMatrixBurst2(ii).num_nodes_unthreholded = length(nodes_in_network(ii).nodes);
      featMatrixBurst2(ii).num_edges_unthresholded = Kbinary - 60;
      featMatrixBurst2(ii).density = density;
      featMatrixBurst2(ii).realvsRN = realvsRN;
      featMatrixBurst2(ii).eigenvector_mean_ratio = eigenvector_mean_ratio;
      featMatrixBurst2(ii).eigenvector_max_ratio = eigenvector_max_ratio;
      featMatrixBurst2(ii).num_rich_club_nodes = num_rich_club_nodes;
      featMatrixBurst2(ii).num_hubs = num_hubs;
      featMatrixBurst2(ii).num_hubs_RN_ratio = num_hubs_RN_ratio;
      featMatrixBurst2(ii).assortativity_ratio = assortativity_ratio;
      featMatrixBurst2(ii).modularity_ratio = modularity_ratio;
      
      writetable(struct2table(featMatrixBurst2), fullfile(pwd, ...
          'featureMatrixBurst2.xlsx'));
      end 
      
      figHandles = findobj('Type', 'figure');
      for k = 1:length(figHandles)
         printformat = "%s figure %d of burst %d.png";
         str = sprintf(printformat, name, k, ii);
         saveas(figHandles(k),str);
      end
      close all
      end
    end
   

    %save per burst feat_matrix 
    %printformat3 = 'featureMatrixBurstN%dS%g.xlsx';
    %str3 = sprintf(printformat3, new_N, new_ISI_N);
    save('featureMatrixBurst', 'featMatrixBurst'); 
    writetable(struct2table(featMatrixBurst), fullfile(pwd, ...
        'featureMatrixBurst.xlsx'));
%% similarity per recording
    
     %similarity features   
     %node similarity
     added_nodes = cat(1,nodes_in_network.nodes);
     %tbl = tabulate(added_nodes);
     bins = 1:length(binary_matrix)+1;
     [N,edges] = histcounts(added_nodes, bins);
     node_similarity = sum(N(N>1))/numel(added_nodes);
     
     %edge similarity 
     burst_matrix_cell = struct2cell(burst_matrix_struct);
     burst_matrix_array = reshape(cell2mat(burst_matrix_cell),...
         length(binary_matrix),[],size(burst_matrix_cell,3));
     ave_burst_matrix = sum(burst_matrix_array, 3)/size(burst_matrix_array, 3);
     ave_burst_matrix(logical(eye(size(ave_burst_matrix)))) = 0; 
     %remove 1s from diagnal
     
     edge_similarity = mean(nonzeros(ave_burst_matrix(:))); 
     %mean of non zero elemnts
    
     featMatrixsimilarity(file).ID = files(file).name(4:14); 
     featMatrixsimilarity(file).DIV = str2num(files(file).name(19:20)); 
     featMatrixsimilarity(file).genotype = files(file).name(1:2);
     %spike features
     featMatrixsimilarity(file).ave_firing_rate = ave_firing_rate;
     featMatrixsimilarity(file).node_similarity = node_similarity;
     featMatrixsimilarity(file).edge_similarity = edge_similarity;
 
      
     save('featureMatrixSimilarity', 'featMatrixsimilarity'); 
     writetable(struct2table(featMatrixsimilarity), fullfile(pwd,...
           'featureMatrixSimilarity.xlsx'));
     
    
      
%% concatenated bursts metrics per recording
     %concatenated features
     for cc = 1:ii
         concatenated_burst_struct(cc).burst_spikes =...
             concatenated_burst_struct(cc).burst_spikes';
     end 
     concatenated_burst = cat(1,concatenated_burst_struct.burst_spikes);
     %amplify signal during burst
     concatenated_burst = concatenated_burst';
     [row, col] = find(concatenated_burst);
     for j = 1:length(col)
         if size(concatenated_burst,2)-col(end)>5
            concatenated_burst(row(j), col(j):col(j)+5) = 1;
            %add 1s to increase signal of network burst
         elseif size(concatenated_burst,2)-col(end)<5
             concatenated_burst(row(end), col(end):end-1) = 1;
         end
     end
     concatenated_burst(:,1) = 0;
     %show raster plot of spikes once amplified
     rasterconcprintformat = 'Raster Plot Concatenated Bursts DIV %d';
     RasterConcTitle = sprintf(rasterconcprintformat, file);
     figure ('visible','off')
     rastPlot(concatenated_burst');
     xlabel('Time (ms)');
     title(RasterConcTitle);
     %[burst_lengths] = Burst.T_end
     %yy = 1:size(concatenated_burst,1);
     %line(
     
     
     conc_ass_matrix = corrcoef(concatenated_burst');
     conc_ass_matrix(isnan(conc_ass_matrix) == 1) = 0;
     conc_ass_matrix = (conc_ass_matrix + conc_ass_matrix.')/2;
     
     conc_thr_matrix = weight_thresholding(conc_ass_matrix);
     %plot heatmap and mea network of concatenated and average burst
     %matrices
     figure ('visible','off')
     conc_ass_heatmap = heatmap(conc_thr_matrix,'Colormap', jet); 
     title('Thresholded Concatenated Burst Network');
     G_conc = plottoMEA(conc_thr_matrix, 'Thresholded Concatenated Burst Network'); 
     [sOut1,tOut1] = findedge(G_conc);
     nnwc = unique(vertcat(unique(sOut1), unique(tOut1)));
     num_nnwc = length(nnwc);
     
     
     thr_ave_burst_matrix = weight_thresholding(ave_burst_matrix);
     figure ('visible','off')
     ave_burst_ass_heatmap = heatmap(thr_ave_burst_matrix,'Colormap', jet); 
     title('Average Burst Network');
     G_ave_matrix = plottoMEA(thr_ave_burst_matrix, 'Average Burst Network');
     [sOut2,tOut2] = findedge(G_ave_matrix);
     nnab = unique(vertcat(unique(sOut2), unique(tOut2)));
     num_nnab = length(nnab);
     
     %burst metrics
     concatenated_time = size(concatenated_burst,2); %in ms
     tot_bursts_DIV = ii;
     
     if isempty(nnwc) ==1 && isempty(nnab) ==1
      break;
     elseif length(nnwc) > 6 && length(nnab) <6
      
     
      %density: kden = density %N = number of vertices %K = number of edges
      [kdenconc,Nconc,Kconc] = density_und(conc_thr_matrix);
       
      
      density_conc = (Kconc)/((num_nnwc*(num_nnwc-1))/2);
     
         
     [realvsRN, eigenvector_mean_ratio, eigenvector_max_ratio,...
         num_rich_club_nodes, num_hubs ,num_hubs_RN_ratio, assortativity_ratio,...
         modularity_ratio] = nullstats(conc_thr_matrix);
     
        %features for conc burst matrix vs random
      featMatrixDIV(file).ID = files(file).name(4:14); 
      featMatrixDIV(file).DIV = str2num(files(file).name(19:20)); 
      featMatrixDIV(file).genotype = files(file).name(1:2);
     %spike features
      featMatrixDIV(file).ave_firing_rate = ave_firing_rate;
      featMatrixDIV(file).node_similarity = node_similarity;
      featMatrixDIV(file).edge_similarity = edge_similarity ;
      %featMatrixDIV(file).mutual_info = ;
      %featMatrixDIV(file).n_repeating_patterns = ;
     %concatenated features
      featMatrixDIV(file).concatenated_time = concatenated_time;
      featMatrixDIV(file).tot_bursts_DIV = tot_bursts_DIV;
     %graph theory features of concatenated bursts
      featMatrixDIV(file).density_conc = density_conc;
      featMatrixDIV(file).active_nodes_conc = length(nnwc);
      %featMatrixDIV(file).ave_weight = ave_weight_conc;
      featMatrixDIV(file).degree_realvsRN_conc = realvsRN;
      %featMatrixDIV(file).deg_distribution_test = hrealvsRN;
      %featMatrixDIV(file).pvaluerealvsRN = pdistributions;
      %featMatrixDIV(file).size_max_connected_ratio = ;
      %featMatrixDIV(file).num_connected_ratio = ;
      %featMatrixDIV(file).mean_deg_ratio = ;
      %featMatrixDIV(file).sd_deg_ratio = ;
      %featMatrixDIV(file).max_deg_ratio = ;
      %featMatrixDIV(file).local_edge_density_ratio = ;
      featMatrixDIV(file).eigenvector_mean_ratio_conc = eigenvector_mean_ratio;
      featMatrixDIV(file).eigenvector_max_ratio_conc = eigenvector_max_ratio;
      featMatrixDIV(file).num_rich_club_nodes_conc = num_rich_club_nodes;
      featMatrixDIV(file).num_hubs_conc = num_hubs;
      featMatrixDIV(file).num_hubs_RN_ratio_conc = num_hubs_RN_ratio;
      featMatrixDIV(file).assortativity_ratio_conc = assortativity_ratio;
      featMatrixDIV(file).modularity_ratio_conc = modularity_ratio;
      
       save('featureMatrixConcatenated', 'featMatrixDIV'); 
       if isnan(featMatrixDIV(file).assortativity_ratio_ave_burst) == 0
       writetable(struct2table(featMatrixDIV), fullfile(pwd,...
           'featureMatrixConcatenated.xlsx'));
       end
      elseif length(nnwc) < 6 && length(nnab) > 6
         
   
      %ave deg ave matrix
     [kdenaveburst,Naveburst,Kaveburst] = density_und(thr_ave_burst_matrix); 
     
     density_ave = (Kaveburst)/((num_nnab*(num_nnab-1))/2);
     
     [realvsRN2, eigenvector_mean_ratio2, eigenvector_max_ratio2,...
         num_rich_club_nodes2, num_hubs2 ,num_hubs_RN_ratio2, assortativity_ratio2,...
         modularity_ratio2] = nullstats(thr_ave_burst_matrix);
     
      %features for ave burst matrix vs random
      featMatrixaveburstDIV(file).ID = files(file).name(4:14); 
      featMatrixaveburstDIV(file).DIV = str2num(files(file).name(19:20)); 
      featMatrixaveburstDIV(file).genotype = files(file).name(1:2);
     %spike features
      featMatrixaveburstDIV(file).ave_firing_rate = ave_firing_rate;
      featMatrixaveburstDIV(file).node_similarity = node_similarity;
      featMatrixaveburstDIV(file).edge_similarity = edge_similarity ;
      %featMatrixDIV(file).mutual_info = ;
      %featMatrixDIV(file).n_repeating_patterns = ;
     %concatenated features
      featMatrixaveburstDIV(file).concatenated_time = concatenated_time;
      featMatrixaveburstDIV(file).tot_bursts_DIV = tot_bursts_DIV;
     %graph theory features of concatenated bursts
      featMatrixaveburstDIV(file).density_ave_burst = density_ave;
      featMatrixaveburstDIV(file).active_nodes_ave_burst = length(nnab);
      %featMatrixaveburstDIV(file).ave_weight_ave_burst = ave_weight_ave_burst;
      featMatrixaveburstDIV(file).degree_realvsRN_ave_burst = realvsRN2;
      featMatrixaveburstDIV(file).eigenvector_mean_ratio_ave_burst = eigenvector_mean_ratio2;
      featMatrixaveburstDIV(file).eigenvector_max_ratio_ave_burst = eigenvector_max_ratio2;
      featMatrixaveburstDIV(file).num_rich_club_nodes_ave_burst = num_rich_club_nodes2;
      featMatrixaveburstDIV(file).num_hubs_ave_burst = num_hubs2;
      featMatrixaveburstDIV(file).num_hubs_RN_ratio_ave_burst = num_hubs_RN_ratio2;
      featMatrixaveburstDIV(file).assortativity_ratio_ave_burst = assortativity_ratio2;
      featMatrixaveburstDIV(file).modularity_ratio_ave_burst = modularity_ratio2;
      
      save('featureMatrixConcatenated', 'featMatrixaveburstDIV'); 
      writetable(struct2table(featMatrixaveburstDIV), fullfile(pwd,...
          'featMatrixaveburstDIV.xlsx'));
     
     elseif length(nnwc) > 6 && length(nnab) > 6
        

      %graph metrics for concatenated burst
      %density: kden = density %N = number of vertices %K = number of edges
      [kdenconc,Nconc,Kconc] = density_und(conc_thr_matrix);
       
      density_conc = (Kconc)/((num_nnwc*(num_nnwc-1))/2);
    
         
     [realvsRN, eigenvector_mean_ratio, eigenvector_max_ratio,...
         num_rich_club_nodes, num_hubs ,num_hubs_RN_ratio, assortativity_ratio,...
         modularity_ratio] = nullstats(conc_thr_matrix);
     
        %features for conc burst matrix vs random
      featMatrixDIV(file).ID = files(file).name(4:14); 
      featMatrixDIV(file).DIV = str2num(files(file).name(19:20)); 
      featMatrixDIV(file).genotype = files(file).name(1:2);
     %spike features
      featMatrixDIV(file).ave_firing_rate = ave_firing_rate;
      featMatrixDIV(file).node_similarity = node_similarity;
      featMatrixDIV(file).edge_similarity = edge_similarity ;
      %featMatrixDIV(file).mutual_info = ;
      %featMatrixDIV(file).n_repeating_patterns = ;
     %concatenated features
      featMatrixDIV(file).concatenated_time = concatenated_time;
      featMatrixDIV(file).tot_bursts_DIV = tot_bursts_DIV;
     %graph theory features of concatenated bursts
      featMatrixDIV(file).density_conc = density_conc;
      featMatrixDIV(file).active_nodes_conc = length(nnwc);
      %featMatrixDIV(file).ave_weight = ave_weight_conc;
      featMatrixDIV(file).degree_realvsRN_conc = realvsRN;
      %featMatrixDIV(file).deg_distribution_test = hrealvsRN;
      %featMatrixDIV(file).pvaluerealvsRN = pdistributions;
      %featMatrixDIV(file).size_max_connected_ratio = ;
      %featMatrixDIV(file).num_connected_ratio = ;
      %featMatrixDIV(file).mean_deg_ratio = ;
      %featMatrixDIV(file).sd_deg_ratio = ;
      %featMatrixDIV(file).max_deg_ratio = ;
      %featMatrixDIV(file).local_edge_density_ratio = ;
      featMatrixDIV(file).eigenvector_mean_ratio_conc = eigenvector_mean_ratio;
      featMatrixDIV(file).eigenvector_max_ratio_conc = eigenvector_max_ratio;
      featMatrixDIV(file).num_rich_club_nodes_conc = num_rich_club_nodes;
      featMatrixDIV(file).num_hubs_conc = num_hubs;
      featMatrixDIV(file).num_hubs_RN_ratio_conc = num_hubs_RN_ratio;
      featMatrixDIV(file).assortativity_ratio_conc = assortativity_ratio;
      featMatrixDIV(file).modularity_ratio_conc = modularity_ratio;
      
       save('featureMatrixConcatenated', 'featMatrixDIV'); 
       writetable(struct2table(featMatrixDIV), fullfile(pwd,...
           'featureMatrixConcatenated.xlsx'));
         
          
       %turn ave_burst_matrix into binary for comparison 
  
   
      
      %ave deg ave matrix
     [kdenaveburst,Naveburst,Kaveburst] = density_und(thr_ave_burst_matrix);
    
     density_ave = (Kaveburst)/((num_nnab*(num_nnab-1))/2);
     
     [realvsRN2, eigenvector_mean_ratio2, eigenvector_max_ratio2,...
         num_rich_club_nodes2, num_hubs2 ,num_hubs_RN_ratio2, assortativity_ratio2,...
         modularity_ratio2] = nullstats(thr_ave_burst_matrix);
     
      %features for ave burst matrix vs random
      featMatrixaveburstDIV(file).ID = files(file).name(4:14); 
      featMatrixaveburstDIV(file).DIV = str2num(files(file).name(19:20)); 
      featMatrixaveburstDIV(file).genotype = files(file).name(1:2);
     %spike features
      featMatrixaveburstDIV(file).ave_firing_rate = ave_firing_rate;
      featMatrixaveburstDIV(file).node_similarity = node_similarity;
      featMatrixaveburstDIV(file).edge_similarity = edge_similarity ;
      %featMatrixDIV(file).mutual_info = ;
      %featMatrixDIV(file).n_repeating_patterns = ;
     %concatenated features
      featMatrixaveburstDIV(file).concatenated_time = concatenated_time;
      featMatrixaveburstDIV(file).tot_bursts_DIV = tot_bursts_DIV;
     %graph theory features of concatenated bursts
      featMatrixaveburstDIV(file).density_ave_burst = density_ave;
      featMatrixaveburstDIV(file).active_nodes_ave_burst = length(nnab);
      %featMatrixaveburstDIV(file).ave_weight_ave_burst = ave_weight_ave_burst;
      featMatrixaveburstDIV(file).degree_realvsRN_ave_burst = realvsRN2;
      featMatrixaveburstDIV(file).eigenvector_mean_ratio_ave_burst = eigenvector_mean_ratio2;
      featMatrixaveburstDIV(file).eigenvector_max_ratio_ave_burst = eigenvector_max_ratio2;
      featMatrixaveburstDIV(file).num_rich_club_nodes_ave_burst = num_rich_club_nodes2;
      featMatrixaveburstDIV(file).num_hubs_ave_burst = num_hubs2;
      featMatrixaveburstDIV(file).num_hubs_RN_ratio_ave_burst = num_hubs_RN_ratio2;
      featMatrixaveburstDIV(file).assortativity_ratio_ave_burst = assortativity_ratio2;
      featMatrixaveburstDIV(file).modularity_ratio_ave_burst = modularity_ratio2;
      
      save('featureMatrixave_burst', 'featMatrixaveburstDIV'); 
      
      if isnan(featMatrixaveburstDIV(file).assortativity_ratio_ave_burst) == 0
        writetable(struct2table(featMatrixaveburstDIV), fullfile(pwd, 'featMatrixaveburstDIV.xlsx'));
      end
      
     end
     
    %save all plots as in folder 
      figHandles6 = findobj('Type', 'figure');
      for kk = 1:length(figHandles6)
         printformat6 = "%s figure %d .png";
         str6 = sprintf(printformat6, name, kk);
         saveas(figHandles6(kk),str6);
      end
      close all     
    
    
end 

    % find similarity between each average matrix over the DIVs


end


