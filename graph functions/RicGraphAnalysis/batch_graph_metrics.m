function batch_graph_metrics

files = dir('*.mat');  % where your .mat files are 

for file = 1:length(files)
    data = load(files(file).name);
    name = files(file).name;
    
    %density
    %kden = density
    %N = number of vertices
    %K = number of edges
    [kden,N,K] = density_und(data.burstAdjMatrix);

    %average, max degree
    [deg] = degrees_und(data.burstAdjMatrix);
    ave_deg = mean(deg);
    max_deg = max(deg);
    %clustering coefficient 
    %C is a vector
    clustering_coefficient =clustering_coef_bu(data.burstAdjMatrix);

    %average path length
    %average_path_length = ave_path_length(file);
    %number of hubs

    %Rich Club Topology
    %R = vector of rich-club coefficients for levels 1 to klevel.
    %Nk = number of nodes with degree>k
    %Ek = number of edges remaining in subgraph with degree>k
    [R,Nk,Ek] = rich_club_bu(data.burstAdjMatrix);

    %topology of functional network in burst
    
    featMatrix(file).batch = files(file).name(4:11); 
    featMatrix(file).ID = files(file).name(4:14); 
    featMatrix(file).DIV = str2num(files(file).name(19:20)); 
    featMatrix(file).genotype = files(file).name(1:2);
    featMatrix(file).burst_number = (files(file).name(54:55));
    featMatrix(file).density = kden;
    featMatrix(file).ave_deg = ave_deg;
    featMatrix(file).max_deg = max_deg;
    featMatrix(file).clustering_coefficient = clustering_coefficient;
    %featMatrix(file).average_path_length = average_path_length;
    featMatrix(file).richclub = R;

end 
    save('featMatrixM5', 'featMatrix'); 
end
