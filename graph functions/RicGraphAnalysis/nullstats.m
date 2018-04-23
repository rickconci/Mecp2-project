function [realvsRN, eigenvector_mean_ratio, eigenvector_max_ratio, num_rich_club_nodes, num_hubs ,num_hubs_RN_ratio, assortativity_ratio, modularity_ratio] = nullstats(adj_matrix)

%plot adjacency matrix
G = graph(adj_matrix,'OmitSelfLoops');

[sOut,tOut] = findedge(G);
nodesinnetwork = unique(vertcat(unique(sOut), unique(tOut)));
num_nodes_real = length(nodesinnetwork);

%average, max degree
[degreal] = degrees_und(adj_matrix);
ave_deg_real = mean(degreal);
std_deg_real = std(degreal);
var_de_real = std_deg_real^2;
max_deg_real = max(degreal);

%crate a null-model with Watts Strogtz rewiring algorithm
WS = WattsStrogatz(num_nodes_real,round(ave_deg_real/2),0.15);
RN = WattsStrogatz(num_nodes_real, round(ave_deg_real/2), 1);
WS_adj = full(adjacency(WS));
RN_adj = full(adjacency(RN));


%compare degree distributions between burstadjmatrix, WS and random
[degWS] = degrees_und(WS_adj);
ave_deg_WS = mean(degWS);
std_deg_WS = std(degWS);
var_de_WS = std_deg_WS^2;
max_deg_WS = max(degWS);
      
[degRN] = degrees_und(RN_adj);
ave_deg_RN = mean(degRN);
std_deg_RN = std(degRN);
var_de_RN = std_deg_RN^2;
max_deg_RN = max(degRN);
      
figure
h1 = histogram(degreal);
title("degree distribution of network burst")

figure
h2 = histogram(degRN);
title("degree distribution of random burst")


[realvsWS, hrealvsWS] = ranksum(degreal,degWS);
[realvsRN, hrealvsRN]  = ranksum(degreal,degRN);
[WSvsRN, hWSvsRN] = ranksum(degWS,degRN);

%distributions = [degreal', full(degRN)']; %full(degWS)',
%pdistributions = kruskalwallis(distributions);

%size_max_connected_ratio
%[CIJscore,sn] = score_wu(adj_matrix,s);


%num_connected_ratio, 


%local_edge_density_ratio, 


%eigenvector_ratio
eigreal = eigenvector_centrality_und(adj_matrix);
eigRN = eigenvector_centrality_und(RN_adj);
eigenvector_mean_ratio = mean(eigreal)/mean(eigRN);
eigenvector_max_ratio = max(eigreal)/max(eigRN);

%number of hubs
[hub_nodes_real, nodesinnetwork_real]  = findHubs(adj_matrix);
%[hub_nodes_WS, nodesinnetwork_WS]  = findHubs(WS_adj);
[hub_nodes_RN, nodesinnetwork_RN]  = findHubs(RN_adj);
num_hubs = length(hub_nodes_real);
num_hubs_RN_ratio = length(hub_nodes_real)/length(hub_nodes_RN);       
      

%Rich Club Topology
%R = vector of rich-club coefficients for levels 1 to klevel.
%Nk = number of nodes with degree>k
%Ek = number of edges remaining in subgraph with degree>k
[RWS,NWS,EkWS] = rich_club_bu(WS_adj, max_deg_WS);
[Rrnd,Nrnd,Ekrnd] = rich_club_bu(RN_adj, max_deg_RN);
if mean(nonzeros(adj_matrix(:))) < 1
    [Rreal] = rich_club_wu(adj_matrix,max_deg_real);
    Ekrndratio = (Ekrnd/max(Ekrnd));
    if length(Ekrndratio) < length(Rreal)
        dif1 = length(Rreal) - length(Ekrndratio);
        Ekrndratio(length(Ekrndratio) + 1:length(Ekrndratio)+ dif1) = 0;
    elseif length(Ekrndratio) > length(Rreal)
        dif2 = length(Ekrndratio) - length(Rreal);
        Ekrndratio(1: 1+dif2) = [];
    end
    for r = 1:length(Rreal)
        richclubratioRN(r) = Rreal(r)/Ekrndratio(r);
    end
    [sortedratio, Iratio] = sort(richclubratioRN);
    for nratio = 1:length(sortedratio)
        if sortedratio(nratio) >= 1
           break
        end
        [richclubnodes] = Iratio(nratio:end);
    end 
    try 
        num_rich_club_nodes = length(richclubnodes);
    catch 
        warning("'Undefined function or variable 'richclubnodes'.")
         num_rich_club_nodes = 0;
    end
elseif mean(nonzeros(adj_matrix(:))) == 1
    [Rreal,Nkreal,Ekreal] = rich_club_bu(adj_matrix,max_deg_real);
    if length(Rrnd) < length(Rreal)
        dif1 = length(Rreal) - length(Rrnd);
        Rrnd(length(Rrnd) + 1:length(Rrnd)+ dif1) = 0;
    elseif length(Rrnd) > length(Rreal)
        dif2 = length(Rrnd) - length(Rreal);
        Rrnd(1: 1+dif2) = [];
    end
    for r = 1:length(Rreal)
        try 
            richclubratioRN(r) = Rreal(r)/Rrnd(r);
        catch 
            warning("Index exceeds matrix dimensions.")
            richclubnodes = [];
            num_rich_club_nodes = 0;
            break
        end 
    end
    [sortedratio, Iratio] = sort(richclubratioRN);
    for nratio = 1:length(sortedratio)
        if sortedratio(nratio) >= 1
           break
        end
        [richclubnodes] = Iratio(nratio:end);
    end
    try 
        num_rich_club_nodes = length(richclubnodes);
    catch 
        warning("'Undefined function or variable 'richclubnodes'.")
         num_rich_club_nodes = 0;
    end
end

%Assortativity
if mean(nonzeros(adj_matrix(:))) < 1
    assortativity_real = assortativity_wei(adj_matrix,0);
    assortativity_RN = assortativity_bin(RN_adj,0);
    assortativity_ratio = assortativity_real/assortativity_RN;
elseif mean(nonzeros(adj_matrix(:)))== 1 
    assortativity_real = assortativity_bin(adj_matrix,0);
    assortativity_RN = assortativity_bin(RN_adj,0);
    assortativity_ratio = assortativity_real/assortativity_RN;
end

%modularity ratio
[Cireal,Qreal]=modularity_und(adj_matrix,1.1);
[CiRN,QRN]=modularity_und(RN_adj,1.1);
modularity_ratio = Qreal/QRN;

end

