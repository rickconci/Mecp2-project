function [hub_nodes, nodes] = findHubs(adj_matrix)

%ass_matrix(isnan(ass_matrix) == 1) = 0;
%mean_of_burst = mean(ass_matrix(:));
%ass_matrix(ass_matrix == 1) = mean_of_burst;
%plot as heatmap 
%figure
%ass_heatmap = heatmap(ass_matrix,'Colormap', jet, 'Title',... 
%'Association Matrix');

%adj_matrix = ass_matrix;

%% find nodes with highest centrality measures and highest degree
G = graph(adj_matrix,'OmitSelfLoops');
figure %('visible','off')
h1 = plot(G,'Layout','force', 'EdgeColor','k');
ax.Visible = 'off';


nodes = unique(G.Edges(:, 1));


%% node with max degree
[deg] = degrees_und(adj_matrix);
[Maxdeg, Indexdeg] = max(deg);
top10deg = prctile(deg,90);
[sorteddeg, Ideg] = sort(deg);
for ndeg = 1:length(sorteddeg)
    if sorteddeg(ndeg) >= top10deg
        break
    end
    degree_hubs2 = Ideg(ndeg:end);
end 
%add max to original plot
highlight(h1,Indexdeg,'NodeColor','green');
%plot bar graph over nodes
%figure
%degreebarchart = bar(sorteddeg);
%set(gca, 'XTickLabel', Ideg);
%title("Degree Centrality");

%plot degree distribution
%figure
%degree_distribution = histogram(deg);
%title("Degree Distribution")

%% closeness centrality 
closeness = centrality(G,'closeness');
top10close = prctile(closeness,90);
[sortedclose, Iclose] = sort(closeness);
for nclo = 1:length(sortedclose)
    if sortedclose(nclo) >= top10close
       break
    end
    closeness_hubs = Iclose(nclo:end);
end 
%plot centrality as bar graphs
%figure
%closebarchart = bar(sortedclose);
%set(gca, 'XTickLabel', Iclose);
%title("Closeness Centrality");
%
%figure
%h3 = plot(G,'Layout','force', 'EdgeColor','k');
%ax.Visible = 'off';
%h3.NodeCData = closeness;
%colormap jet
%colorbar
%title('Closeness Centrality Scores - Unweighted')

%% betweenness centrality
BC = betweenness_bin(adj_matrix);
[Maxbetw, Indexbetw] = max(BC);
top10betw = prctile(BC,90);
[sortedbetw, Ibetw] = sort(BC);
for nbet = 1:length(sortedbetw)
    if sortedbetw(nbet) >= top10betw
        break
    end
    betw_hubs = Ibetw(nbet:end);
end 
%highlight max node in original graph
highlight(h1,Indexbetw, 'NodeColor','cyan');

%plot centrality as bar graphs
%figure
%betwbarchart = bar(sortedbetw);
%set(gca, 'XTickLabel', Ibetw);
%title("Betweenness Centrality");

%plot with color map 
%figure
%h2 = plot(G,'Layout','force', 'EdgeColor','k');
%ax.Visible = 'off';
%h2.NodeCData = BC;
%colormap jet
%colorbar
%title("Betweenness Centrality with Color");

%% eigenvector centrality 
v = eigenvector_centrality_und(adj_matrix);
[Maxeig, Indexeig] = max(v);
top10eig = prctile(v,90);
[sortedeig, Ieig] = sort(v);
for neig = 1:length(sortedeig)
    if sortedeig(neig) >= top10eig
        break
    end
    eig_hubs = Ieig(neig:end);
end

%highlight max eigenvector node on original graph
highlight(h1, Indexeig, 'NodeColor', 'magenta');

%plot eig values as bar chart
%figure
%eigbarchart = bar(sortedeig);
%set(gca, 'XTickLabel', Ieig);
%title("Eigenvector Centrality");

%plot graph with eig color map 
%figure
%h3 = plot(G,'Layout','force', 'EdgeColor','k');
%ax.Visible = 'off';
%h3.NodeCData = v;
%colormap jet
%colorbar
%title("Eigenvector Centrality with Color");

try
    hub_nodes = intersect(intersect(intersect(degree_hubs2, closeness_hubs),eig_hubs), betw_hubs);
catch 
   %if (strcmp(ME.identifier,"'Undefined function or variable 'degree_hubs2'"))
   warning("'Undefined function or variable 'degree_hubs2'")
   hub_nodes = 0;
end




%highlight minimum spanning tree
[T,p] = minspantree(G);
highlight(h1,T,'EdgeColor','g','LineWidth',1.5)


end
