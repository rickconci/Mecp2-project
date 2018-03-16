function findHubs(adj_matrix)

%find nodes with highest centrality measures and highest degree
figure("Visible", "off");
G = graph(full(adj_matrix));
h = plot(G,'Layout','force', 'EdgeColor','k');
axis tight manual
ax = gca;
ax.Visible = 'off';
ax.NextPlot = 'replaceChildren';
ax.Visible = 'off';

%node with max degree
[deg] = degrees_und(adj_matrix);
[Maxdeg, Indexdeg] = max(deg);
highlight(h,Indexdeg,'NodeColor','green');

%node with max betweenness centrality
BC = betweenness_bin(adj_matrix);
[Maxbet, Indexbet] = max(BC);
highlight(h,Indexbet, 'NodeColor','cyan');

%node with max eigenvector centrality 
v = eigenvector_centrality_und(adj_matrix);
[Maxeig, Indexeig] = max(v);
highlight(h, Indexeig, 'NodeColor', 'magenta');

%highlight minimum spanning tree
[T,p] = minspantree(G);
highlight(h,T,'EdgeColor','r','LineWidth',1.5)

end
