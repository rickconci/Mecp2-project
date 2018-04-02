function plottoMEA(adj_matrix)

x = 1:8;
y = 1:8;
[XX, YY] = meshgrid(x,y);

P = [XX(:), YY(:)];

[x,y] = adjacency_plot_und(adj_matrix, P);
plot(x,y)


end
