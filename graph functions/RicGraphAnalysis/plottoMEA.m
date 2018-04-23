function Graph = plottoMEA(adj_matrix, Title)
G = graph(adj_matrix,'OmitSelfLoops');
if length(adj_matrix) == 60
    figure ('visible','off')

    x_raw= repelem(1:8, 8);
    x_raw(1) = [];
    x_raw(7) = [];
    x_raw(55) = [];
    x_raw(61) = [];
    x = x_raw;

    y_rep = repmat(1:8, 8);
    y = y_rep(1,:);
    y(1) = [];
    y(7) = [];
    y(55) = [];
    y(61) = [];

    colormap jet
    [degreal] = degrees_und(adj_matrix);
    nSizes = 2*sqrt(degreal-min(degreal)+0.2);
    nColors = degreal;
    h1 = plot(G,'Layout','force', 'EdgeColor','k', 'MarkerSize',nSizes,'NodeCData',nColors,'EdgeAlpha',0.1);
    h1.XData = x;
    h1.YData = y;
    colorbar
    ax.Visible = 'off'; 
    title(sprintf(Title));
end

Graph = G;
end
