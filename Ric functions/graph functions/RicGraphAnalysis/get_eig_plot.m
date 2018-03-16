function eig_plot = get_eig_plot(eig_values)

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

z = eig_values.^2;
color = [0 0 0];
sf = 35;

BubblePlot(x,y,z,color,sf);


end
