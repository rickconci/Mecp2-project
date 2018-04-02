function real_adj_matrix = get_real_adj_matrix(ass_matrix)

%remove 1s and turn into to the mean of the data
ass_matrix(isnan(ass_matrix) == 1) = 0;
mean_of_burst = mean(ass_matrix(:));
ass_matrix(ass_matrix == 1) = mean_of_burst;
%plot as heatmap 
%figure
%ass_heatmap = heatmap(ass_matrix,'Colormap', jet, 'Title', 'Association Matrix');

% turn into vector
assvec = ass_matrix(:);

%to test significance, assume that a covariance with random results would
%give a normal distribution with mean 0, and std same as burst distribution
%(same number of data points)
%compare distribution to normal with mean of 0 and same std as
%distribution and zero values which are within the 95th % of the normal
%distribution 

%calculate the standard deviation 
std_of_burst = std(assvec);

%plot the distribution of the burst correlation values as histogram
%figure
%histogram(unique(assvec));
%xlabel('correlation coefficient');
%ylabel('frequency');

a = std_of_burst; %standard deviation
b = 0; %set mean 
normal_dist = a.*randn(length(assvec),1) + b;
%hold on 
%histogram(normal_dist);


%add a line to deptict the 95th percentile of the normal distribution 
%hold on;
%line([prctile(normal_dist, 95),prctile(normal_dist, 95)], ylim, 'LineWidth', 1, 'Color', 'k');
%hold off
for ii = 1:length(assvec)
    if assvec(ii)< prctile(normal_dist, 95)
        assvec(ii) = 0;
    else
        assvec(ii) = 1;
    end
end

%reshape the vector so as to become a matrix again

[row, col] = size(ass_matrix);
real_adj_matrix = reshape(assvec, [row,col]);


%plot matrix as heatmap
%figure
%adj_heatmap = heatmap(real_adj_matrix, 'Colormap', jet, 'Title', 'Adjacency Matrix');

%turn NaN into 0s for future analysis
real_adj_matrix(isnan(real_adj_matrix) == 1) = 0;


end 
