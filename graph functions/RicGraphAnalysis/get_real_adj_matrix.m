function real_adj_matrix = get_real_adj_matrix(ass_matrix)

%remove 1s and turn into to the mean of the data
ass_matrix(isnan(ass_matrix) == 1) = 0;
ass_matrix(logical(eye(size(ass_matrix)))) = 0;
%mean_of_burst = mean(ass_matrix(:));
%ass_matrix(logical(eye(size(ass_matrix)))) = mean_of_burst;

% turn into vector
%assvec = ass_matrix(:);


h1 = heatmap(ass_matrix,'Colormap', jet);

[str] = strengths_und(ass_matrix);
str(str==0) = 0.1;

normalised_weights = zeros(length(str));
num_edges = zeros(1, length(str));
null_int = zeros(1, length(str));
thr_matrix = zeros(length(str));

for jj = 1:length(str)
    normalised_weights(:, jj) = ass_matrix(:, jj)/str(jj);
    thr_matrix = normalised_weights;
    num_edges(1,jj) = nnz(ass_matrix(:, jj));
    null_int(1,jj) = 1/(num_edges(1,jj)-1);
    null_int(null_int==-1) = 0; 
    thr_matrix(thr_matrix(:,jj) <= null_int(1,jj)) = 0;
end

alpha = zeros(length(ass_matrix));
fun = @(x,k) (1-x).^(k-2);
for cc = 1:length(ass_matrix)
    for rr = 1:length(ass_matrix)
        alpha(rr,cc) = 1 - (num_edges(cc)-1)*integral(@(
    end
end


thr_local_matrix = (thr_matrix+thr_matrix')/2;
thr_global_matrix = normalised_weights;
h2 = heatmap(thr_matthr_local_matrixrix,'Colormap', jet);
figure
h3 = heatmap((burst_matrix1 - thr_matthr_local_matrixrix),'Colormap', jet);

G_thr = plottoMEA(thr_matrix, 'thr_matrix');
G_burst = plottoMEA(burst_matrix1, 'burst_matrix1');
%to test significance, assume that a covariance with random results would
%give a normal distribution with mean 0, and std same as burst distribution
%(same number of data points)
%compare distribution to normal with mean of 0 and same std as
%distribution and zero values which are within the 95th % of the normal
%distribution 

%calculate the standard deviation 
%std_of_burst = std(assvec);

%plot the distribution of the burst correlation values as histogram
%figure
%histogram(unique(assvec));
%xlabel('correlation coefficient');
%ylabel('frequency');

%a = std_of_burst; %standard deviation
%b = 0; %set mean 
%normal_dist = a.*randn(length(assvec),1) + b;
%hold on 
%histogram(normal_dist);


%add a line to deptict the 95th percentile of the normal distribution 
%hold on;
%line([prctile(normal_dist, 95),prctile(normal_dist, 95)], ylim, 'LineWidth', 1, 'Color', 'k');
%hold off
%for ii = 1:length(assvec)
%    if assvec(ii)< prctile(normal_dist, 95)
%        assvec(ii) = 0;
%    elseif assvec(ii)> prctile(normal_dist, 95)
%        assvec(ii) = 1;
%    end
%
%reshape the vector so as to become a matrix again

%[row, col] = size(ass_matrix);
%real_adj_matrix = reshape(assvec, [row,col]);


%plot matrix as heatmap
%figure
%adj_heatmap = heatmap(real_adj_matrix, 'Colormap', jet, 'Title', 'Adjacency Matrix');

%turn NaN into 0s for future analysis
%real_adj_matrix(isnan(real_adj_matrix) == 1) = 0;


end 
