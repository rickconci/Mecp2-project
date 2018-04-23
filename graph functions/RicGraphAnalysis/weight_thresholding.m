function thr_matrix = weight_thresholding(ass_matrix)

ass_matrix(isnan(ass_matrix) == 1) = 0;
ass_matrix(logical(eye(size(ass_matrix)))) = 0;

weights = nonzeros(ass_matrix);
max_weight = max(weights);

%figure
%h1 = histogram(weights);


thr_matrix = ass_matrix;
thr_matrix(thr_matrix<(0.3*max_weight)) = 0;
thr_matrix(thr_matrix>0) = 1;

%figure
%h2 = heatmap((thr_matrix),'Colormap', jet);
%figure
%h2 = heatmap((ass_matrix - thr_matrix),'Colormap', jet);


end
