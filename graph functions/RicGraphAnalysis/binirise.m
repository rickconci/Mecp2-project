function binary_matrix = binirise(ass_matrix)

ass_matrix(ass_matrix>0) = 1;
ass_matrix(ass_matrix<=0) = 0;
binary_matrix = ass_matrix;

end
