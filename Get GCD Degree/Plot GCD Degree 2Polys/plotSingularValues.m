
function plotSingularValues(arr_SingularValues)
%
% % Inputs
%
% arr_SingularValues : Array of Singular values of each S_{k}

figure_name = sprintf('%s : All Singular Values of S_{k} ');
figure('name', figure_name);
hold on

for i = 1:1:length(arr_SingularValues)
   
    % Get vector of singular values of S_{i}
    vSingularValues = arr_SingularValues{i};
    
    vec_i = i.*ones(length(vSingularValues));
    
    plot(vec_i, log10(vSingularValues), '*')
    
end
hold off


end