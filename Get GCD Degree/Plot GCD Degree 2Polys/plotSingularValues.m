
function plotSingularValues(arr_SingularValues, myLimits, limits)
%
% % Inputs
%
% arr_SingularValues : Array of Singular values of each S_{k}
%
% myLimits 
%
% limits

%
myLowerLimit = myLimits(1);
myUpperLimit = myLimits(2);

%
nSubresultants = myUpperLimit - myLowerLimit + 1;

%
lowerLimit = limits(1);
upperLimit = limits(2);


figure_name = sprintf([mfilename ' : All Singular Values of S_{k}']);
figure('name', figure_name);
hold on

for i = 1 : 1 : nSubresultants
   
    k = myLowerLimit + (i-1);

    % Get vector of singular values of S_{i}
    vSingularValues = arr_SingularValues{i};
    
    vec_k = k.*ones(length(vSingularValues));
    
    plot(vec_k, log10(vSingularValues), '*')
    
end
vline(lowerLimit);
vline(upperLimit);
hold off


end