
function plotMinimumSingularValues(vMinimumSingularValues, myLimits, limits)
%
% % Inputs 
% 
% vMinimumSingularValues :
%
% myLimits
%
% limits

myLowerLimit = myLimits(1);
myUpperLimit = myLimits(2);

lowerLimit = limits(1);
upperLimit = limits(2);



x_vec = myLowerLimit : 1 : myUpperLimit;

figure_name = sprintf('%s : Minimum Singular Values of S_{k}',mfilename);
figure('name',figure_name);
hold on
plot(x_vec, log10(vMinimumSingularValues), '-s', 'DisplayName', 'Singular Values')
ylabel('log \sigma(k)')
xlabel('k')
legend(gca,'show');
vline(lowerLimit)
vline(upperLimit)
hold off

end