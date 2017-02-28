
function plotMinimumSingularValues(vMinimumSingularValues, limits)
%
% % Inputs 
% 
% vMinimumSingularValues :
%
% 

figure_name = sprintf('%s : Minimum Singular Values of S_{k}',mfilename);
figure('name',figure_name);
hold on
xlim([1 +inf]);
plot(log10(vMinimumSingularValues),'-s','DisplayName','Singular Values')
ylabel('log \sigma(k)')
xlabel('k')
legend(gca,'show');
hold off

end