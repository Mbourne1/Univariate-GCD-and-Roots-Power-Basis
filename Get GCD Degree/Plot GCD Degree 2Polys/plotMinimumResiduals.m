
function plotMinimumResiduals(vMinimumResidual, limits)
%
% % Inputs
%
% vMinimumResidual : Vector containing minimum residuals
%
% limits : [lowerLimit : upperLimit]

lowerLimit = limits(1);
upperLimit = limits(2);

% % Plot Residuals
figure_name = sprintf('%s : Residuals',mfilename);
figure('name',figure_name);
hold on
xlim([1 +inf]);
plot(log10(vMinimumResidual),'-s','DisplayName','SVD')
ylabel('log r(k)')
xlabel('k')
legend(gca,'show');
hold off
end
