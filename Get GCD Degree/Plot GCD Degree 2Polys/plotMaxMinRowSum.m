function plotMaxMinRowSum(vMaxRowNormR1, vMinRowNormR1, limits)
%
% % Inputs
%
% vMaxRowNormR1
%
% vMinRowNormR1
%
% limits = [lowerLimit upperLimit]

[lowerLimit] = limits(1);
[upperLimit] = limits(2);

x = lowerLimit : 1 : upperLimit;

% Plot Graph of ratio of max : min row sum in R1 from the QR decompositions.
figure_name = sprintf('%s : Max:min Row Norms',mfilename);
figure('name',figure_name)

vRatio_MaxMin_RowNorm_R = (vMinRowNormR1 ./vMaxRowNormR1);

plot(x,log10(vRatio_MaxMin_RowNorm_R),'red-s');
hold on

legend('Min:Max Row Norms of Rows in R1 from the QR decomposition of S_{k}');
title('Min:Max Row Norms of Rows in R1 from the QR Decomposition of S_{k} (Original)');
xlim([1 upperLimit]);
xlabel('k')
ylabel('log_{10} Min/Max Row Norm')

vline(lowerLimit,'b','');
vline(upperLimit,'b','');

hold off
end