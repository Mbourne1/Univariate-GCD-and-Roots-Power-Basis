function plotMinimumSingularValues(vMinimumSingularValues, limits_k, limits_t, rank_range)
%
% % Inputs 
% 
% vMinimumSingularValues : (Vector)
%
% limits_k : (Int) (Int) : Limits on the degree of the GCD. Defines the range
% of k values 
%
% limits_t : [Int Int] : Prior computed limits for the degree of the GCD
%
% rank_range : [Float Float] 
%


% 
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

%
lowerLimit_t = limits_t(1);
upperLimit_t = limits_t(2);

%
rank_range_low = rank_range(1);
rank_range_high = rank_range(2);

%
x_vec = lowerLimit_k : 1 : upperLimit_k;

figure_name = sprintf('%s : Minimum Singular Values of S_{k}',mfilename);
figure('name',figure_name);
hold on
plot(x_vec, log10(vMinimumSingularValues), '-s', 'DisplayName', 'Singular Values')
ylabel('log \sigma(k)')
xlabel('k')
legend(gca,'show');

% Plot vertical lines
vline(lowerLimit_t)
vline(upperLimit_t)

% plot horizontal lines
hline(rank_range_low);
hline(rank_range_high);

hold off

end