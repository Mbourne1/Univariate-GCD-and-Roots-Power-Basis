
function plotRowNorms(arr_RowNorms, limits)
% 
% % Inputs
%
% arr_RowNorms
%
% limits : 
lower_lim = limits(1);
upper_lim = limits(2);

figure_name = sprintf('%s : Diag Norm',mfilename);
figure('name',figure_name)
hold on
for i = 1:1:length(arr_RowNorms)
    
    
    % get vector of row norms
    vec_RowNorms = arr_RowNorms{i};
    vec_i = i.*ones(length(vec_RowNorms));
    
    plot(vec_i, log10(vec_RowNorms),'*');

end

hold off

xlabel('k')
ylabel('log10 Row Norm of R1 from QR decomposition of S_{k}')
title(sprintf('log10 Row Norm of R1 from the QR decomposition of each subresultant S_{k}'));
xlim([1 upper_lim]);

vline(lower_lim,'b','');
vline(upper_lim,'b','');
hold off
end





