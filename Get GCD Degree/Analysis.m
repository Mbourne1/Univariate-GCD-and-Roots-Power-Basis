
function [max_Delta_MaxMin_Diag_R, indexMaxChange_Ratio_DiagsR] = Analysis(vRatio_MaxMin_Diagonals_R)

% Get Degree by max:min Diagonals

vRatio_MaxMin_Diagonals_R = sanitize(vRatio_MaxMin_Diagonals_R);

% Get the change in the ratios of diagonal elements from one subresultant
% to the next.
vDelta_MaxMin_Diag_R = abs(diff(log10(vRatio_MaxMin_Diagonals_R)));

% Get the maximum change in diag ratio and its index
[max_Delta_MaxMin_Diag_R, indexMaxChange_Ratio_DiagsR] = max(vDelta_MaxMin_Diag_R);


end