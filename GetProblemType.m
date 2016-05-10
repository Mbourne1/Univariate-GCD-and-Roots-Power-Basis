function [t] = GetProblemType(vMinimumSingularValues,degree_limits)
% Get the problem type, dependent on the vector of singular values from the
% series s_{k}
%
% Get the type of problem.
% Problem Type.
% Singular      : All Subresultants S_{k} are Singular, and rank deficient
% NonSingular   : All Subresultants S_{k} are Non-Singular, and full rank
% Mixed         : Some Subresultants are Singular, others are Non-Singular.

% Get the name of the function which called this function
[St,~] = dbstack();
calling_function = St(2).name;

% Set upper and lower limits of the degree of the GCD.
lower_lim = degree_limits(1);
upper_lim = degree_limits(2);

global SETTINGS

% Get the differences between minimum singular value S_{k} and S_{k+1}
vDeltaMinSingVal = abs(diff(log10(vMinimumSingularValues)));

% Get the index of the largest change in minimum singular values
[maxChangeSingularValues, indexMaxChange] = max(vDeltaMinSingVal);

fprintf([mfilename ' : ' sprintf('Largest Change in Singular Values : %4.5e \n' ,abs(maxChangeSingularValues))]);
fprintf([mfilename ' : ' sprintf('Threshold : %4.5e \n', SETTINGS.THRESHOLD)]);

% if the largest of the changes is below a threshold, then the largest
% change is not significant, and cannot determine the degree of the GCD. so
% all subresultants are either singular or nonsingular. The GCD is either 0
% or k where k is the upper limit.
if  abs(maxChangeSingularValues) < SETTINGS.THRESHOLD
    

    % % Rank Deficient = Singular = GCD = maximum
    % % Full Rank = Non-Singular => GCD = 0
    % maxChange is insignificant
    % Get the average minimum singular value
    avgMinSingularValue = log10(norm(vMinimumSingularValues));
    
    if  avgMinSingularValue < SETTINGS.THRESHOLD
        % If all singular values are close to zero, then the Sylvester matrices
        % are all rank deficient, and are all singular
        % gcd is min(m,n)
        fprintf([calling_function ' : ' sprintf('All Rank Deficient \n')]);
        t = upper_lim;
    else
        
        fprintf([calling_function ' : ' sprintf('All Full Rank \n')]);
        t = 0;
    end
else
    % maxChange is signifcant
    fprintf([calling_function ' : ' 'min < Deg(GCD) < max \n'])
    t = lower_lim + indexMaxChange - 1;
end



end