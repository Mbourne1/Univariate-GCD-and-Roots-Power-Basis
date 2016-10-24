function [t] = GetGCDDegree_MultipleSubresultants(vMinimumSingularValues,degree_limits,degree_limits_comp)
% Get the problem type, dependent on the vector of singular values from the
% series s_{k}
%
% % Inputs
%
% vMinimumSingularValues
%
% degree_limits
%
% bool_PossCoprime 
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
lower_lim_comp = degree_limits_comp(1);
upper_lim_comp = degree_limits_comp(2);



% Set global variables
global SETTINGS

% Get the index of the largest change in minimum singular values
[maxChangeSingularValues, indexMaxChange] = Analysis(vMinimumSingularValues);

MessageToConsole( sprintf('Largest Change in Singular Values : %4.5e' ,abs(maxChangeSingularValues)));
MessageToConsole( sprintf('Threshold : %4.5e', SETTINGS.THRESHOLD));


if (abs(maxChangeSingularValues) < SETTINGS.THRESHOLD)
    bool_significant_change = false;
else
    bool_significant_change = true;
end


% if the largest of the changes is below a threshold, then the largest
% change is not significant, and cannot determine the degree of the GCD. so
% all subresultants are either singular or nonsingular. The GCD is either 0
% or k where k is the upper limit.
if  (bool_significant_change == false)
    

    % % Rank Deficient = Singular = GCD = maximum
    % % Full Rank = Non-Singular => GCD = 0
    % maxChange is insignificant
    % Get the average minimum singular value
    avgMinSingularValue = log10(mean(vMinimumSingularValues));
    MessageToConsole( sprintf('Threshold to dermine rank : %e \n', SETTINGS.THRESHOLD_RANK));
    MessageToConsole( sprintf('Average of Singular values: %2.4f \n',avgMinSingularValue) );
    
  
    
    if  avgMinSingularValue < SETTINGS.THRESHOLD_RANK 
        % If all singular values are close to zero, then the Sylvester matrices
        % are all rank deficient, and are all singular
        % gcd is min(m,n)
        t = upper_lim_comp;
        fprintf([mfilename ' : ' calling_function ' : ' sprintf('All Rank Deficient \n')]);
        fprintf([mfilename ' : ' calling_function ' : ' sprintf('t = %i \n',t)]);
    
    elseif lower_lim_comp > 1
            t = lower_lim_comp;
            fprintf([mfilename ' : ' calling_function ' : ' sprintf('t = %i \n',t)]);
        
    else
        t = 0;
        fprintf([mfilename ' : ' calling_function ' : ' sprintf('All Full Rank \n')]);
        fprintf([mfilename ' : ' calling_function ' : ' sprintf('t = %i \n',t)]);
        
    end
else
    t = lower_lim_comp + indexMaxChange - 1;
    % maxChange is signifcant
    fprintf([mfilename ' : ' calling_function ' : ' 'min < Deg(GCD) < max \n'])
    fprintf([mfilename ' : ' calling_function ' : ' sprintf('t = %i \n',t)])
    
end



end

function [] = MessageToConsole(str)
[St,~] = dbstack();
calling_function = St(3).name;


fprintf([mfilename ' : ' calling_function ' : ' str '\n'])
end