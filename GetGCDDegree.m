function t = GetGCDDegree(fx,gx)
% Get the degree of the GCD d(x) of f(x) and g(x), by Sylvester matrix method.
%
%   Inputs
%
%   fw : Coefficients of polynomial f(x)
%
%   gw : Coefficietns of polynomial g(x)

global PLOT_GRAPHS

% Get Degree of polynomial f(x)
m = GetDegree(fx);

% Get degree of polynomial g(x)
n = GetDegree(gx);

% Initialise a vector to store the minimum singular values for each
% S_{k}.
vMinimumSingularValues = zeros(1,min(m,n));

% Initialise a vector to store the minimum distnaces for each S_{k}.
vMinimumDistance = zeros(1,min(m,n));

matrix = [];

% Set the initial value of k = 1.
k = 1;

% Build the Sylvester Matrix

C_f = BuildT1(fx,n-k);
C_g = BuildT1(gx,m-k);

Sk = [C_f C_g];

% Get QR Decomposition of Sk
[Q,R] = qr(Sk);


% For each possible value of k, k = 0,...,min(m,n)
for k = 1:1:min(m,n)
    
    % If not the first subresultant, build by removing rows and columns.
    if k > 1
        
        % update C_f and C_g by removing rows and columns
        C_f = C_f(1:m+n-k+1,1:n-k+1);
        C_g = C_g(1:m+n-k+1,1:m-k+1);
        
        % Update S_{k}
        Sk = [C_f C_g];
        
        % Perform QR Decomposition of Sk, by QR delete.
        % Remove the last column
        [Q,R] = qrdelete(Q,R,m+n+2-((2*k)-2),'col');
        
        % Remove last column of C_{1}(f)
        [Q,R] = qrdelete(Q,R,n+2-k,'col');
        
        % Remove last row
        [Q,R] = qrdelete(Q,R,m+n+2-k,'row');
        
    end
    
    
    % Get the diagonal values in the matrix R_{k} from the QR decomposition
    % of S_{k}
    vDiagsR = diag(R);
    
    % Get the number of diagonals in R_{k}
    nDiagsR = size(vDiagsR,1);
    
    % Save the diagonal values into a matrix of all diagonals of R_{k} for
    % each k
    matrix = [matrix ;k.*ones(nDiagsR,1) vDiagsR];
    
    % Add to the vector of minimum Singular values from SVD of S_{k}.
    vSingularValues = svd(Sk);
    
    % Get the minimum Singular value from SVD of S_{k}
    vMinimumSingularValues(k) = min(vSingularValues);
    
    % Add to the vector of minimum distances
    vMinimumDistance(k) = GetMinDistance(Sk);
    
end

% Normalise the vector of minimal singular values

vMinimumSingularValues = Normalise(vMinimumSingularValues);

% If only one subresultant exists, use an alternative method.
if min(m,n) == 1 % If only one Subresultant Exists
    
    t = GetRank_One_Subresultant(vSingularValues);
    return;
    
end

% Get the type of problem.
% Problem Type.
% Singular      : All Subresultants S_{k} are Singular, and rank deficient
% NonSingular   : All Subresultants S_{k} are Non-Singular, and full rank
% Mixed         : Some Subresultants are Singular, others are Non-Singular.

[ProblemType,t] = GetProblemType(vMinimumSingularValues);

fprintf('Problem Type : %s \n',ProblemType)

% Print the degree of the GCD
fprintf('The computed degree of the GCD is %i \n',t);

%% Plot data
switch PLOT_GRAPHS
    case 'y'
        
        figure('name','Get Degree : QR Scatter data')
        hold on
        scatter(matrix(:,1),log10(matrix(:,2)))
        hold off
        
        % Plot the minimum singular values of S_{k} k=1,...,min(m,n)
        figure('name','Get Degree : Plot Min Sing val')
        hold on
        plot(log10(vMinimumSingularValues),'-o')
        hold off
        
        % Plot the minimum distances of S_{k} k=1,...,min(m,n)
        figure('name','Get Degree : Plot Min Distance')
        hold on
        plot(log10(vMinimumDistance) ,'-s')
        hold off
        
    case 'n'
    otherwise
        error('error: plot_graphs either y or n')
        
end

end


function t = GetRank_One_Subresultant(vMinSingVal)
% Given the vector
% Get the rank, where only one subresultant exists.
global PLOT_GRAPHS

global THRESHOLD

% Only one subresultant
fprintf('Only one subresultant exists. \n')
switch PLOT_GRAPHS
    case 'y'
        figure('name','GetDegree - One Subresultant - Singular Values of S1')
        hold on
        title('Singular values of S_{1}')
        plot(log10(vMinSingVal))
        hold off
    case 'n'
end

% Get the differences between singular values of S_{1}
vDiffSingularValues = diff(log10(vMinSingVal));

% Get the index of the largest change in singular values of S_{1}
[deltaSingularValues, ~] = max(vDiffSingularValues);

% If the change is smaller than the predefined threshold value, then plot
% is considered 'flat'.
if deltaSingularValues < THRESHOLD
    
    % The subresultant is of full rank, in which case t = 0
    t = 0;
    fprintf('The only Subresultant S_{1} appears to be of NonSingular. \n');
    return
    
else % val > threshold
    
    % The subresultant S_{1} is rank deficient, in which case t = 1
    t = 1;
    fprintf('The only Subresultant S_{1} appears to be Singular \n');
    return
    
end
end


function [ProblemType , t] = GetProblemType(vMinimumSingularValues)
% Get the problem type, dependent on the vector of singular values from the
% series s_{k}
%
% Get the type of problem.
% Problem Type.
% Singular      : All Subresultants S_{k} are Singular, and rank deficient
% NonSingular   : All Subresultants S_{k} are Non-Singular, and full rank
% Mixed         : Some Subresultants are Singular, others are Non-Singular.



global THRESHOLD

min_mn = length(vMinimumSingularValues);

% Get the differences between minimum singular value S_{k} and S_{k+1}
vDeltaMinSingVal = abs(diff(log10(vMinimumSingularValues)));
vDeltaMinSingVal(isnan(vDeltaMinSingVal)) = 0 ;
vDeltaMinSingVal(vDeltaMinSingVal==inf)=0;

% Get the index of the largest change in minimum singular values
[maxChangeSingularValues, indexMaxChange] = max(vDeltaMinSingVal);

if  abs(maxChangeSingularValues) < THRESHOLD
    
    % maxChange is insignificant
    % Get the average minimum singular value
    avgMinSingularValue = log10(norm(vMinimumSingularValues));
    
    if  avgMinSingularValue > -6
        % if all singular values are close to zero, then rank deficient, degree of
        % gcd is min(m,n)
        ProblemType = 'NonSingular';
        t = min_mn;
    else
        % if all singular values are not close to zero, then full rank, degree
        % of gcd is 0
        ProblemType = 'Singular';
        t = 0;
    end
else
    % maxChange is signifcant
    ProblemType = 'Mixed';
    t = indexMaxChange;
end

end


