function t = GetGCDDegree2(fx,gx,deg_limits)
%  GetGCDDegree(fx,gx)
% 
% Get the degree of the GCD d(x) of f(x) and g(x), by Sylvester matrix method.
%
% Inputs.
%
% fw : Coefficients of polynomial f(w)
%
% gw : Coefficietns of polynomial g(w)
%
% nDistinctRoots : Number of distinct roots in previous computation.
%
% Outputs.
%
% t : Degree of GCD of f(x) and g(x)



% Global Variables
global PLOT_GRAPHS

% Get Degree of polynomial f(x)
m = GetDegree(fx);

% Get degree of polynomial g(x)
n = GetDegree(gx);

% If the number of distinct roots in f(x) is one, then the degree of the
% GCD of f(x) and f'(x) = m-1 = n.
lower_lim = deg_limits(1);
upper_lim = deg_limits(2);

if (upper_lim == lower_lim)
    t = upper_lim;
    return;
end

% if the lower limit is zero, ie - deg(GCD) = 0, then set to 1, since we
% are only interested in Sylvester matrices S_{1},...
if lower_lim == 0
    lower_lim = 0;
    possible_coprime = true;
else 
    possible_coprime = false;
end

% Get the number of subresultants which must be constructed.
nSubresultants = upper_lim - lower_lim +1 ;

% Initialise a vector to store the minimum singular values for each
% S_{k}.

vMinimumSingularValues = zeros(1,nSubresultants);


% Initialise a vector to store the minimum distnaces for each S_{k}.
vMinimumDistance = zeros(1,nSubresultants);
vMax_Diag_R1 = zeros(1,nSubresultants);
vMin_Diag_R1 = zeros(1,nSubresultants);

vMax_R1_RowNorm = zeros(1,nSubresultants);
vMin_R1_RowNorm = zeros(1,nSubresultants);

matrix = [];

% Set the initial value of k to be the lower limit of the possible degree
% of the GCD.
k = lower_lim;

% Build the Sylvester Matrix
C_f = BuildT1(fx,n-k);
C_g = BuildT1(gx,m-k);
Sk = [C_f C_g];

% Get QR Decomposition of S_k(f,g)
[Q,R] = qr(Sk);

% For each possible value of k, k = 0,...,min(m,n)
for k = lower_lim:1:upper_lim
    
    i = k - lower_lim + 1;
    
    % If not the first subresultant, build by removing rows and columns.
    if i > 1
        
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
    
    [Q,R] = qr(Sk);
    
    % Take absolute values of R_{k}
    abs_R = abs(R);

    % Get number of rows in R1_{k}
    [nRowsR1,~] = size(diag(abs_R));

    % Obtain R1 the top square of the |R| matrix.
    R1 = abs_R(1:nRowsR1,1:nRowsR1);    
    
    % Get the diagonal values in the matrix R_{k} from the QR decomposition
    % of S_{k}
    vDiagsR1 = diag(R1);
    vDiagsR1_norm = vDiagsR1 ./ norm(vDiagsR1);
    
    % Get the number of diagonals in R_{k}
    nDiagsR1 = size(vDiagsR1,1);
    
    % Save the diagonal values into a matrix of all diagonals of R_{k} for
    % each k
    matrix = [matrix ;k.*ones(nDiagsR1,1) vDiagsR1];
    
    % Add to the vector of minimum Singular values from SVD of S_{k}.
    vSingularValues = svd(Sk);
    
    % Get the minimum Singular value from SVD of S_{k}
    vMinimumSingularValues(i) = min(vSingularValues);
    
    % Add to the vector of minimum distances
    vMinimumDistance(i) = GetMinDistance(Sk);
    
    % Get maximum and minimum row diagonals of R1
    vMax_Diag_R1(i) = max(vDiagsR1);
    vMin_Diag_R1(i) = min(vDiagsR1);
    
    % Get Norms of each row in the matrix R1
    vR1_RowNorms = sqrt(sum(R1.^2,2))./norm(R1);
    
    % Get maximum and minimum row norms of rows of R1.
    vMax_R1_RowNorm(i) = max(vR1_RowNorms);
    vMin_R1_RowNorm(i) = min(vR1_RowNorms);
    
end

% Normalise the vector of minimal singular values

% Get ratio of max:min diagonals of R1
vRatio_MaxMin_Diagonals_R = vMax_Diag_R1 ./ vMin_Diag_R1;
vRatio_MaxMin_Diagonals_R = sanitize(vRatio_MaxMin_Diagonals_R);


vRatio_MaxMin_RowNorm_R = vMax_R1_RowNorm ./ vMin_R1_RowNorm;


% If only one subresultant exists, use an alternative method.
if (upper_lim -lower_lim == 0 ) % If only one Subresultant Exists
     t = GetRank_One_Subresultant(vSingularValues);
    return;
end

% Get the type of problem.
% Problem Type.
% Singular      : All Subresultants S_{k} are Singular, and rank deficient
% NonSingular   : All Subresultants S_{k} are Non-Singular, and full rank
% Mixed         : Some Subresultants are Singular, others are Non-Singular.
%% Plot data
switch PLOT_GRAPHS
    case 'y'
        
        x = lower_lim:1:upper_lim;
        
        figure_name = sprintf('%s : Plot QR Scatter Data',mfilename);
        figure('name',figure_name)
        hold on
        scatter(matrix(:,1),log10(matrix(:,2)))
        hold off
        
        % Plot the minimum singular values of S_{k} k=1,...,min(m,n)
        figure_name = sprintf('%s : Plot Min Sing Val',mfilename);
        figure('name',figure_name)
        hold on
        plot(x,log10(vMinimumSingularValues),'-o')
        hold off
        
        % Plot the minimum distances of S_{k} k=1,...,min(m,n)
        figure_name = sprintf('%s : Plot Min Distance',mfilename);
        figure('name',figure_name)
        hold on
        plot(x,log10(vMinimumDistance) ,'-s')
        hold off
        
        figure_name = sprintf('%s : Max:min Row Diagonals',mfilename);
        figure('name',figure_name)
        
        plot(x,log10(vRatio_MaxMin_Diagonals_R),'red-s');
        hold on
        axis([lower_lim,upper_lim,0,inf])
        legend('Max:Min diag element of subresultant S_{k}');
        title('Max:Min diagonal elements of R1 from the QR decomposition of S_{k} (Original)');
        ylabel('log_{10} max:min diag element')
        hold off
        
        % Plot Graph of ratio of max : min row sum in R1 from the QR decompositions.
        figure_name = sprintf('%s : Max:min Row Norms',mfilename);
        figure('name',figure_name)
        
        plot(x,log10(vRatio_MaxMin_RowNorm_R),'red-s');
        hold on
        axis([lower_lim,upper_lim,0,inf])
        legend('Max:Min Row Norms of Rows in R1 from the QR decomposition of S_{k}');
        title('Max:Min Row Norms of Rows in R1 from the QR Decomposition of S_{k} (Original)');
        hold off
        
        
    case 'n'
    otherwise
        error('error: plot_graphs either y or n')
        
end

[ProblemType,t] = GetProblemType(vMinimumSingularValues,lower_lim);


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


function [ProblemType , t] = GetProblemType(vMinimumSingularValues,lower_lim)
% Get the problem type, dependent on the vector of singular values from the
% series s_{k}
%
% Get the type of problem.
% Problem Type.
% Singular      : All Subresultants S_{k} are Singular, and rank deficient
% NonSingular   : All Subresultants S_{k} are Non-Singular, and full rank
% Mixed         : Some Subresultants are Singular, others are Non-Singular.



global THRESHOLD

min_mn = lower_lim + length(vMinimumSingularValues) - 1;

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
    
    if  avgMinSingularValue < -6
        % If all singular values are close to zero, then rank deficient, degree of
        % gcd is min(m,n)
        ProblemType = 'NonSingular';
        fprintf('NonSingular \n')
        t = min_mn;
    else
        % if all singular values are not close to zero, then full rank, degree
        % of gcd is 0
        ProblemType = 'Singular';
        fprintf('Singular \n')
        t = 0;
    end
else
    % maxChange is signifcant
    ProblemType = 'Mixed';
    fprintf('Mixed \n')
    t = lower_lim + indexMaxChange - 1;
end

end


