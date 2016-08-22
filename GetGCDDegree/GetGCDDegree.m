function t = GetGCDDegree(fx,gx,degree_limits)
% GetGCDDegree(fx,gx)
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


% Get Degree of polynomial f(x)
m = GetDegree(fx);

% Get degree of polynomial g(x)
n = GetDegree(gx);

% If the number of distinct roots in f(x) is one, then the degree of the
% GCD of f(x) and f'(x) = m-1 = n.
lower_lim = degree_limits(1);
upper_lim = degree_limits(2);

% Get the number of subresultants which must be constructed.
nSubresultants = upper_lim - lower_lim +1 ;

% Initialise a vector to store the minimum singular values for each
% S_{k}.
vMinimumSingularValues = zeros(1,nSubresultants);


% Initialise a vector to store the minimum distnaces for each S_{k}.
vMinimumResidual = zeros(1,nSubresultants);
vMaxDiagR1      = zeros(1,nSubresultants);
vMinDiagR1      = zeros(1,nSubresultants);
vMaxRowNormR1   = zeros(1,nSubresultants);
vMinRowNormR1   = zeros(1,nSubresultants);

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
    
    % Take absolute values of R_{k}
    abs_R = abs(R);
    
    % Get number of rows in R1_{k}
    [nRowsR1,~] = size(diag(abs_R));
    
    % Obtain R1 the top square of the |R| matrix.
    R1 = abs_R(1:nRowsR1,1:nRowsR1);
    
    % Get the diagonal values in the matrix R_{k} from the QR decomposition
    % of S_{k}
    vDiagsR1 = diag(R1);
    
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
    vMinimumResidual(i) = GetMinDistance(Sk);
    
    % Get maximum and minimum row diagonals of R1
    vMaxDiagR1(i) = max(vDiagsR1);
    vMinDiagR1(i) = min(vDiagsR1);
    
    % Get Norms of each row in the matrix R1
    vR1_RowNorms = sqrt(sum(R1.^2,2))./norm(R1);
    
    % Get maximum and minimum row norms of rows of R1.
    vMaxRowNormR1(i) = max(vR1_RowNorms);
    vMinRowNormR1(i) = min(vR1_RowNorms);
    
end

% Determine whether the two polynomials can be corpime from the lower
% bound.
if lower_lim == 1 
    can_be_coprime = true;
else
    can_be_coprime = false;
end

% If only one subresultant exists, use an alternative method.
if (upper_lim == lower_lim ) % If only one Subresultant Exists
    
    if (can_be_coprime)
        % Use the singular values from the only subresultant S_{1} to determine
        % if S_{1} is full rank or rank deficient.
        
        t = GetGCDDegree_OneSubresultant(vSingularValues);
    else
        % Since can not be corpime, and only one subresultant exists
       
        t = lower_lim;
        fprintf([mfilename ' : ' sprintf('One subresultant : t = %i \n',t)]);
    end
    
else
    % Get the type of problem.
    % Problem Type.
    % Singular      : All Subresultants S_{k} are Singular, and rank deficient
    % NonSingular   : All Subresultants S_{k} are Non-Singular, and full rank
    % Mixed         : Some Subresultants are Singular, others are Non-Singular.
    PlotGraphs()
    [t] = GetGCDDegree_MultipleSubresultants(vMaxDiagR1./vMinDiagR1,degree_limits);
    
end





end






