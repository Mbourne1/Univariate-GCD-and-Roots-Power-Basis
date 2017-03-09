function t = GetGCDDegree_2Polys(fx, gx, limits)
% GetGCDDegree(fx,gx)
%
% Get the degree of the GCD d(x) of f(x) and g(x), by Sylvester matrix method.
%
% Inputs.
%
% fx : Coefficients of polynomial f(w)
%
% gx : Coefficietns of polynomial g(w)
%
% degree_limits :
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

myLimits = [0 min(m,n)];

myLowerLimit = myLimits(1);
myUpperLimit = myLimits(2);

% Get the number of subresultants which must be constructed.
nSubresultants = myUpperLimit - myLowerLimit +1 ;

% Initialise a vector to store the minimum distnaces for each S_{k}.
vMinimumResidual = zeros(1,nSubresultants);



arr_SingularValues = cell(nSubresultants,1);

arr_R = cell(nSubresultants,1);
arr_R1 = cell(nSubresultants,1);
arr_Q = cell(nSubresultants,1);
arr_Tf = cell(nSubresultants,1);
arr_Tg = cell(nSubresultants,1);
arr_Sk = cell(nSubresultants,1);

% Set the initial value of k to be the lower limit of the possible degree
% of the GCD.
k = myLowerLimit;

% Build the Sylvester Matrix
arr_Tf{1} = BuildT1(fx,n-k);
arr_Tg{1} = BuildT1(gx,m-k);
arr_Sk{1} = [arr_Tf{1} arr_Tg{1}];

% Get QR Decomposition of S_k(f,g)
[arr_Q{1},arr_R{1}] = qr(arr_Sk{1});



% For each possible value of k, k = 0,...,min(m,n)
for i = 2 : 1 : nSubresultants
    
    k = myLowerLimit + (i-1);
    
    % update C_f and C_g by removing rows and columns
    arr_Tf{i} = arr_Tf{i-1}(1:m+n-k+1,1:n-k+1);
    arr_Tg{i} = arr_Tg{i-1}(1:m+n-k+1,1:m-k+1);
    
    % Update S_{k}
    arr_Sk{i} = [arr_Tf{i} arr_Tg{i}];
    
    % Perform QR Decomposition of Sk, by QR delete.
    % Remove the last column
    [Q,R] = qrdelete(arr_Q{i-1}, arr_R{i-1}, m+n+2-((2*k)-2), 'col');
    
    % Remove last column of C_{1}(f)
    [Q,R] = qrdelete(Q,R,n+2-k,'col');
    
    % Remove last row
    [Q,R] = qrdelete(Q,R,m+n+2-k,'row');
    
    arr_Q{i} = Q;
    arr_R{i} = R;
    
    
    
    % Add to the vector of minimum distances
    vMinimumResidual(i) = GetMinDistance(arr_Sk{i});
    
    
    
    
    
end

% Metric used to compute the degree of the GCD
% This value is changed in the SetGlobalVariables file
%   * Singular Values
%   * R1 Row Norms
%   * R1 Row Diagonals
%   * Residuals

global SETTINGS

fprintf('Metric to compute GCD degree : %s \n', SETTINGS.METRIC);

switch SETTINGS.METRIC
    
    case 'Singular Values'
        
        % Initialise a vector to store the minimum singular values for each
        % S_{k}.
        vMinimumSingularValues = zeros(nSubresultants,1);
        
        for i = 1:1:nSubresultants
            % Add to the vector of minimum Singular values from SVD of S_{k}.
            arr_SingularValues{i} = svd(arr_Sk{i});
            % Get the minimum Singular value from SVD of S_{k}
            vMinimumSingularValues(i) = min(arr_SingularValues{i});
            
        end
        
        if (SETTINGS.PLOT_GRAPHS)
            plotSingularValues(arr_SingularValues,  myLimits, limits);
            plotMinimumSingularValues(vMinimumSingularValues, myLimits, limits);
        end
        
        metric = vMinimumSingularValues;
        
    case 'Residuals'
        
        metric = vMinimumResidual;

        if (SETTINGS.PLOT_GRAPHS)
            plotMinimumResiduals(vMinimumResidual, myLimits, limits);
        end

    case 'R1 Row Diagonals'
        
        arr_DiagonalsR1 = cell(nSubresultants,1);
        vMaxDiagR1      = zeros(1,nSubresultants);
        vMinDiagR1      = zeros(1,nSubresultants);
        
        for i = 1:1:nSubresultants
            
            % Take absolute values of R_{k}
            abs_R = abs(arr_R{i});
            
            % Get number of rows in R1_{k}
            [~,nColsR1] = size(abs_R);
            
            % Obtain R1 the top square of the |R| matrix.
            arr_R1{i} = abs_R(1:nColsR1,1:nColsR1);
            
            % Get the diagonal values in the matrix R_{k} from the QR decomposition
            % of S_{k}
            arr_DiagonalsR1{i} = diag(arr_R1{i});
            
            % Get maximum and minimum row diagonals of R1
            vMaxDiagR1(i) = max(arr_DiagonalsR1{i});
            vMinDiagR1(i) = min(arr_DiagonalsR1{i});
            
        end
        
        if (SETTINGS.PLOT_GRAPHS)
            % Plot all diagonals of R1_{k} for k = 1,...,min(m,n)
            plotRowDiagonals(arr_DiagonalsR1, myLimits, limits)
            plotMaxMinRowDiagonals(vMaxDiagR1, vMinDiagR1, myLimits, limits);
        end

        metric = vMinDiagR1./vMaxDiagR1;
        
    case 'R1 Row Norms'
        
        % Initialise some vectors
        arr_R1_RowNorms = cell(nSubresultants,1);
        vMaxRowNormR1   = zeros(1,nSubresultants);
        vMinRowNormR1   = zeros(1,nSubresultants);
        
        for i = 1:1:nSubresultants
            % Take absolute values of R_{k}
            abs_R = abs(arr_R{i});
            
            % Get number of rows in R1_{k}
            [~,nColsR1] = size(abs_R);
            
            % Obtain R1 the top square of the |R| matrix.
            arr_R1{i} = abs_R(1:nColsR1,1:nColsR1);
            % Get Norms of each row in the matrix R1
            arr_R1_RowNorms{i} = sqrt(sum(arr_R1{i}.^2,2))./norm(arr_R1{i});
            
            % Get maximum and minimum row norms of rows of R1.
            vMaxRowNormR1(i) = max(arr_R1_RowNorms{i});
            vMinRowNormR1(i) = min(arr_R1_RowNorms{i});
        end
        
        if (SETTINGS.PLOT_GRAPHS)
            plotRowNorms(arr_R1_RowNorms, myLimits, limits);
            plotMaxMinRowSum(vMaxRowNormR1, vMinRowNormR1, myLimits, limits);
        end
        
        metric = vMinRowNormR1./vMaxRowNormR1;
        
    otherwise
        error('err');
end


% If only one subresultant exists, use an alternative method.
if (myUpperLimit == myLowerLimit ) % If only one Subresultant Exists
    
    if (myLowerLimit == 1)
        % Use the singular values from the only subresultant S_{1} to determine
        % if S_{1} is full rank or rank deficient.
        
        t = GetGCDDegree_OneSubresultant(metric);
    else
        % Since can not be corpime, and only one subresultant exists
        
        t = myLowerLimit;
        fprintf([mfilename ' : ' sprintf('One subresultant : t = %i \n',t)]);
    end
    
else
    % Get the type of problem.
    % Problem Type.
    % Singular      : All Subresultants S_{k} are Singular, and rank deficient
    % NonSingular   : All Subresultants S_{k} are Non-Singular, and full rank
    % Mixed         : Some Subresultants are Singular, others are Non-Singular.
    
    [t] = GetGCDDegree_MultipleSubresultants(metric, myLimits);
    
end





end






