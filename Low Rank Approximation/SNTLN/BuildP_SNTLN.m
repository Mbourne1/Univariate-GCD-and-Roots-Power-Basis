
function P = BuildP_SNTLN(m,n,k,alpha,theta,idx_col)
% The matrix P is such that a column c_{t} can be expressed as P_{t}[f;g]
% Given a column of the Sylvester matrix S_{t}(f,g), obtain a decomposition
% so that it is expressed as a matrix vector product where the vector
% contains only coefficients of f and g.
%
% % Inputs
%
% m : Degree of polynomial f(x)
%
% n : Degree of polynomial g(x)
%
% k : Degree of polynomial d(x)
%
% alpha : Optimal value of \alpha
%
% theta : Optimal value of \theta
% 
% idx_col : Index of optimal column S_{k}(f,g)
%
% % Outputs
%
% P : Matrix P

nCols_Tf = n-k+1;

if idx_col <= nCols_Tf
    
    P2 = zeros(m+n-k+1,n+1);
    
    % P1 is a diagonalisation of a vector given by [zeros; ones; zeros]
    num_zero_rows_top = idx_col-1;
    num_zero_rows_bottom = (m+n-k+1) - (m+1) - num_zero_rows_top;
    P1 = ...
        [
        zeros(num_zero_rows_top,m+1);
        diag(ones(m+1,1));
        zeros(num_zero_rows_bottom,m+1);
        ];
    
    
    % Suppose the column is from the second partitno, then P has the structure
    % [0 P2].
else
    P1 = zeros(m+n-k+1,m+1);
    opt_col_n = idx_col - (n-k+1);
    
    num_zero_rows_top = opt_col_n-1;
    num_zero_rows_bottom = (m+n-k+1) - (n+1) - num_zero_rows_top;
    P2 = ...
        [
        zeros(num_zero_rows_top,n+1);
        diag(ones(n+1,1));
        zeros(num_zero_rows_bottom,n+1)];
    
end

% Get thetas corresponding to polynomial f(x)
th_f = diag(GetWithThetas(ones(m+1,1),theta));

% Get thetas corresponding to polynomial g(x)
th_g = diag(GetWithThetas(ones(n+1,1),theta));

% Build the matrix P
P = [P1*th_f alpha.*P2*th_g];


end