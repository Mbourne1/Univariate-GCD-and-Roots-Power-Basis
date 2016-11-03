function P = BuildP(idx_Col,m,n,k)
% BuildP(idx_Col,m,n,t)
%
% Build the matrix P_{t}, where h_{t} = P_{t}z
%
% % Inputs.
%
% idx_Col : Index of column removed from Sylvester matrix S_{k}, and will
% be constructed using the vector P.
%
% m : Degree of polynomial f(x)
%
% n : Degree of polynomial g(x)
%
% k : Index of Sylvester matrix S_{k}

if idx_Col <= n-k+1 % Column is in first partition of S_{k}(f,g)
    
    i = idx_Col;
    
    P = ...
        [
        zeros(i-1,m+1)      zeros(i-1,n+1);
        eye(m+1,m+1)        zeros(m+1,n+1);
        zeros(n-k-i+1,m+1)  zeros(n-k-i+1,n+1);
        ];
    
else % Column is in second partition of S_{k}(f,g)
    
    % Get index with respect to the second partiton of the sylvester matrix
    % only.
    i = idx_Col - (n-k+1);
    
    P = ...
        [
        zeros(i-1,m+1)      zeros(i-1,n+1);
        zeros(n+1,m+1)      eye(n+1,n+1)
        zeros(m-k-i+1,m+1)  zeros(m-k-i+1,n+1);
        ];
    
end

end