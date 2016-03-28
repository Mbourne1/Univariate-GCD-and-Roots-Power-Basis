function [h] = PolyAdd(f,g)

% if nCols > 1, not a column vector
if size(f,2) >1 || size(g,2) >1
   error('Not a column vector') 
end

% Get degree of polynomial f(x)
m = GetDegree(f);


% Get degree of polynomial g(x)
n = size(g);


if m < n
    f = [f; zeros(n-m,1)];
elseif n < m
    g = [g; zeros(m-n,1)];
end

h = f+g;

end