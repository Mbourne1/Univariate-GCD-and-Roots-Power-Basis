function [h] = PolyAdd(f,g)

% if nCols > 1, not a column vector
if size(f,2) >1 || size(g,2) >1
   error('Not a column vector') 
end

% Get degree of polynomial f(x)
[r,~] = size(f);
m = r-1;

% Get degree of polynomial g(x)
[r,~] = size(g);
n = r-1;

if m < n
    f = [f; zeros(n-m,1)];
elseif n < m
    g = [g; zeros(m-n,1)];
end

h = f+g;

end