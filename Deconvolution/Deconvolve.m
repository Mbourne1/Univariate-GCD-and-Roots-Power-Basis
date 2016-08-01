function [h] = Deconvolve(f,g)

% Get degree of polynomial f(x).
m = GetDegree(f);

% Get degree of polynomial g(x).
n = GetDegree(g);

% Build the matrix C(g
C_g = BuildT1(g,m-n);

% Solve C(g)*h = f
h = SolveAx_b(C_g,f);

end
