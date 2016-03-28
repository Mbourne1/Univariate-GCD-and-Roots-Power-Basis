function [h] = Deconvolve(q)
% Deconvolve the set of polynomials q{i}

% Get the number of polynomials in q
[~,c] = size(q);
nPolys = c;

for i = 1:1:nPolys-1
   h{i} = Deconvolve_Separate(q{i},q{i+1}) ;
end

end

function [h] = Deconvolve_Separate(f,g)

% Get degree of polynomial f(x).
m = GetDegree(f);

% Get degree of polynomial g(x).
n = GetDegree(g);

% Build the matrix C(g
C_g = BuildT1(g,m-n);

% Solve C(g)*h = f
h = SolveAx_b(C_g,f);

end
