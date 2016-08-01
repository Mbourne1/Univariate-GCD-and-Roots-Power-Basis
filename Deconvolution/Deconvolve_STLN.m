function [hx] = Deconvolve_STLN(f,g)
% Get the polynomial h = f/g

% Get degree of polynomial f(x).
m = GetDegree(f);

% Get degree of polynomial g(x).
n = GetDegree(g);

% Build the matrix C(g
C_g = BuildT1(g,m-n);

% Solve C(g)*h = f
hx = SolveAx_b(C_g,f);

% Build the vector z(x), corresponding to perturbations of coefficients of
% g(x).
z = zeros(n+1,1);

% Build the Cauchy matrix of coefficients of z(x)
C_z = BuildT1(z,m-n);

% Build the matrix C(h) where C(h) * g = C(g) * h
Y = BuildT1(hx,n);

% Build the matrix H = [Y (C+E) -I]
[ Y (C_g + C_z) -eye(m+1)]

% condition
condition(ite) = norm(res_vec) ./ norm(f + s 

end
