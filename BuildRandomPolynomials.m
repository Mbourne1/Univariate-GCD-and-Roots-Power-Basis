function [f_roots,g_roots] = BuildRandomPolynomials(m,n,t,intvl_low,intvl_high)
% Build two random polynomials f(x) and g(x) with a gcd d(x)
% Note this function does not allow multiple roots
% Note this function allows roots to be 'close'.
% define t, the degree of the gcd
%
% Input
%
% m : Degree of Polynomial f(x)
%
% n : Degree of polynomial g(x)
%
% t : Degree of GCD(f,g) d(x)
%
% intvl_low : lowest root
%
% intvl_high : highest root


global SEED
rng(SEED)


% for each degree t, get a root in the interval [0,1]
a = intvl_low ;
b = intvl_high;

d_roots = zeros(t,2);

roots = a + (b-a).*rand(t,1);
d_roots = [roots ones(t,1)];

% % % Build set of roots of f

% include the roots of d
f_roots = d_roots;
roots = a + (b-a) .* rand(m-t,1);
f_roots = [d_roots; roots ones(m-t,1)]

% Build set of roots of g
g_roots = d_roots;
roots = a + (b-a) .* rand(n-t,1);
g_roots = [d_roots; roots ones(n-t,1)]

f_roots;
g_roots;


