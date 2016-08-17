function [h] = Deconvolve(f,g)
% DECONVOLVE Given two polynomails f(x,y) and g(x,y), computes the
% polynomail division f(x)/g(x) = h(x)
%
% % Inputs.
%
% f : Coefficients of polynomial f(x)
%
% g : Coefficients of polynomial g(x)
%
% % Outputs.
%
% h : Coefficients of polynomial h(x)


% Get degree of polynomial f(x).
m = GetDegree(f);

% Get degree of polynomial g(x).
n = GetDegree(g);

% Build the matrix C_{m-n}(g)
C_g = BuildT1(g,m-n);

% Solve C_{m-n}(g)*h = f
h = SolveAx_b(C_g,f);

end
