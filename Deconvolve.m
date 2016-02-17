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

% Get degree of polynomial f(x)
[r,~] = size(f);
m = r-1;

C_g = BuildC(g,m);

h = pinv(C_g)*f;

end

function C_g = BuildC(g,m)
% Build the matrix C such that C(f)*h = g, where h is an unknown vector.
%
% Inputs
%
% g : Input polynomial
%
% m : Degree of polynomial f

% Get the degree of polynomial g(x)
[r,~] = size(g);
n = r-1;

% For each column of C(g)
for j = 0:1:m-n
   C_g(j+1:j+1+n,j+1) = g;
end

end