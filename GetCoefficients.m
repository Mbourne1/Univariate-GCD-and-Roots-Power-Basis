function [fx] = GetCoefficients(roots_f)
% Given the roots and multiplicities of a polynomial f(x), calculate the
% coefficients
%
%   Inputs.
%
%   roots_f : Matrix consisting of rows of [root multiplicity] pairs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the number of distinct roots
[nRoots,~] = size(roots_f);

% Initialise the polynomials
Poly = 1;

% for each distinct root
for i = 1:1:nRoots
    
    % Get the root
    r = roots_f(i,1);
    
    % Get the multiplicity
    m = roots_f(i,2);
    
    % Multiply poly by the factor (x-r) m times
    for j = 1:1:m
        temp = [1;-r];
        Poly = conv(Poly,temp);
    end
    
end

fx = Poly;

fx = flipud(Poly);
end