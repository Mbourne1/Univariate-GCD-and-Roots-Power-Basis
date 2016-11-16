function [root_mult_array] = o_roots_matlab(fx)
% Given a vector of coefficients of a polynomail in Bernstein form. Return
% the set of roots given by zheng MultRoots() function.
%
% % Inputs.
%
% fx : Column vector of coefficients of the polynomial f(x) 
%
% % Outputs.
%
% roots_mult_array : Roots and multiplicities of f(x) as calculated by 
% the zheng MultRoots function.
%
%
%


% Build the vector of corresponding binomial coefficients
roots_calc = roots(flipud(fx));

[nEntries,~] = size(roots_calc);

root_mult_array = [roots_calc(:,1) ones(nEntries,1)];

% Printout roots to screen
PrintRoots(root_mult_array,'MATLAB METHOD');




end

