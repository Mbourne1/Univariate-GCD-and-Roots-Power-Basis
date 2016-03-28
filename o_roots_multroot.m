function [root_mult_array] = o_roots_multroot(fx)
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
addpath('multroot')
addpath('multroot/multroot')
roots_calc = multroot(fx);


root_mult_array = [roots_calc(:,1) roots_calc(:,2)];

% Printout roots to screen
PrintRoots(root_mult_array,'MULTROOTS METHOD');




end

