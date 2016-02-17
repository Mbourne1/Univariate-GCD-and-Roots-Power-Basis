function [roots_Bb] = o_roots_multroot(fx)
% Given a vector of coefficients of a polynomail in Bernstein form. Return
% the set of roots given by zheng MultRoots() function.
%
% % Inputs.
%
%   fx : Column vector of coefficients of the polynomial f(x) 
%
% % Outputs.
%
%   roots_Bb : Roots of f(x) as calculated by the zheng MultRoots function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% Get the degree of polynomial f(x) 
[r,~] = size(fx);
m = r-1;

% Build the vector of corresponding binomial coefficients
addpath('multroot')
addpath('multroot/multroot')
roots_calc = multroot(fx);

% convert roots wrt t, Bernstein basis
roots_Bb = [roots_calc(:,1) roots_calc(:,2)];

% Printout roots to screen
PrintRoots(roots_Bb,'MULTROOTS METHOD');




end

