function [f_noisy,noise_vec] = Noise(fx,mu)
% Given a polynomial f(x), add noise to its coefficients at the level
% specified by mue.
%
%   Inputs;
%
%   fx : 
%
%   mu :
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise the global variables
global SEED

% Get the degree of the polynomial f(x)
[r,~] = size(fx);
m = r-1;

% Set the random number generator.
rng(SEED)

% Generate a vector of random numbers between -1 and 1
a = -1;
b = 1;
r = a + (b-a).*rand(m+1,1);

% Get the vector of noise
noise_vec = fx.*(rand.*mu);

% Add noise to exact coefficients
f_noisy = fx + noise_vec;


end