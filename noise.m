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

% Set the random number generator.
rng(SEED)

% Get the degree of the polynomial f(x)
m = GetDegree(fx);

% Generate a vector of random numbers between -1 and 1
a = -1;
b = 1;
r = a + (b-a).*rand(m+1,1);

% Get the vector of noise
noise_vec = fx.*(r.*mu);

% Add noise to exact coefficients
f_noisy = fx + noise_vec;


end