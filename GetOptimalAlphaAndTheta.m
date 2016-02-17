function [alpha, theta] = GetOptimalAlphaAndTheta(fx,gx)

% Ensure that f(x) and g(x) are column vectors
if size(fx,2) >1  || size(gx,2)>2
    error('f(x) and g(x) must be column vectors')
end

f = [1 -1 0  0];

% Get degree of polynomial f
[r,~] = size(fx);
m = r-1;

% Get degree of polynomial g
[r,~] = size(gx);
n = r-1;

% Build the first partiion

Part1 = [ones(m+1,1) zeros(m+1,1) -1.*(0:1:m)' zeros(m+1,1)];

Part2 = [ones(n+1,1) zeros(n+1,1) -1.*(0:1:n)' -1.*ones(n+1,1)];

Part3 = [zeros(m+1,1) -1.*ones(m+1,1) (0:1:m)' zeros(m+1,1)];

Part4 = [zeros(n+1,1) -1.*ones(n+1,1) (0:1:n)' ones(n+1,1)];

lambda_vec = fx;
mu_vec = gx;
rho_vec = fx;
tau_vec = gx;

A = -[Part1; Part2; Part3; Part4];


b = -[log10(abs(lambda_vec)); log10(abs(mu_vec)); -log10(abs(rho_vec));-log10(abs(tau_vec))];


x = linprog(f,A,b);

try
    theta = 10^x(3);
    alpha = 10^x(4);
    %fprintf('Optimal theta 1 and theta 2 given by: \n  theta_{1}: %0.5e \n  theta_{2}: %0.5e',theta1,theta2)
catch
    fprintf('Failed to optimize\n')
    theta = 1;
    alpha = 1; 

end

end