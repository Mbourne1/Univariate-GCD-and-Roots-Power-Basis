function [] = o_roots_myNewMethod(ex_num,noise,seed)

global PLOT_GRAPHS 
PLOT_GRAPHS = 'y';

global SEED 
SEED = seed;

global MAX_ERROR_SNTLN 
MAX_ERROR_SNTLN = 1e-13;

global MAX_ITE_SNTLN
MAX_ITE_SNTLN = 500;

% Given Input Polynomial f
fx_exact = Examples_Roots(ex_num);

% Add Noise to the coefficients of f
fx = Noise(fx_exact,noise);

gx = Differentiate(fx);

display(fx);
display(gx);

% Get degree of polynomial f
[r,~] = size(fx);
m = r - 1;

% Get degree of polynomial g
n = m - 1;

% Get the degree t of the GCD of f(x) and g(x)
t = GetDegree(fx,gx);

% Build the Sylvester Matrix
Cf = BuildC1(fx,n,t);
Cg = BuildC1(gx,m,t);
St = [Cf Cg];

display(St);

% Get optimal column with minimal distance Ax = b
[~, col] = GetMinDistance(St);

% Get the matrix A_{t}
At = St;
At(:,col) = [];
% Get the column c_{t}
ct = St(:,col);

display(At)
display(ct)

% Obtain low rank approximation
fx_output = STLN(fx,t);
gx_output = Differentiate(fx_output);

C1_noisy = BuildC1(fx,n,1);
C2_noisy = BuildC1(gx,m,1);
S1_noisy = [C1_noisy C2_noisy];

C1_lowrank = BuildC1(fx_output,n,1);
C2_lowrank = BuildC1(gx_output,m,1);
S1_lowrank = [C1_lowrank C2_lowrank];

vSingularValues_lowRank = svd(S1_lowrank);
vSingularValues_Noisy = svd(S1_noisy);

vSingularValues_lowRank = vSingularValues_lowRank./vSingularValues_lowRank(1)
vSingularValues_Noisy = vSingularValues_Noisy./vSingularValues_Noisy(1)

figure('name','Singular Values')
hold on
plot(log10(vSingularValues_lowRank),'DisplayName','Low Rank Approx');
plot(log10(vSingularValues_Noisy),'DisplayName','Noisy');
legend(gca,'show');
hold off

% Get the distance of f(x) from the exact polynomial
fprintf('Distance of f(x) noisy from f(x) exact: \n')
norm(fx - fx_exact) ./ norm(fx_exact)



% Get the distance of f(x) + \delta f(x) from the exact polynomial
fprintf('Distance of f(x) +  delta f(x) from f(x) exact: \n')
norm(fx_output - fx_exact) ./ norm(fx_exact)
