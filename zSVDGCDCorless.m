function [] = zSVDGCDCorless()

ex_num = '2';

switch ex_num
    case '1'
        
        p = [1; 5.503; 9.765; 7.647; 2.762; 0.37725];
        q = [1; -2.993; -0.7745; 2.0070; 0.7605];
    case '2'
        p = [1; 5.503; 9.765; 7.647; 2.762; 0.37725];
        PrintCoefficientsBivariate(p,'f')
        el = 1e-5;
        noise_vec = rand(size(p)) .* el;
        q = p + noise_vec
        
end
[r,~] = size(p);
m = r - 1;

[r,~] = size(q);
n = r - 1;

% Build the Sylvester Matrix
C1 = BuildC1(p,n,1);
C2 = BuildC1(q,m,1);
S_noisy = [C1 C2];
S_exact = [C1 C1];


vSingularValues_noisy = svd(S_noisy)
vSingularValues_exact = svd(S_exact)

rank(S_exact)
rank(S_noisy)

% Plot Results
figure('name','Singular Values')
plot(log10(vSingularValues_noisy))
hold on
plot(log10(vSingularValues_exact))
title('Plot Singular Values of S_{1}(f,g)')
hold off

%% Analysis of results
% Distance of S_{noisy} from S_{exact}
norm(S_exact - S_noisy) ./ norm(S_exact)





end