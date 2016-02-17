function [fx_output,gx_output,dx_output, ux_output, vx_output , alpha_output, theta_output] = o1(fx,gx)
% given two polynomials f(x) and g(x) return the gcd d(x)
global BOOL_PREPROC
global PLOT_GRAPHS
global BOOL_SNTLN


% Get degree of polynomial f(x)
[r,~] = size(fx);
m = r - 1;

% Get degree of polynomial g(x)
[r,~] = size(gx) ;
n = r - 1;

% Initialise useful vectors
vecm = (0:1:m)';
vecn = (0:1:n)';

%% Preprocess the polynomials f(x) and g(x) 
switch BOOL_PREPROC
    case 'y'
        % Apply Preprocessing
        
        % Get geometric mean
        lambda  = geomean(abs(fx(fx~=0)));
        mu      = geomean(abs(gx(gx~=0)));
        
        % Divide f(x) and g(x) by geometric mean
        fx_n = fx./ lambda;
        gx_n = gx./ mu;
        
        % Get opitmal values of alpha and theta
        [alpha, theta] = GetOptimalAlphaAndTheta(fx_n,gx_n);
        
        % Obtain f(w) and g(w) from f(x) and g(x)
        fw = fx_n.*(theta.^vecm);
        gw = gx_n.*(theta.^vecn);
        
        switch PLOT_GRAPHS
            case 'y'
                % Plot the coefficients of fx
                figure('name','Coefficients of f(x)')
                plot(fx,'-s','DisplayName','Coefficients Before Processing')
                hold on
                plot(fw,'-o','DisplayName','Coefficients After Processing')
                legend(gca,'show');
                hold off
                
                % Plot the coefficients of g(x)
                figure('name','Coefficients of g(x)')
                plot(gx,'-s','DisplayName','Coefficients Before Processing')
                hold on
                plot(gw,'-o','DisplayName','Coefficients After Processing')
                legend(gca,'show');
                hold off
            case 'n'
            otherwise
                error('err : plot graphs either y or n')
        end
        
    case 'n'
        % Dont apply preprocessing
        fx_n = fx;
        gx_n = gx;
        fw = fx;
        gw = gx;
        
        lambda = 1;
        mu = 1;
        
        theta = 1;
        alpha = 1;
        
    otherwise
        error('bool_preproc either y or n');
end

%%
% Get the degree of the gcd
t = getDegree(fw,alpha.*gw);

% Given the degree t, get the optimal column for removal from S_{t}(f,g)
C_f_preprocessed = BuildC1(fw,n,t);
C_g_preprocessed = BuildC1(gw,m,t);
St_preprocesed = [C_f_preprocessed alpha.*C_g_preprocessed];

% Get the minimum distance between any one of the columns of S_{k}(f,g) and the
% remaining columns of S_{k}(f,g)
[~,opt_col] = getMinDistance(St_preprocesed);

%%
% Perform SNTLN to obtain low rank approximation
switch BOOL_SNTLN
    case 'y'
        % Get f(x)+\delta x and g(x) + \delta x for which a low rank
        % approximate Sylvester matrix is obtained.
        [fx_new,gx_new,alpha_new,theta_new] = SNTLN(fx,gx,t,opt_col,lambda,mu,alpha,theta);
        
        fx_n = fx_new;
        gx_n = gx_new;
        theta = theta_new;
        alpha = alpha_new;
        fw_new = fx_new .* (theta.^vecm);
        gw_new = gx_new .* (theta.^vecn);
        
    case 'n'
        
        fw_new = fw;
        gw_new = gw;
    otherwise
        error('error')
end

C_f_LowRankApprox = BuildC1(fw_new,n,t);
C_g_LowRankApprox = BuildC1(gw_new,m,t);
St_LowRankApprox = [C_f_LowRankApprox alpha.*C_g_LowRankApprox];

%%

% Get the Sylvester matrix of the unprocessed
C_f_unproc = BuildC1(fx,n,t);
C_g_unproc = BuildC1(gx,m,t);
St_Unproc = [C_f_unproc C_g_unproc];

%% 
% Get Quotients from either the low rank approximation or the Sylvester 
% Subresultant of the preprocessed polynomials f(\omega) and g(\omega)

switch BOOL_SNTLN
    case 'y'
        St = St_LowRankApprox;
    case 'n'
        St = St_preprocesed;
end


% Get Quotient polynomials u(x) and v(x) such that 
% f(x)/u(x) = g(x)/v(x) = d(x)

% Get the matrix A_{t}(f,g), S_{t} with the optimal column removed
At = St;
At(:,opt_col) = [];

% Get the column c_{t} removed from S_{t}
ct = St(:,opt_col) ;

% Get the vector x from A_{t}x = c_{t}
x_ls = pinv(At) * ct;

% insert 1 into x_ls in the position corresponding to the col
x = [x_ls(1:opt_col-1); -1 ; x_ls(opt_col:end)];

% Split x into v(\omega) and u(\omega)
vw = x(1:n-t+1);
uw = -x(n-t+2:end);

% Get u(x) and v(x)
ux = uw./(theta.^((0:1:m-t)'));
vx = vw./(theta.^((0:1:n-t)'));

%%
% Get the GCD d(x) by the APF 

C_u = BuildC1(uw,t,0);
C_v = BuildC1(vw,t,0);

C1 = ...
    [
    C_u;
    C_v;
    ];

% Build the Right hand side vector [f(\omega) ; alpha.* g(\omega)]
rhs_vec = [fw;alpha.*gw];

% Calculate d(\omega)
dw = pinv(C1)*rhs_vec;

% Get d(x)
dx = dw./(theta.^(0:1:t))';

%% 
% Get the singular values for
% 1 : The Sylvester Subresultant for unprocessed f(x) and g(x)
% 2 : The Sylvester Subresultant for preprocessed f(\omega) g(\omega)
% 3 : The Sylvester Subresultant which is a low rank approximation of 2.

sing_vals_unproc = svd(St_Unproc);
sing_vals_unproc = sing_vals_unproc ./ sing_vals_unproc(1);

sing_vals_preproc = svd(St_preprocesed);
sing_vals_preproc = sing_vals_preproc ./ sing_vals_preproc(1);

sing_vals_LowRankApprox = svd(St_LowRankApprox);
sing_vals_LowRankApprox = sing_vals_LowRankApprox ./ sing_vals_LowRankApprox(1);


%%
% Plot Singular values of unproc, preproc, lowrank approx
switch PLOT_GRAPHS
    case 'y'
        figure('name','Singular Values')
        hold on
        plot(log10(sing_vals_unproc),'-s','DisplayName','Unprocessed')
        plot(log10(sing_vals_preproc),'-o','DisplayName','Preprocessed')
        plot(log10(sing_vals_LowRankApprox),'-*','DisplayName','Low Rank Approx')
        legend(gca,'show');
        hold off
    case 'n'
    otherwise
        error('error: plot_graphs either y or n')
end

%% Outputs
fx_output = fx_n;
gx_output = gx_n;
dx_output = dx;
ux_output = ux;
vx_output = vx;
theta_output = theta;
alpha_output = alpha;


end