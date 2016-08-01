function [fx_new,gx_new,alpha_new,theta_new] = ...
    SNTLN(fx_n,gx_n,alpha,theta,t)
% Given two input polynomials and the degree of their GCD, Obtain the Low
% Rank Approximation Sylvester Matrix

% Globals
global SETTINGS


% Get degree of polynomial f(x)
m = GetDegree(fx_n);

% Get degree of polynomial g(x)
n = GetDegree(gx_n);

% Set the initial iteration number
ite = 1;

% Create the identity matrix I, such that S*I = S
I = eye(m+n-2*t+2);


% Get f(w) from f(x)
fw = GetWithThetas(fx_n,theta);

% Get partial f(w) wrt \theta
fw_wrt_theta = Differentiate_wrt_theta(fw,theta);

% Get partial f(w) wrt \alpha
fw_wrt_alpha = zeros(m+1,1);

% Get g(w)
gw = GetWithThetas(gx_n,theta);

% Get partial g(w) wrt \theta
gw_wrt_theta = Differentiate_wrt_theta(gw,theta);

% Get partial g(w) wrt \alpha
gw_wrt_alpha = gw;

% Get the matrix T, used to construct S_{t}, where S_{t} = DTQ.
T = BuildT(fw,alpha.*gw,t);

% Get the index of the optimal colummn for removal
[~,opt_col] = GetMinDistance(T);

% Create the matrix M, such that S_{t}*M = A_{t}
M = I;
M(:,opt_col) = [];

% Create the matrix e, such that S_{t}*e = c_{t}
e = I(:,opt_col);

% Get the partial T wrt alpha
T_wrt_alpha = BuildT(fw_wrt_alpha,gw_wrt_alpha,t);

% Get partial T wrt theta
T_wrt_theta = BuildT(fw_wrt_theta, alpha.*gw_wrt_theta,t);

% Get the Matrix A_{t} = Sk with column c_{t} removed
At = T*M;

% Get the column c_{t}, the removed column from S_{t}
ct = T*e;

% Get the column partial c_{t} wrt \alpha
ct_wrt_alpha = T_wrt_alpha(:,opt_col);

% Get the column partial c_{t} wrt \theta
ct_wrt_theta = T_wrt_theta(:,opt_col);

% Initialise the vector of perturbations z, corresponding to polynomials
% f and g.
z = zeros(m+n+2,1);
z_fx = zeros(m+1,1);
z_gx = zeros(n+1,1);

% Build the matrix N where N has the same structure as T
N = BuildT(z_fx,alpha.*z_gx,t);

% Get the partial derivative of N_{t}
N_wrt_theta = zeros(m+n-t+1, m+n-(2*t)+2);
N_wrt_alpha = zeros(m+n-t+1, m+n-(2*t)+2);

% Get the column h_{t} - s.t h_{t} has the same structure as c_{t}
ht = N*e;
h_wrt_alpha     = N_wrt_alpha*e;
h_wrt_theta     = N_wrt_theta*e;

% Build the matrix T + N
TN = T + N;

TN_wrt_alpha = T_wrt_alpha + N_wrt_alpha;
TN_wrt_theta = T_wrt_theta + N_wrt_theta;

% Calculate the initial estimate of x - the vector whcih contains the
% coefficients of the quotient polynomials u and v.
x_ls = SolveAx_b(At,ct);

% Build the matrix P_{t} where P_{t}*[f;g] = c_{t}
Pt = BuildP(m,n,t,alpha(ite),theta(ite),opt_col);

% Build the matrix Y such that Y*[f;g] = A*x
Yt = BuildY(m,n,t,x_ls,opt_col,alpha,theta);

% Create the matrix C for input into iteration

H_z = GetH_z(opt_col,n,t,Yt,Pt,alpha(ite));

H_x = TN*M;

H_alpha = TN_wrt_alpha*M*x_ls - ...
    (ct_wrt_alpha + h_wrt_alpha);

H_theta = TN_wrt_theta*M*x_ls - ...
    (ct_wrt_theta + h_wrt_theta);

C = [H_z H_x H_alpha H_theta];

% Define the starting vector for the iterations for the LSE problem.
start_point     =   ...
    [...
    z;...
    x_ls;...
    alpha(ite);...
    theta(ite);
    ];

yy  = start_point;

%%
% Create matrix E.
E = eye(2*m+2*n-2*t+5);

% Set the initial value of vector p to be zero
% f = -(yy - start_point);
f = zeros(2*m+2*n-2*t+5,1);

% Get the initial residual vector
res_vec = (ct+ht) - (At*x_ls);

% Get residual
condition = norm(res_vec) ./ norm(ct);


while condition(ite) > SETTINGS.MAX_ERROR_SNTLN && ite < SETTINGS.MAX_ITE_SNTLN
    
    % Perfome LSE Calculation min|Ey-f| subject to Cy=g
    y_lse = LSE(E,f,C,res_vec);
    
    % Increment the iteration number
    ite = ite + 1;
    
    % Add the small changes found in LSE problem to existing values
    yy = yy + y_lse;
    
    % obtain the small changes
    delta_zk        = y_lse(1:m+n+2,1);
    delta_xk        = y_lse((m+n+3):(2*m+2*n-2*t+3),1);
    delta_alpha     = y_lse(2*m+2*n-2*t+4);
    delta_theta     = y_lse(2*m+2*n-2*t+5);
    
    % Update variables z_{k}, x_{k}, where z_{k} are perturbations in the
    % coefficients of f and g. x_{k} is the solution vector, containing
    % coefficients u and v.
    z = z + delta_zk;
    
    %x_ls = x_ls + delta_xk;
    x_ls = SolveAx_b(TN*M,ct+ht);
    
    % Update \alpha and \theta
    alpha(ite) = alpha(ite-1) + delta_alpha;
    theta(ite) = theta(ite-1) + delta_theta;
    
    % Obtain polynomials in modified basis f(w,\theta)
    fw = GetWithThetas(fx_n,theta(ite));
    gw = GetWithThetas(gx_n,theta(ite));
    
    % Construct the subresultant matrix of T.
    T = BuildT(fw,alpha(ite).*gw,t);
    
    % Get new c_{t}
    ct = T*e;
    
    % Calculate the partial derivatives of fw and gw with respect to alpha
    fw_wrt_alpha    = zeros(m+1,1);
    gw_wrt_alpha    = gw;
    
    % Calculate the partial derivatives of fw and gw with respect to theta
    fw_wrt_theta = Differentiate_wrt_theta(fw,theta(ite));
    gw_wrt_theta = Differentiate_wrt_theta(gw,theta(ite));
    
    % Calculate the partial derivative of T wrt alpha
    T_wrt_alpha = BuildT(fw_wrt_alpha,gw_wrt_alpha,t);
            
    % Calculate the partial derivative of T wrt theta
    T_wrt_theta = BuildT(fw_wrt_theta, alpha(ite).*gw_wrt_theta,t);
    

    % Calculate the derivatives of c_{k} with respect to \alpha and \theta
    ct_wrt_alpha     = T_wrt_alpha*e;
    ct_wrt_theta     = T_wrt_theta*e;
    
    % Create the vector of structured perturbations zf and zg applied
    % to F and G.
    z_fx      = z(1:m+1);
    z_gx      = z(m+2:end);
    
    % Calculate the subresultant matrix of the structured perturbations.
    z_fw = GetWithThetas(z_fx,theta(ite));
    z_gw = GetWithThetas(z_gx,theta(ite));
    
    % Calculate the derivatives of z_fw and z_gw with repect to alpha.
    partial_zfw_wrt_alpha    = zeros(m+1,1);
    partial_zgw_wrt_alpha    = z_gw;
    
    % Calculate the derivatives of z_fw and z_gw with respect to theta.
    partial_zfw_wrt_theta = Differentiate_wrt_theta(z_fw,theta(ite));
    partial_zgw_wrt_theta = Differentiate_wrt_theta(z_gw,theta(ite));   
        
    % Build the Coefficient Matrix N, of structured perturbations, with
    % same structure as T.
    N = BuildT(z_fw,alpha(ite).*z_gw,t);
    
    % Build the Sylvester matrix 
    N_wrt_alpha = BuildT(partial_zfw_wrt_alpha, partial_zgw_wrt_alpha,t);
    
    % Build the Syvlester matrix
    N_wrt_theta = BuildT(partial_zfw_wrt_theta, alpha(ite).*partial_zgw_wrt_theta,t);
    
    % update h - the vector of structured perturbations equivalent to ck
    ht = N*e;
    
    % Calculate the derivative of h with respect to alpha
    h_alpha = N_wrt_alpha*e;
    
    % Calculate the derivative of h with respect to theta
    h_theta = N_wrt_theta*e;
    
    % Update T+N
    TN = T + N;
    
    % Update parital derivative of T+N wrt alpha
    TN_alpha = T_wrt_alpha + N_wrt_alpha;
    
    % Update parital derivative of T+N wrt theta
    TN_theta = T_wrt_theta + N_wrt_theta;
    
    % Calculate the matrix DY where Y is the Matrix such that E_{k}x = Y_{k}z.
    Yt = BuildY(m,n,t,x_ls,opt_col,alpha(ite),theta(ite));
    
    % Calculate the matrix P where P is the matrix such that c = P[f;g]
    Pt = BuildP(m,n,t,alpha(ite),theta(ite),opt_col);
       
    % Calculate the residual g
    res_vec = (ct+ht) - ((TN * M) * x_ls);
    
    % Create the matrix C. This is made up of four submatrices, HZ, Hx,
    % H_alpha and H_theta.
    
    Hz      = GetH_z(opt_col,n,t,Yt,Pt,alpha(ite));
    
    Hx      = TN*M;
    
    H_alpha = TN_alpha*M*x_ls - (ct_wrt_alpha + h_alpha);
    
    H_theta = TN_theta*M*x_ls - (ct_wrt_theta + h_theta);
    
    C = [Hz,Hx,H_alpha,H_theta];  % the matrix C
    
    % Calculate the normalised residual of the solution.
    condition(ite) = norm(res_vec) ./ norm(ct+ht) ;
    
    % Update fnew - used in LSE Problem.
    f = -(yy-start_point);
    
    
end

fprintf('Required Number of iterations : %i \n',ite)

if condition(ite) < condition(1)
    
    fx_new = fx_n + z_fx;
    gx_new = gx_n + z_gx;
    theta_new = theta(ite);
    alpha_new = alpha(ite);
else
    fx_new = fx_n;
    gx_new = gx_n;
    theta_new = theta(1);
    alpha_new = alpha(1);
end

switch SETTINGS.PLOT_GRAPHS
    case 'y'
        figure('name','SNTLN - Residuals')
        hold on
        plot(log10(condition))
        hold off
        
        figure('name','SNTLN - thetas')
        hold on
        plot(log10(theta))
        hold off
    case 'n'
    otherwise
        err('err');
end

if ite == SETTINGS.MAX_ITE_SNTLN
    % revers
    fx_new = fx_n;
    gx_new = gx_n;
    theta_new = theta(1);
    alpha_new = alpha(1);
    
end


end


function Hz = GetH_z(mincol,n,d,DY,DP,alpha)
if mincol<=(n-d+1)
    Hz=(DY-DP);
else
    Hz=DY-(alpha.*DP);
end
end

function Y = BuildY(m,n,t,x_w,opt_col,alpha,theta)
% Build the Matrix Y such that Y(x)*z = E(z)*x


% Insert 0 into the position of the optimal column in the x vector


x_a = x_w(1:opt_col-1);
x_b = x_w(opt_col:end);

% Get the vector x(w), where x includes thetas
x_w = [x_a ;0 ;x_b];

% Get the vectors u(\omega,\theta) and v(\omega,\theta)

% The matrix Y is such that S_{k}(f,g)*x = Y(x)*[f;g]

% Get the x values corresponding to v_{k}
x_uw = x_w(1:n-t+1);

% Get the x values corresponding to u_{k}
x_vw = x_w(n-t+2:end);


Y2 = BuildY1(x_vw,n,theta);
Y1 = BuildY1(x_uw,m,theta);

Y = [Y1 alpha.*Y2];


end

function Y1 =  BuildY1(zv,m,theta)

% get degree of zv = n-t
n_t = length(zv)-1;

% Initialise the matrix Y1
Y1 = zeros(m+n_t+1,m+1);

% for each column
for i = 0:1:m
    Y1(i+1:n_t+1+i,i+1) = zv;
end


Y1 = Y1 * diag(theta.^(0:1:m));
end

function P = BuildP(m,n,t,alpha,theta,opt_col)
% The matrix P is such that a column c_{t} can be expressed as P_{t}[f;g]
% Given a column of the Sylvester matrix S_{t}(f,g), obtain a decomposition
% so that it is expressed as a matrix vector product where the vector
% contains only coefficients of f and g.

% P has a different structure depending on whether the column is from the
% first or second partition of S_{t}(f,g)

% Suppose the column is from the first partition, then P has the structure
% [P1 0]. There are n-k+1 columns in the first partiton
if opt_col <= n-t+1
    
    P2 = zeros(m+n-t+1,n+1);
    
    % P1 is a diagonalisation of a vector given by [zeros; ones; zeros]
    num_zero_rows_top = opt_col-1;
    num_zero_rows_bottom = (m+n-t+1) - (m+1) - num_zero_rows_top;
    P1 = ...
        [
        zeros(num_zero_rows_top,m+1);
        diag(theta.^(0:1:m));
        zeros(num_zero_rows_bottom,m+1);
        ];
    
    
    % Suppose the column is from the second partitno, then P has the structure
    % [0 P2].
else
    P1 = zeros(m+n-t+1,m+1);
    opt_col_n = opt_col - (n-t+1);
    
    num_zero_rows_top = opt_col_n-1;
    num_zero_rows_bottom = (m+n-t+1) - (n+1) - num_zero_rows_top;
    P2 = ...
        [
        zeros(num_zero_rows_top,n+1);
        diag(theta.^(0:1:n)');
        zeros(num_zero_rows_bottom,n+1)];
    
    
    
end

P = [P1 alpha.*P2];


end