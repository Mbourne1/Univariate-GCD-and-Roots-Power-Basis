function [fx_new,gx_new,alpha_new,theta_new] = SNTLN(fx,gx,t,opt_col,lambda,mu,alpha,theta)
% Given two input polynomials and the degree of their GCD, Obtain the Low
% Rank Approximation Sylvester Matrix

global MAX_ERROR_SNTLN
global MAX_ITE_SNTLN
global PLOT_GRAPHS



% Get degree of polynomial f(x)
[r,~] = size(fx);
m = r - 1;

% Get degree of polynomial g(x)
[r,~] = size(gx);
n = r - 1;

% Set the initial iteration number
ite = 1;

% Initialise useful vectors
vecm = (0:1:m)';
vecn = (0:1:n)';
vecnk = (0:1:n-t)';
vecmk = (0:1:m-t)';



%%
% Create the identity matrix I, the matrix M formed from I by removing the
% column equivalent to the optimal column for removal from the Sylvester
% subresultant matrix, so that S_{t}(f,g)*M = A_{t}, where A_{t} is the
% Sylvester subresultant matrix with the column removed.
I = eye(m+n-2*t+2,m+n-2*t+2);

M = I;

M(:,opt_col) = [];

% Let e be the column removed from the identity matrix, such that
% S_{t}(f,g) * e gives the column c_{t}, where c_{t} is the optimal column
% removed from the Sylvester subresultant.
e = I(:,opt_col);
%%

fx_n = fx./lambda;
gx_n = gx./mu;

% Get f(w) and g(w)
fw = fx_n .* (theta.^vecm);

% Get partial fw wrt theta
partial_fw_wrt_theta = vecm.* fx .* (theta.^(-1:1:m-1)');
% Get partial fw wrt alpha
partial_fw_wrt_alpha = zeros(m+1,1);

% Get g(w)
gw = gx_n .* (theta.^vecn);
% Get partial fw wrt theta
partial_gw_wrt_theta = vecn.* gx .* (theta.^(-1:1:n-1)');
partial_gw_wrt_alpha = gw;

% Get partial fw wrt alpha
%%

% Get the Sylvester Matrix T
T1 = BuildC1(fw,n,t);
T2 = BuildC1(gw,m,t);
T = [T1 alpha.*T2];

% Get the partial Sk wrt alpha
partial_T1_wrt_alpha = BuildC1(partial_fw_wrt_alpha,n,t);
partial_T2_wrt_alpha = BuildC1(partial_gw_wrt_alpha,m,t);
partial_T_wrt_alpha = [partial_T1_wrt_alpha partial_T2_wrt_alpha];

% Get partial sk wrt theta
partial_T1_wrt_theta = BuildC1(partial_fw_wrt_theta,n,t);
partial_T2_wrt_theta = BuildC1(partial_gw_wrt_theta,m,t);
partial_T_wrt_theta = [partial_T1_wrt_theta alpha.* partial_T2_wrt_theta];

% Get the Matrix Ak - Sk with column ck removed
Ak = T*M;

% Get the column ck - Removed column from Sk
ck = T*e;

% Get the column partial ck_wrt_alpha
Partial_ck_wrt_alpha = partial_T_wrt_alpha(:,opt_col);

% Get the column partial ck_wrt_theta
Partial_ck_wrt_theta = partial_T_wrt_theta(:,opt_col);


% Get the vector of perturbations given by z
z = zeros(m+n+2,1);
zf = z(1:m+1);
zg = z(m+2:end);

% Build the matrix N - s.t Ek has the same structure as T
N1 = BuildC1(zf,n,t);
N2 = BuildC1(zg,m,t);
N = [N1 alpha.* N2];

% Get the partial derivative of Ek
Partial_N_wrt_theta = zeros(m+n-t+1, m+n-(2*t)+2);
Partial_N_wrt_alpha = zeros(m+n-t+1, m+n-(2*t)+2);

% Build the column h - s.t h has the same structure as ck
h = N*e;
Partial_h_wrt_alpha     = Partial_N_wrt_alpha*e;
Partial_h_wrt_theta     = Partial_N_wrt_theta*e;

% Get the partial derivative of Ak wrt alpha
T1_wrt_alpha = zeros(m+n-t+1,n-t+1);
T2_wrt_alpha = BuildC1(gw,m,t);
T_wrt_alpha = [T1_wrt_alpha, T2_wrt_alpha];


TN = T + N;
TN_wrt_alpha = T_wrt_alpha + Partial_N_wrt_alpha;
TN_wrt_theta = partial_T_wrt_theta + Partial_N_wrt_theta;

%%
% Calculate the initial estimate of x - the vector whcih contains the
% coefficients of the quotient polynomials u and v.
[~,n2] = size(Ak);

[Q2,R] = qr(Ak);

R1 = R(1:n2,:);

cd = Q2'*ck;

c1 = cd(1:n2,:);

% EDIT - 11/06/2015
% x_ls includes thetas corresponding to uw and vw
x_ls = R1\c1;

%%
P = BuildP(m,n,t,alpha(ite),theta(ite),opt_col);
test1 = P*[fx_n;gx_n];
test2 = ck;



% Strip the thetas from x vector
x_ls_wrt_omega = x_ls;

% Get the vector of thetas corresponding to x(\omega), with missing value.
thetas = [theta.^vecnk ; theta.^vecmk];
thetas(opt_col,:) = [];
x_ls_wrt_x = x_ls_wrt_omega./thetas;

% Build the matrix Y such that Y*[f;g] = A*x
Y = BuildY(m,n,t,x_ls_wrt_x,opt_col,alpha,theta);
test1 = Y*[fx_n;gx_n];
test2 = Ak*x_ls;

%%

% Create the matrix C for input into iteration

H_z     = GetH_z(opt_col,n,t,Y,P,alpha(ite));

H_x     = TN*M;

H_alpha  = TN_wrt_alpha*M*x_ls - ...
    (Partial_ck_wrt_alpha + Partial_h_wrt_alpha);

H_theta = TN_wrt_theta*M*x_ls - ...
    (Partial_ck_wrt_theta + Partial_h_wrt_theta);

C       = [H_z H_x H_alpha H_theta];

% Define the starting vector for the iterations for the LSE problem.
start_point     =   ...
    [...
    z;...
    x_ls;...
    alpha(ite);...
    theta(ite)
    ];

yy              =   start_point;

%%
% Create matrix E.
E = eye(2*m+2*n-2*t+5);

% Set the initial value of vector p to be zero
p = zeros(2*m+2*n-2*t+5,1);

res_vec = (ck+h) - (Ak*x_ls);

x_wrt_x = x_ls;

z_fx = zeros(m+1,1);
z_gx = zeros(n+1,1);
%%
% Get residual
condition = norm(res_vec);



while condition(ite) > MAX_ERROR_SNTLN && ite < MAX_ITE_SNTLN
    
    y = LSE(E,p,C,res_vec);
    
    % Increment the iteration number
    ite = ite + 1;
    
    % Add the small changes found in LSE problem to existing values
    yy = yy + y;
    
    % obtain the small changes
    delta_zk        = y(1:m+n+2,1);
    delta_xk        = y((m+n+3):(2*m+2*n-2*t+3),1);
    delta_alpha     = y(2*m+2*n-2*t+4);
    delta_theta     = y(2*m+2*n-2*t+5);
    
    % Update variables z_{k}, x_{k}, where z_{k} are perturbations in the
    % coefficients of f and g. x_{k} is the solution vector, containing
    % coefficients u and v.
    z = z + delta_zk;
    
    x_wrt_x = x_wrt_x + delta_xk;
    
    % Update alpha and theta
    alpha(ite) = alpha(ite-1) + delta_alpha;
    theta(ite) = theta(ite-1) + delta_theta;
    
    % Obtain polynomials in modified bersntein basis a_{i}\theta^{i}
    fw = fx_n.*(theta(ite).^vecm);
    gw = gx_n.*(theta(ite).^vecn);
    
    % Construct the subresultant matrix of T.
    T1 = BuildC1(fw,n,t);
    T2 = BuildC1(gw,m,t);
    T = [T1 alpha(ite).*T2];
    
    % Calculate the partial derivatives of fw and gw with respect to alpha
    Partial_fw_wrt_alpha    = zeros(1,m+1);
    Partial_gw_wrt_alpha    = gw;
    
    % Calculate the partial derivatives of fw and gw with respect to theta
    Partial_fw_wrt_theta    = vecm.*fw./theta(ite);
    Partial_gw_wrt_theta    = vecn.*gw./theta(ite);
    
    % Calculate the partial derivative of T wrt alpha
    Partial_T1_wrt_alpha    = 0.*BuildC1(fw,n,t);
    Partial_T2_wrt_alpha    = BuildC1(gw,m,t);
    Partial_T_wrt_alpha     = [Partial_T1_wrt_alpha Partial_T2_wrt_alpha];
    
    % Calculate the partial derivative of T wrt theta
    Partial_T1_wrt_theta    = BuildC1(Partial_fw_wrt_theta,n,t);
    Partial_T2_wrt_theta    = BuildC1(Partial_gw_wrt_theta,m,t);
    Partial_T_wrt_theta     = [Partial_T1_wrt_theta alpha(ite).* Partial_T2_wrt_theta];
    
    % Get new ck
    ck = T*e;
    
    % Calculate the derivatives of c_{k} with respect to \alpha and \theta
    Partial_ck_wrt_alpha     = Partial_T_wrt_alpha*e;
    Partial_ck_wrt_theta     = Partial_T_wrt_theta*e;
    
    % Create the vector of structured perturbations zf and zg applied
    % to F and G.
    z_fx      = z(1:m+1);
    z_gx      = z(m+2:end);
    
    % Calculate the subresultant matrix of the structured perturbations.
    z_fw     = z_fx.*(theta(ite).^vecm);
    z_gw     = z_gx.*(theta(ite).^vecn);
    
    % Calculate the derivatives of z_fw and z_gw with repect to alpha.
    Partial_zfw_wrt_alpha    = zeros(m+1,1);
    Partial_zgw_wrt_alpha    = z_gw;
    
    % Calculate the derivatives of z_fw and z_gw with respect to theta.
    Partial_zfw_wrt_theta    = vecm.*(z_fw./theta(ite));
    Partial_zgw_wrt_theta    = vecn.*(z_gw./theta(ite));
    
    % Build the Coefficient Matrix N, of structured perturbations, with
    % same structure as T.
    N1 = BuildC1(z_fw,n,t);
    N2 = BuildC1(z_gw,m,t);
    N = [N1 alpha(ite).*N2];
    
    Partial_N1_wrt_alpha = BuildC1(Partial_zfw_wrt_alpha,n,t);
    Partial_N2_wrt_alpha = BuildC1(Partial_zgw_wrt_alpha,m,t);
    Partial_N_wrt_alpha = [Partial_N1_wrt_alpha Partial_N2_wrt_alpha];
    
    Partial_N1_wrt_theta = BuildC1(Partial_zfw_wrt_theta,n,t);
    Partial_N2_wrt_theta = BuildC1(Partial_zgw_wrt_theta,m,t);
    Partial_N_wrt_theta = [Partial_N1_wrt_theta alpha(ite).*Partial_N2_wrt_theta];
    
    % update h - the vector of structured perturbations equivalent to ck
    h = N*e;
    
    % Calculate the derivative of h with respect to alpha
    h_alpha = Partial_N_wrt_alpha*e;
    
    % Calculate the derivative of h with respect to theta
    h_theta = Partial_N_wrt_theta*e;
    
    % Update T+N
    TN = T + N;
    
    % Update parital derivative of T+N wrt alpha
    TN_alpha = Partial_T_wrt_alpha + Partial_N_wrt_alpha;
    
    % Update parital derivative of T+N wrt theta
    TN_theta = Partial_T_wrt_theta + Partial_N_wrt_theta;
    
    % Calculate the matrix DY where Y is the Matrix such that E_{k}x = Y_{k}z.
    Y = BuildY(m,n,t,x_wrt_x,opt_col,alpha(ite),theta(ite));
    
    % Calculate the matrix P where P is the matrix such that c = P[f;g]
    P = BuildP(m,n,t,alpha(ite),theta(ite),opt_col);
    
    % Calculate the residual q and vector p.
    rk = (ck+h) - (TN * M * x_wrt_x);
    
    residual(ite) = norm(rk);
    
    % Create matrix E.
    E = eye(2*m+2*n-2*t+5);
    
    % Create the matrix C. This is made up of four submatrices, HZ, Hx,
    % H_alpha and H_theta.
    
    Hz      = GetH_z(opt_col,n,t,Y,P,alpha(ite));
    
    Hx      = TN*M;
    
    H_alpha = TN_alpha*M*x_wrt_x - (Partial_ck_wrt_alpha + h_alpha);
    
    H_theta = TN_theta*M*x_wrt_x - (Partial_ck_wrt_theta + h_theta);
    
    C = [Hz,Hx,H_alpha,H_theta];  % the matrix C
    
    % Calculate the new right hand vector
    ek = ck + h;
    
    % update gnew - used in LSE Problem.
    res_vec = rk;
    
    % Calculate the normalised residual of the solution.
    condition(ite) = norm(rk) ;
    
    % Update fnew - used in LSE Problem.
    p = -(yy-start_point);
    
end

fprintf('Required Number of iterations : %i \n',ite)


fx_new = fx_n + z_fx;
gx_new = gx_n + z_gx;
theta_new = theta(ite);
alpha_new = alpha(ite);


switch PLOT_GRAPHS
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
end

if ite == MAX_ITE_SNTLN
    % revers
    fx_new = fx
    gx_new = gx
    theta_new = theta(1)
    alpha_new = alpha(1)
    
end


end


function Hz = GetH_z(mincol,n,d,DY,DP,alpha)
if mincol<=(n-d+1)
    Hz=(DY-DP);
else
    Hz=DY-(alpha.*DP);
end
end

function Y = BuildY(m,n,t,x,opt_col,alpha,theta)
% Build the Matrix Y such that


% Insert 0 into the position of the optimal column in the x vector
xa = x(1:opt_col-1);
xb = x(opt_col:end);
x = [xa ;0 ;xb];

% Get the vectors u(\omega,\theta) and v(\omega,\theta)

% The matrix Y is such that S_{k}(f,g)*x = Y(x)*[f;g]

% Get the x values corresponding to v_{k}
ux = x(1:n-t+1);

% Get the x values corresponding to u_{k}
vx = x(n-t+2:end);


Y2 = BuildY1(vx,n,theta);
Y1 = BuildY1(ux,m,theta);

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

% Multiply the rows by theta^(i) where i is the row index
Y1 = diag(theta.^(0:1:m+n_t)) * Y1;
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
        zeros(num_zero_rows_bottom,m+1)];
    
    
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
