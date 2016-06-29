function [dx] = GetGCD(ux,vx,fx,gx,t,alpha,theta)
% Get the coefficients of the GCD d(x), of f(x) and g(x), given u(x) and v(x).
% This function has two branches. d(x) can be computed either by utilising
% [u(x), v(x),f(x) and g(x)], or just [u(x) and f(x)]
%
%   Inputs.
%
%   ux      : Coefficients of input polynomial u(x)
%
%   vx      : Coefficients of input polynomial v(x)
%
%   fx      : Coefficients of polynomial f(x)
%
%   gx      : Coefficients of polynomial g(x)
%
%   t       : Degree of the GCD d(x)
%
%   alpha   : Optimal value of \alpha for Sylvester matrix S(f,\alpha g)
%
%   theta   : Optimal value of \theta s.t f(x) -> f(\theta,\omega)
%             and g(x) -> g(\theta,\omega)


GCD_METHOD = 'ux';

switch GCD_METHOD
    case 'ux'
        dx =  GetGCD_ux(ux,fx,t,theta);
    case 'ux and vx'
        dx = GetGCD_ux_and_vx(ux,vx,fx,gx,t,alpha,theta);
end


end

function [dx] = GetGCD_ux(ux,fx,t,theta)

% Get u(w) from u(x)
uw = GetWithThetas(ux,theta);

% Get f(w) from f(x)
fw = GetWithThetas(fx,theta);

% Build the toeplitz matrix C_{m}(u)
C_u = BuildT1(uw,t);

% Build the Right hand side vector [f(\omega) ; alpha.* g(\omega)]
rhs_vec = fw;

% Get the vector d.
d_ls = SolveAx_b(C_u,rhs_vec);

% Calculate d(\omega)
dw = d_ls;

% Get d(x)
dx = GetWithoutThetas(dw,theta);
end


function [dx] = GetGCD_ux_and_vx(ux,vx,fx,gx,t,alpha,theta)

% Get u(w) and v(w)
uw = GetWithThetas(ux,theta);
vw = GetWithThetas(vx,theta);

% Get f(w) and g(w)
fw = GetWithThetas(fx,theta);
gw = GetWithThetas(gx,theta);

% Get the GCD d(x) by the APF
C_u = BuildT1(uw,t);
C_v = BuildT1(vw,t);

C1 = ...
    [
    C_u;
    C_v;
    ];

% Build the Right hand side vector [f(\omega) ; alpha.* g(\omega)]
rhs_vec = [fw;alpha.*gw];

% Get the vector d.
d_ls = SolveAx_b(C1,rhs_vec);

% Calculate d(\omega)
dw = d_ls;

% Get d(x)
dx = GetWithoutThetas(dw,theta);
end