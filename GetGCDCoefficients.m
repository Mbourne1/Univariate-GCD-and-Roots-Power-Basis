function [dx] = GetGCDCoefficients(ux,vx,fx,gx,t)
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


GCD_METHOD = 'ux and vx';

switch GCD_METHOD
    case 'ux'
        dx =  GetGCD_ux(ux,fx,t);
    case 'ux and vx'
        dx = GetGCD_ux_and_vx(ux,vx,fx,gx,t);
end


end

function [dx] = GetGCD_ux(ux,fx,t)

% Build the toeplitz matrix C_{m}(u)
C_u = BuildT1(ux,t);

% Build the Right hand side vector [f(\omega) ; alpha.* g(\omega)]
rhs_vec = fx;

% Get the vector d.
d_ls = SolveAx_b(C_u,rhs_vec);

% Calculate d(\omega)
dx = d_ls;

end


function [dx] = GetGCD_ux_and_vx(ux,vx,fx,gx,t)


% Get the GCD d(x) by the APF
C_u = BuildT1(ux,t);
C_v = BuildT1(vx,t);

C1 = ...
    [
    C_u;
    C_v;
    ];

% Build the Right hand side vector [f(\omega) ; alpha.* g(\omega)]
rhs_vec = [fx;gx];

% Get the vector d.
d_ls = SolveAx_b(C1,rhs_vec);

% Calculate d(\omega)
dx = d_ls;


end