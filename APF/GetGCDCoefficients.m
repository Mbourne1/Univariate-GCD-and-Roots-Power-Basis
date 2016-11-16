function dx = GetGCDCoefficients(ux,vx,fx,gx,k)
% Get the common divisor of degree k, of the polynomials f(x) and g(x).
%
% % Inputs
%
% ux : Coefficients of the polynomial u(x)
%
% vx : Coefficients of the polynomial v(x)
%
% fx : Coefficients of the polynomial f(x)
%
% gx : Coefficients of the polynomial g(x)
%
% % Outputs
%
% dx : Coefficients of the polynomial d(x) of degree k.


GCD_METHOD = 'ux and vx';

switch GCD_METHOD
    case 'ux'
        
        dx =  GetGCD_ux(ux,fx,k);
        
    case 'ux and vx'
        
        dx = GetGCD_ux_and_vx(ux,vx,fx,gx,k);
        
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


function [dx] = GetGCD_ux_and_vx(ux, vx, fx, gx, k)
%
% % Inputs
%
% ux : Coefficients of polynomial u(x)
%
% vx : Coefficients of polynomial v(x)
%
% fx : Coefficients of polynomial f(x)
%
% gx : Coefficients of polynomial g(x)
% 
% k : Degree of common divisor d(x)

% Get the GCD d(x) by the APF
C_u = BuildT1(ux,k);
C_v = BuildT1(vx,k);

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