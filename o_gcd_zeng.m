function [u,v,w] = o_gcd_zeng(fx,gx)
% Given two polynomials, get the gcd using zeng method
% Note Zengs method takes strings as inputs, so must first convert to
% string from vector form




% Get symbolic expression for f(x)
sym_fx = GetSymbolicFromCoefficients(fx);

% Get symbolic expression for g(x)
sym_gx = GetSymbolicFromCoefficients(gx);


%sample
x = sym('x')
y = sym('y')
sym_fx = (y-1)*(y-2)*(x-5)
sym_gx = (y-1)*(y-2)*(x+1)


% Get f(x) as string
string_fx = char(expand(sym_fx));

% Get g(x) as string
string_gx = char(expand(sym_gx));

%string_fx = '1*y^2*x - 3*x*y +2*x -5*y^2 +15*y -10'
%string_gx = '1*y^2*x - 3*x*y +2*x + 1*y^2 -3*y +2'

%string_fx = '-45*x*y - 15*x^3*y - 20*x*y^2 + 27*x*y^3 + 9*x^3*y^3 + 12*x*y^4';
%string_gx = '45*x^2*y^2 + 15*x^2*y^3 - 27*x^2*y^4 - 9*x^2*y^5';

%f_mat = [-10 15 -5 ; 2 -3 1];
%g_mat = [2 -3 1 ; 2 -3 1];
 
tol = 1e-10;
[u,v,w,res,cond] = PolynomialGCD(string_fx,string_gx);

%[u2,v2,w2] = mvGCD(f_mat,g_mat)

end