function Sk = BuildT_3Polys(fx, gx, hx, k)
% Build the Sylvester Subresultant matrix S_{k}(f,g)
%
% % Inputs
%
% fx : Coefficients of the polynomial f(x)
%
% gx : Coefficients of the polynomial g(x)
%
% hx : Coefficients of the polynomial h(x)
%
% k : Index of Sylvester Subresultant matrix to be constructed.


% Get degree of polynomial f(x)
m = GetDegree(fx);

% Get degree of polynomial g(x)
n = GetDegree(gx);

% Get the degree of polynomial h(x)
o = GetDegree(hx);

C_f1 = BuildT1(fx,n-k);
C_f2 = BuildT1(fx,o-k);

C_g = BuildT1(gx,m-k);
C_h = BuildT1(hx,m-k);

block = blkdiag(C_f1,C_f2);
col = [C_g ; C_h];

Sk = [block col];


end