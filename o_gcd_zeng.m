function [u,v,w] = o_gcd_zeng(fx,gx)
% Given two polynomials, get the gcd using zeng method
% Note Zengs method takes strings as inputs, so must first convert to
% string from vector form

% Initialise the symbol x
x = sym('x');


sym_fx = GetSymbolicFromCoefficients(fx);
sym_gx = GetSymbolicFromCoefficients(gx);

string_fx = char(expand(sym_fx));
string_gx = char(expand(sym_gx));


[u,v,w] = PolynomialGCD(string_fx,string_gx);

end