function [] = o_gcd_zeng(ex_num,tol)
% Given two polynomials, get the gcd using Zeng PolynomialGCD method.

%Initialise symbolic variables.
x = sym('x');
y = sym('y');


switch ex_num
    case 1
        % Polynomials in x only
        sym_fx = (x-1)*(x-2) *(x+1);
        sym_gx = (x-1)*(x-2) *(x-5);
        
    case 2
        % Polynomials in y only
        sym_fx = (y-1)*(y-2) *(y+1); 
        sym_gx = (y-1)*(y-2) *(y-5);
        
    case 3
        % Polynomials in x and y where GCD(f,g) = f = g.
        sym_fx = (x-1)*(y-2);
        sym_gx = (x-1)*(y-2);
        
    case 4
        % Polynomials in x and y where GCD(f,g) = (x-1)(y-1)
        sym_fx = (x-1)*(y-1) * (x-3) * (x - 2);
        sym_gx = (x-1)*(y-1) * (x-5) * (x - 7);
end

% Get f(x) as string
string_fx = char(expand(sym_fx));

% Get g(x) as string
string_gx = char(expand(sym_gx));

display(string_fx)
display(string_gx)

[gcd,v,w,res,cond] = PolynomialGCD(string_fx,string_gx,tol);

display(gcd)

end