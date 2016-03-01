function [fx] = Examples_Intersection_Implicit(ex_num)
% Given the example number, return the set of coefficients for polynomials
% f(x) and g(x).

% roots_fx : A matrix where each row consists of [root multiplicity]
% pairs.

% roots_gx : A matrix where each row consists of [root multiplicity]
% pairs.

switch ex_num
    case '1'
        
        roots_fx = ...
            [
            1   1;
            2   2;
            1.7 1;
            ];
    case '2'
        roots_gx = ...
            [
            1   1;
            2   2;
            1.5 1;
            ];
        
end

% Get the polynomial coefficients of f(x) and g(x) as a column vector. Note
% that the first entry is the column is a_{0}x^0 and the last is a_{m}x^{m}
fx = GetCoefficients(roots_fx);





end