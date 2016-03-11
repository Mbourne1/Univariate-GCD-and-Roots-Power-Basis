function [roots_fx] = Examples_Univariate_Implicit(ex_num)
% Given the example number, return the set of coefficients for polynomials
% f(x) and g(x).

% roots_fx : A matrix where each row consists of [root multiplicity]
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
        roots_fx = ...
            [
            1   1;
            2   2;
            1.5 1;
            ];
        
end







end