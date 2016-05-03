function root_mult_arr = Examples_Roots(ex_num)
% Given an example number, return a polynomial in terms of its roots and
% their corresponding multiplicities.
%
% Input
% 
% Example Number
%
% Outputs.
%
% root_mult_arr : Array of roots and corresponding multiplicities in f(x)

switch ex_num
    case '1'
        root_mult_arr = ...
            [
            0.1 7;
            0.9 12;
            ];
    case '2'
        root_mult_arr = ...
            [
            0.1 1;
            0.2 2;
            0.9 3;
            ];
        
    case '3'   
    root_mult_arr = ...
            [
                6   10;
                -3  3;
                -7  1;
            ];
    case '4'
        root_mult_arr = ...
            [
            1   1
            2   2
            3   3
            4.5 1
            ];
    case '5'
        root_mult_arr = ...
            [
            1.05467 1
            2.24587 2
            5.54743 3
            1.75647 2
            ];
    case '6'
        root_mult_arr = ...
            [
            1.05467 1
            2.24587 2
            5.54743 3
            1.75647 2
            2.56478 5
            1.15445 1
            ];
    case 'Custom' 
        prompt = 'Enter the degree of Polynomial f(x) :';
        m = input(prompt);
        
        intvl_low = -1;
        intvl_high = +1;
        
        root_mult_arr = BuildRandomPolynomial(m,intvl_low,intvl_high);
        
    otherwise 
        error('err')
end




end