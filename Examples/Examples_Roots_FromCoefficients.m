
function fx_exact = Examples_Roots_FromCoefficients(ex_num)

switch ex_num
    
    case '1'
        
        
    
%     case '1'
%         fx_exact = [1; 1; -2; -2; 1; 1];
%     case '2'
%         fx_exact = [-1; -1; 1; 1];
%     case '3'
%         fx_exact = [-5; 1; -5; 1];
%     case '4'
%         fx_exact = [48; -32; 0; 0; 1];
%     case '5' %2 - 2*x - 6*x^4 + 6*x^5 + 6*x^8 - 6*x^9 -2*x^12 + 2*x^13
%         fx_exact = [2; -2; 0; 0; -6; 6; 0; 0; 6; -6; 0; 0; -2; 2];
%     case '6'
%         % intersection of 
%         % x^2 + (y+5)^2 - 5^2
%         % 'x^2-y'
%         fx_exact = [0; 0; 11; 0; 1];
    otherwise 
        error('err')
end

end