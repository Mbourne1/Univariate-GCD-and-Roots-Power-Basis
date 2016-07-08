
function [fx,gx,dx_exact] = GCD_Examples_FromCoefficients(ex_num)
switch ex_num
    case '1'
        
        fx = [1; 1; -2; -2; 1; 1];
        gx = [-2; 1; 4; -2; -2; 1];
        dx_exact = [1; 0; -2; 0; 1];
    case '2'
        
        fx = [-1; -1; 1; 1];
        gx = [2; -1; -2; 1];
        dx_exact = [-1; 0; 1];
        
    case '3'
        fx = [-5; 1; -5; 1];
        gx = [-2; 1; -2; 1];
        dx_exact = [1; 0; 1];
        
    case '4'
        fx = [6.8; -17.6; 16.5; -6.7; 1];
        gx = [6;   -16  ; 15.5; -6.5; 1];
        
end

end
