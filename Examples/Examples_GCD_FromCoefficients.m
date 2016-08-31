
function [fx,gx,dx] = Examples_GCD_FromCoefficients(ex_num)
switch ex_num
    case '1'
        
        syms x
        
        
        f = (x-2)^7 * (x-3)^12;
        g = diff(f,x);
        d = (x-2)^6 * (x-3)^11;
        
        
        % Display symbolic polynomials f and g
        display(f)
        display(g)
        
        % Get degree of f(x,y), g(x,y) and d(x,y)
        m = double(feval(symengine, 'degree', f));
        n = double(feval(symengine, 'degree', g));
        t = double(feval(symengine, 'degree', d));
        
        % Get coefficients of f(x,y), g(x,y) and d(x,y)
        fx = double(rot90(coeffs(f,[x],'All'),2));
        gx = double(rot90(coeffs(g,[x],'All'),2));
        dx = double(rot90(coeffs(d,[x],'All'),2));
        
        fx = fx';
        gx = gx';
        dx = dx';
        
    case '2'
        
        fx = [-1; -1; 1; 1];
        gx = [2; -1; -2; 1];
        dx = [-1; 0; 1];
        
    case '3'
        fx = [-5; 1; -5; 1];
        gx = [-2; 1; -2; 1];
        dx = [1; 0; 1];
        
    case '4'
        fx = [6.8; -17.6; 16.5; -6.7; 1];
        gx = [6;   -16  ; 15.5; -6.5; 1];
        
end



end
