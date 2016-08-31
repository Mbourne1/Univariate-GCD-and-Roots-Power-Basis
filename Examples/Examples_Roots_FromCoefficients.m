
function fx_exact = Examples_Roots_FromCoefficients(ex_num)

x = sym('x');

switch ex_num
    
    case '1'
        
        f = (x - 0.5)^4 * (x + 0.75)^7;
    
    case '2'
        
        f = (x - 1.5)^4 * (x + 0.75)^7 * (x - 10.1)^3;
        
    case '3'
        
        f = (x - 1.5)^4 * (x + 0.75)^7 * (x - 10.1)^3 * (x-0.17523547)^5;
        
    otherwise 
        error('err')
end

display(f);

m = double(feval(symengine, 'degree', f));

fx_exact = flipud(double(coeffs(f,x,'All'))');

end
