function [] = Test_Deconvolution
% Test the different methods of deconvolving the polynomials f_{i}(x), to
% form the set of polynomials h_{i} where h_{i} = f{i}/f_{i+1}
%
%
%
%
%


% Input f_{i} polynomials
x = sym('x');

ex_num = '2';

switch ex_num
    case '1'
        
        f{1} = (x-2)^7 * (x-3)^12;
        f{2} = (x-2)^6 * (x-3)^11;
        f{3} = (x-2)^5 * (x-3)^10;
        f{4} = (x-2)^4 * (x-3)^9;
        f{5} = (x-2)^3 * (x-3)^8;
        f{6} = (x-2)^2 * (x-3)^7;
        f{7} = (x-2)^1 * (x-3)^6;
        f{8} = (x-3)^5;
        f{9} = (x-3)^4;
        f{10} = (x-3)^3;
        f{11} = (x-3)^2;
        f{12} = (x-3);
        f{13} = 1;
        
        vMult = [7 , 12];
        
    case '2'
        
        
        vMult = [ 2 3 3 3 4 8 ];
        
        factor(1) = (x-2);
        factor(2) = (x-3.2789);
        factor(3) = (x-1.589);
        factor(4) = (x-0.7213);
        factor(5) = (x-1.5432);
        factor(6) = (x+5.72);
        
        highest_pwr = max(vMult);
        
        
        for i = 0:1:highest_pwr
            
            mults = ((vMult - i) + abs(vMult-i)) ./2;  
            
            f{i+1} = prod(factor.^(mults));
            M(i+1) = double(feval(symengine, 'degree', (f{i+1})));
        end
        
        
end

% Get the degree structure of the polynomials h_{i}
deg_struct_h = diff(M);
% Get the degree structure of the polynomials w_{i}
deg_struct_w = diff([deg_struct_h 0]);
% Get the multiplicities of the roots.
vMultiplicities = find(deg_struct_w~=0);


for i = 1:1:length(f)-1
    h{i} = f{i} / f{i+1};
end

% %
% %
% Get coefficients vectors
for i = 1:1:length(f)
    try
        arr_fx{i,1} = sym2poly(f{i})';
        arr_hx{i,1} = sym2poly(h{i})';
    catch
        arr_fx{i,1} = 1;
    end
    
end


% %
% %
% %
% Testing deconvolution batch method with constraints

LineBreakLarge()
fprintf('Deconvolve with the constraint that all h_{i} for i between m_{j} and m_{j+1} are equal \n')
fprintf('Staggered Staircase Method \n\n')

arr_hx_test_1 = Deconvolve_Batch_Constrained(arr_fx,vMultiplicities);

for i = 1:1:length(arr_hx_test_1)
    %display(arr_hx_test_1{i} );
end

% Compare each computed h{i} with actual h_{i}
for i = 1:1:length(arr_hx_test_1)
    
    exact = arr_hx{i};
    comp = arr_hx_test_1{i};
    
    err_measure = norm(exact - comp) ./ norm(exact);
    fprintf([mfilename ': ' sprintf('%i : Error : %2.4e \n', i,err_measure)]);
end


% %
% %
% %
% Testing standard deconvolution batch method
LineBreakLarge();
fprintf('Deconvolve Staircase method \n')
fprintf('Deconvolution Batch \n')
arr_hx_test_2 = Deconvolve_Batch(arr_fx);
for i = 1:1:length(arr_hx_test_2)
    %display(arr_hx_test_2{i});
end

% Compare each computed h{i} with actual h_{i}
for i = 1:1:length(arr_hx_test_2)
    
    exact = arr_hx{i};
    comp = arr_hx_test_2{i};
    
    err_measure = norm(exact - comp) ./ norm(exact);
    fprintf([mfilename ': ' sprintf('%i : Error : %2.4e \n', i,err_measure)]);
end




% %
% %
% %
% Testing deconvolution
LineBreakLarge();
fprintf('Deconvolve - Each deconvolution is independent \n');
fprintf('Deconvolution Batch Separate \n');
arr_hx_test_3 = Deconvolve_Batch_Separate(arr_fx);
for i = 1:1:length(arr_hx_test_3)
    %display(arr_hx_test_3{i});
end
% Compare each computed h{i} with actual h_{i}
for i = 1:1:length(arr_hx_test_3)
    
    exact = arr_hx{i};
    comp = arr_hx_test_3{i};
    
    err_measure = norm(exact - comp) ./ norm(exact);
    fprintf([mfilename ': ' sprintf('%i : Error : %2.4e \n', i,err_measure)]);
end


end