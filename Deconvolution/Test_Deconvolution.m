function [] = Test_Deconvolution(ex_num,emin,bool_preproc)
% TEST_DECONVOLUTION Test the different methods of deconvolving the set of
% polynomials f_{i}(x), to form the set of polynomials h_{i} 
% where h_{i} = f{i}/f_{i+1}
%
% % Inputs
%
% ex_num : Example number (String)
%
% emin : Lower noise level
%
% emax : Upper noise level
%
% % Outputs.
%
%
% Example
%
% >> Test_Deconvolution('1',1e-10,1e-8)

% Set settings pertaining to this test
global SETTINGS
SETTINGS.PLOT_GRAPHS = 'n';
SETTINGS.MAX_ERROR_DECONVOLUTIONS = 1e-13;
SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS = 100;
SETTINGS.SEED = 1024;
SETTINGS.PREPROC_DECONVOLUTIONS = bool_preproc;


% Input f_{i} polynomials
x = sym('x');

switch ex_num
    case '1'
        
        % Create set of factors whose multiplicities are defined in vMult
        factor(1,1) = (x+1.0177);
        factor(2,1) = (x-0.5296);
        
        % Set multiplicity of each factor
        vMult = [9 ; 14];
        
    case '2'
                                
        % Create set of factors whose multiplicities are defined in vMult
        factor(1,1) = (x-2);
        factor(2,1) = (x-3.2789);
        factor(3,1) = (x-1.589);
        factor(4,1) = (x-0.7213);
        factor(5,1) = (x-1.5432);
        factor(6,1) = (x+5.72);
        
        % Set multiplicitiy of each factor
        vMult = [ 1; 3; 4; 4; 5; 12 ];
        
    case '3'
        factor(1,1) = (x-2);
        factor(2,1) = (x-3.2789);
        factor(3,1) = (x-1.589);
        vMult = [2; 4; 12 ];
        
    case '4'
        factor(1,1) = (x-0.56897);
        factor(2,1) = (x+1.24672);
        factor(3,1) = (x+0.56921);
        vMult = [3; 6; 9];
      
end

% Get highest power of any factor
highest_pwr = max(vMult);

% %
% %
% Generate polynomials f_{0}(x) ,..., f_{m}(x) = 1. Where each f_{i+1}(x) is
% the f_{i+1} = GCD(f_{i},f'_{i}).
arr_sym_f = cell(highest_pwr+1,1);
vDeg_fx = zeros(highest_pwr+1,1);

for i = 0:1:highest_pwr
    
    % Get the multiplicities of the roots of f_{i+1}
    mults = ((vMult - i) + abs(vMult-i)) ./2;
    
    % Get the symbolic polynomial f_{i+1}
    arr_sym_f{i+1} = prod(factor.^(mults));
    
    % Get the degree of polynomial f_{i+1}(x)
    vDeg_fx(i+1) = double(feval(symengine, 'degree', (arr_sym_f{i+1})));
end

% Display Polynomial f(x) in symbolic form 
display(arr_sym_f{1})

% Get the degree structure of the polynomials h_{i} where h_{i} =
% f_{i-1}(x)/f_{i}(x)
vDeg_arr_hx = diff(vDeg_fx);

% Get the degree structure of the polynomials w_{i} where w_{i} =
% h_{i-1}/h_{i}
vDeg_arr_wx = diff([vDeg_arr_hx; 0]);

% Get the multiplicities of the roots.
vMultiplicities = find(vDeg_arr_wx~=0);

% Get the sequence of polynomials h_{i}(x) in symbolic form
sym_arr_h = cell(length(arr_sym_f)-1,1);
for i = 1:1:length(arr_sym_f)-1
    sym_arr_h{i} = arr_sym_f{i} / arr_sym_f{i+1};
end

% %
% %
% Get coefficients vectors of f_{i}(x) and h_{i}(x)
nPolys_arr_fx = size(arr_sym_f,1);
nPolys_arr_hx = size(arr_sym_f,1) - 1;

arr_fx = cell(nPolys_arr_fx,1);
arr_hx = cell(nPolys_arr_hx,1);

for i = 1:1:nPolys_arr_fx
    if i <= nPolys_arr_hx
        arr_fx{i,1} = sym2poly(arr_sym_f{i})';
        arr_hx{i,1} = sym2poly(sym_arr_h{i})';
    else
        arr_fx{i,1} = 1;
    end
    
end

% %
% %
% %
% Add noise to the coefficients of f_{i}(x)
arr_fx_noisy = cell(nPolys_arr_fx,1);

for i = 1:1:nPolys_arr_fx
    arr_fx_noisy{i,1} = Noise(arr_fx{i},emin);
end


%--------------------------------------------------------------------------
% %
% %
% %
% Testing deconvolution separate
LineBreakLarge();
fprintf([mfilename 'Deconvolution Separate'  '\n']);

arr_hx_test_4 = Deconvolve_Separate(arr_fx_noisy);
err_deconvolveSeparate = GetErrors(arr_hx_test_4,arr_hx);

%--------------------------------------------------------------------------
% %
% %
% %
% Testing standard deconvolution batch method
LineBreakLarge();
fprintf([mfilename ' : ' 'Deconvolve Batch' '\n']);

arr_hx_deconvolveBatch = Deconvolve_Batch(arr_fx_noisy);
err_deconvolveBatch = GetErrors(arr_hx_deconvolveBatch,arr_hx);

%--------------------------------------------------------------------------
% %
% %
% %
% Testing standard deconvolution batch method with STLN
LineBreakLarge();
fprintf([mfilename ' : ' 'Deconvolve batch with STLN' '\n'])

arr_hx_deconvolveBatchWithSTLN = Deconvolve_Batch_With_STLN(arr_fx_noisy);
err_deconvolveBatchWithSTLN = GetErrors(arr_hx_deconvolveBatchWithSTLN,arr_hx);


% -------------------------------------------------------------------------
% %
% %
% %
% Testing deconvolution batch method with constraints

LineBreakLarge()
fprintf([mfilename ' : ' 'Deconvolvve Batch Constrained' '\n'])

arr_hx_batchConstrained = Deconvolve_Batch_Constrained(arr_fx_noisy,vMultiplicities);
err_batchConstrained = GetErrors(arr_hx_batchConstrained,arr_hx);


% -------------------------------------------------------------------------
% %
% %
% %
% Testing deconvolution batch method with constraints and low rank approx

LineBreakLarge()
fprintf([mfilename ' : ' 'Deconvolve Batch Constrained With STLN'])

arr_hx_batchConstrainedWithSTLN = Deconvolve_Batch_Constrained_With_STLN(arr_fx_noisy,vMultiplicities);
err_batchConstrainedWithSTLN = GetErrors(arr_hx_batchConstrainedWithSTLN,arr_hx);




%--------------------------------------------------------------------------

% %
% Plot Errors

figure_name = sprintf([mfilename ' : ' 'Plotting errors of h(x)']);
figure('name',figure_name)
hold on
plot(log10(err_deconvolveSeparate),'-.','DisplayName','Separate')
plot(log10(err_deconvolveBatch),'-*','DisplayName','Batch')
plot(log10(err_deconvolveBatchWithSTLN),'-s','DisplayName','Batch with STLN')
plot(log10(err_batchConstrained),'-s','DisplayName','Batch Constrained')
plot(log10(err_batchConstrainedWithSTLN),'-s','DisplayName','Batch Constrained with STLN')
xlabel('k');
ylabel('log_{10} error h_{i}(x)');
xlim([ 1 nPolys_arr_hx])
legend(gca,'show');
hold off

% ------------------------------------------------------------------------

Disp_error(err_deconvolveSeparate,'Separate');
Disp_error(err_deconvolveBatch,'Batch');
Disp_error(err_deconvolveBatchWithSTLN, 'Batch STLN');
Disp_error(err_batchConstrained, 'Batch Constrained');
Disp_error(err_batchConstrainedWithSTLN, 'Batch Constrained STLN');




% %
% %
% %
% Print to file

A = ...
    [
    norm(err_deconvolveSeparate),...
    norm(err_deconvolveBatch),...
    norm(err_deconvolveBatchWithSTLN),...
    norm(err_batchConstrained),...
    norm(err_batchConstrainedWithSTLN)
    ];

fileID = fopen('Deconvolution/Test_Deconvolution.txt','a');
fprintf(fileID,'%s %s %6.2e %6.2e %6.2e %6.2e %6.2e %6.2e  \n',ex_num,bool_preproc,emin,A);
fclose(fileID);

end

function Disp_error(v_err,name)

err = norm(v_err);

fprintf([sprintf('Error in %s method : %e',name,err) '\n']);

end

function [vErrors] =  GetErrors(arr_hx_comp , arr_hx_exact)
% GETERRORS Get the distance between each h

% Initialise vector to store errors.
vErrors = zeros(size(arr_hx_comp,1),1);

% Compare each computed h{i} with actual h_{i}
for i = 1:1:size(arr_hx_comp,1)
    
    exact = arr_hx_exact{i}./ arr_hx_exact{i,1}(1,1);
    comp = arr_hx_comp{i}./arr_hx_comp{i,1}(1,1);
    
    vErrors(i) = norm(exact - comp) ./ norm(exact);
    
    
end
end