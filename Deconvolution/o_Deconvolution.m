function [] = o_Deconvolution(ex_num, emin, bool_preproc)
% O_DECONVOLUTION Test the different methods of deconvolving the set of
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
% >> o_Deconvolution('1',1e-10,'y')

% Add path for examples
restoredefaultpath; savepath
addpath (...
    'Build Matrices',...
    'Deconvolution',...
    'Formatting',...
    'Preprocessing'...
    );

addpath(genpath('../Examples'));


% Set settings pertaining to this test
global SETTINGS
SETTINGS.PLOT_GRAPHS = true;
SETTINGS.MAX_ERROR_DECONVOLUTIONS = 1e-20;
SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS = 20;
SETTINGS.SEED = 1024;
SETTINGS.PREPROC_DECONVOLUTIONS = bool_preproc;



% % Get the factor array and multiplicity vector
[factor_mult_arr_f] = Deconvolution_Examples_Univariate(ex_num);

factor = factor_mult_arr_f(:,1);
vMult = double(factor_mult_arr_f(:,2));

% Get highest power of any factor
highest_pwr = max(vMult);

% %
% %
% Generate polynomials f_{0}(x) ,..., f_{m}(x) = 1. Where each f_{i+1}(x) is
% the f_{i+1} = GCD(f_{i},f'_{i}).
arr_sym_fx = cell(highest_pwr+1, 1);
vDeg_fx = zeros(highest_pwr+1, 1);

for i = 0:1:highest_pwr
    
    % Get the multiplicities of the roots of f_{i+1}
    mults = ((vMult - i) + abs(vMult-i)) ./2;
    
    % Get the symbolic polynomial f_{i+1}
    arr_sym_fx{i+1} = prod(factor.^(mults));
    
    % Get the degree of polynomial f_{i+1}(x)
    vDeg_fx(i+1) = double(feval(symengine, 'degree', (arr_sym_fx{i+1})));
end

% Display Polynomial f(x) in symbolic form
display(arr_sym_fx{1})

% Get the degree structure of the polynomials h_{i} where h_{i} =
% f_{i-1}(x)/f_{i}(x)
vDeg_arr_hx = diff(vDeg_fx);

% Get the degree structure of the polynomials w_{i} where w_{i} =
% h_{i-1}/h_{i}
vDeg_arr_wx = diff([vDeg_arr_hx; 0]);

% Get the multiplicities of the roots.
vMultiplicities = find(vDeg_arr_wx~=0);

% Get the sequence of polynomials h_{i}(x) in symbolic form
sym_arr_h = cell(length(arr_sym_fx)-1, 1);

for i = 1 : 1 : length(arr_sym_fx)-1
    sym_arr_h{i} = arr_sym_fx{i} / arr_sym_fx{i+1};
end

% %
% %
% Get coefficients vectors of f_{i}(x) and h_{i}(x)

% Get the number of polynomials in the array of f_{i}(x)
nPolys_arr_fx = size(arr_sym_fx, 1);

% Get the number of polynomials in the solution array h_{i}(x)
nPolys_arr_hx = size(arr_sym_fx, 1) - 1;

% Initialise the arrays f_{i}(x) and h_{i}(x)
arr_fx = cell(nPolys_arr_fx , 1);
arr_hx_exact = cell(nPolys_arr_hx , 1);

for i = 1:1:nPolys_arr_fx
    if i <= nPolys_arr_hx
        
        arr_fx{i,1} = sym2poly(arr_sym_fx{i})';
        arr_hx_exact{i,1} = sym2poly(sym_arr_h{i})';
        
    else
        
        arr_fx{i,1} = 1;
        
    end
    
end

% %
% %
% %
% Add noise to the coefficients of f_{i}(x)
arr_fx_noisy = cell(nPolys_arr_fx, 1);

for i = 1 : 1 : nPolys_arr_fx
    
    arr_fx_noisy{i,1} = AddNoiseToPoly(arr_fx{i},emin);
    
end


%--------------------------------------------------------------------------
% %
% %
% %
% Testing deconvolution separate
LineBreakLarge();
fprintf([mfilename 'Deconvolution Separate'  '\n']);

arr_hx_test_4 = Deconvolve_Separate(arr_fx_noisy);
vErrors_separate = GetErrors(arr_hx_test_4,arr_hx_exact);

%--------------------------------------------------------------------------
% %
% %
% %
% Testing standard deconvolution batch method
LineBreakLarge();
fprintf([mfilename ' : ' 'Deconvolve Batch' '\n']);

arr_hx_deconvolveBatch = Deconvolve_Batch(arr_fx_noisy);
vErrors_batch = GetErrors(arr_hx_deconvolveBatch,arr_hx_exact);

%--------------------------------------------------------------------------
% %
% %
% %
% Testing standard deconvolution batch method with STLN
LineBreakLarge();
fprintf([mfilename ' : ' 'Deconvolve batch with STLN' '\n'])

arr_hx_deconvolveBatchWithSTLN = Deconvolve_Batch_With_STLN(arr_fx_noisy);
vErrors_batch_with_STLN = GetErrors(arr_hx_deconvolveBatchWithSTLN,arr_hx_exact);


% -------------------------------------------------------------------------
% %
% %
% %
% Testing deconvolution batch method with constraints

LineBreakLarge()
fprintf([mfilename ' : ' 'Deconvolvve Batch Constrained' '\n'])

arr_hx_batchConstrained = Deconvolve_Batch_Constrained(arr_fx_noisy,vMultiplicities);
vErrors_batch_constrained = GetErrors(arr_hx_batchConstrained,arr_hx_exact);


% -------------------------------------------------------------------------
% %
% %
% %
% Testing deconvolution batch method with constraints and low rank approx

LineBreakLarge()
fprintf([mfilename ' : ' 'Deconvolve Batch Constrained With STLN \n'])

arr_hx_batchConstrainedWithSTLN = Deconvolve_Batch_Constrained_With_STLN(arr_fx_noisy,vMultiplicities);
vErrors_batch_constrained_with_STLN = GetErrors(arr_hx_batchConstrainedWithSTLN,arr_hx_exact);




%--------------------------------------------------------------------------

% %
% Plot Errors
if(SETTINGS.PLOT_GRAPHS)
    
        figure_name = sprintf([mfilename ' : ' 'Plotting errors of h(x)']);
        figure('name',figure_name)
        hold on
        plot(log10(vErrors_separate),'-.','DisplayName','Separate')
        plot(log10(vErrors_batch),'-*','DisplayName','Batch')
        plot(log10(vErrors_batch_with_STLN),'-s','DisplayName','Batch with STLN')
        plot(log10(vErrors_batch_constrained),'-s','DisplayName','Batch Constrained')
        plot(log10(vErrors_batch_constrained_with_STLN),'-s','DisplayName','Batch Constrained with STLN')
        xlabel('k');
        ylabel('log_{10} error h_{i}(x)');
        xlim([ 1 nPolys_arr_hx])
        legend(gca,'show');
        hold off
    
end
% ------------------------------------------------------------------------

Disp_error(vErrors_separate,'Separate');
Disp_error(vErrors_batch,'Batch');
Disp_error(vErrors_batch_with_STLN, 'Batch STLN');
Disp_error(vErrors_batch_constrained, 'Batch Constrained');
Disp_error(vErrors_batch_constrained_with_STLN, 'Batch Constrained STLN');




% %
% %
% %
% Print to file

A = ...
    [
    norm(vErrors_separate),...
    norm(vErrors_batch),...
    norm(vErrors_batch_with_STLN),...
    norm(vErrors_batch_constrained),...
    norm(vErrors_batch_constrained_with_STLN)
    ];


my_error.separate = norm(vErrors_separate);
my_error.batch = norm(vErrors_batch);
my_error.batchSTLN = norm(vErrors_batch_with_STLN);
my_error.batchConstrained = norm(vErrors_batch_constrained);
my_error.batchConstrainedSTLN = norm(vErrors_batch_constrained_with_STLN);

PrintToResultsFile(ex_num,bool_preproc,emin,my_error);


end

function [] = PrintToResultsFile(ex_num, bool_preproc, emin, my_error)
%
% % Inputs
%
% ex_num :Example Number
%
% bool_preproc


fullFileName = sprintf('Deconvolution/Results/Results_o_deconvolutions%s.txt',datetime('today'));

% If file already exists append a line
if exist(fullFileName, 'file')
    
    fileID = fopen(fullFileName,'a');
    WriteNewLine()
    fclose(fileID);
    
else % File doesnt exist so create it
    
    fileID = fopen( fullFileName, 'wt' );
    WriteHeader()
    WriteNewLine()
    fclose(fileID);
    
end

    function WriteNewLine()
        
        %
        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s \n',...
            datetime('now'),...
            ex_num,...
            bool_preproc,...
            emin,...
            my_error.separate,...
            my_error.batch,...
            my_error.batchSTLN ,...
            my_error.batchConstrained ,...
            my_error.batchConstrainedSTLN....
            );
        
    end

    function WriteHeader()
        
        fprintf(fileID,'DATE,EX_NUM,BOOL_PREPROC,NOISE,err_separate,err_batch,err_batch_STLN,err_constrained,err_constrained_STLN \n');
        
    end


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