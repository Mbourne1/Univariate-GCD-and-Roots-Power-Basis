function [GM_fx, GM_gx, GM_hx, alpha, theta] = Preprocess_3Polys(fx, gx, hx)
% Preprocess(fx,gx)
%
% Preprocess the input polynomials to obtain geometric means of the
% coefficients of f(x) and g(x), and optimal values of alpha and theta such
% that the ratio of entries in sylvester matrix S(f,g) is minimized.
%
% Inputs.
%
% fx : Vector of coefficients of polynomial f(x)
%
% gx : Vector of coefficients of polynomial g(x)
%
% hx : Vector of coefficients of polynomial h(x)
%
% % Outputs
%
% [GM_fx, GM_gx, GM_hx] : Geometric mean of entries of f(x), g(x) and h(x)
%
% alpha :
%
% theta :

global SETTINGS

% Get the degree of polynomials f(x), g(x) and h(x)
m = GetDegree(fx);
n = GetDegree(gx);
o = GetDegree(hx);


% Get Mean of entries of f(x) in C_{n-k}(f)
GM_fx = GetMean(fx);

% Get Mean of entries of g(x) in C_{m-k}(g)
GM_gx = GetMean(gx);

% Get Mean of entries of g(x) in C_{m-k}(g)
GM_hx = GetMean(hx);


% Divide f(x) and g(x) by geometric mean
fx_n = fx./ GM_fx;
gx_n = gx./ GM_gx;
hx_n = hx./ GM_hx;

switch SETTINGS.BOOL_ALPHA_THETA
    case 'y'
        
        % Get opitmal values of alpha and theta
        
        [alpha1, theta1] = GetOptimalAlphaAndTheta_3Polys(fx_n, gx_n, hx_n);
        display(alpha1)
        display(theta1)
        
%         alpha2 = GetOptimalAlpha_3Polys(fx_n, gx_n, hx_n);
%         display(alpha2)
%         
%         theta2 = GetOptimalTheta_3Polys(fx_n, alpha2.*gx_n, alpha2.*hx_n);
%         display(theta2)
        
        alpha = alpha1;
        theta = theta1;
        
        
        % Obtain f(w) and g(w) from f(x) and g(x)]
        fw = GetWithThetas(fx_n, theta);
        gw = GetWithThetas(gx_n, theta);
        hw = GetWithThetas(hx_n, theta);
        
        F_max = max(fx_n);
        F_min = min(fx_n);
        G_max = max(gx_n);
        G_min = min(gx_n);
        H_max = max(hx_n);
        H_min = min(hx_n);
        
        PrintToFile(F_max,F_min,G_max,G_min,m,n,'1','1');
        
        %%
        F_max = max(fw);
        F_min = min(fw);
        G_max = max(gw);
        G_min = min(gw);
        H_max = max(hw);
        H_min = min(hw);
        
        PrintToFile(F_max,F_min,G_max,G_min,m,n,alpha,theta);
        
        % Plot the unprocessed and preprocessed coefficients of
        % f(x), f(w), g(x) and g(w).
        
%         switch SETTINGS.PLOT_GRAPHS
%             case 'y'
%                 PlotCoefficients(fx,fw,'f');
%                 PlotCoefficients(gx,gw,'g');
%                 PlotCoefficients(hx,hw,'h');
%             case 'n'
%             otherwise
%                 error('err')
%         end
    case 'n'
        
        % Dont apply preprocessing
        theta = 1;
        alpha = 1;
        
    otherwise
        error('bool_preproc either y or n');
end

end



function [] = PrintToFile(F_max,F_min,G_max,G_min,m,n,alpha,theta)

global SETTINGS

fullFileName = 'Preprocessing/Results_Preprocessing.txt';


if exist(fullFileName, 'file')
    fileID = fopen(fullFileName,'a');
    fprintf(fileID,'%s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s, \t %s \n',...
        datetime('now'),...
        SETTINGS.PROBLEM_TYPE,...
        SETTINGS.EX_NUM,...
        int2str(m),...
        int2str(n),...
        SETTINGS.MEAN_METHOD,...
        F_max,...
        F_min,...
        G_max,...
        G_min,...
        SETTINGS.BOOL_ALPHA_THETA,...
        alpha,...
        theta...
        );
    fclose(fileID);
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
  uiwait(msgbox(warningMessage));
end

end