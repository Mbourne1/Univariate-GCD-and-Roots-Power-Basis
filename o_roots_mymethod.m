function [root_mult_array] = o_roots_mymethod(f)
% Given the polynomial f(x) calculate its real roots by square free 
% decomposition.
% 
% Inputs.
% 
% fx : Coefficients of polynomial f(x)
%
% Outputs.
% 
% root_mult_array : output two columns, the first containing the root, the
%                   second its corresponding multiplicity.


% Initialise an iteration counter
ite = 1;

% Initialise an array 'q' which stores the gcd outputs from each gcd
% calculation
fx{1} = f;

% let vGCD_Degree store the degrees corresponding to the array of
% GCDs stored in q.
M(1) = GetDegree(fx{1});

% Let theta_vec store all theta values used in each iteration.
vTheta(1) = 1;

% Get the number of distinct roots of f_{1}. Since this is unknown at this
% time, set number of distinct roots to be m_{1} = deg(f_{1}).
d(1) = GetDegree(fx);

% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.

while length(fx{ite})-1 > 0

    % if degree of f(ite_num) is greater than one
    if M(ite) > 1
        
        fprintf('GCD Calculation Loop iteration = %i \n', ite );
        fprintf('Compute GCD of f_{%i} and derivative f_{%i}\n\n',ite,ite);
        
        
        % Perform GCD computation
        
        % Get upper and lower bounds of the next GCD computation
        % M_{i+1} > M_{i} - d_{i-1}
        try
            lower_lim = M(ite)-d(ite-1);
            upper_lim = M(ite)-1;
            fprintf('Minimum degree of f_{%i}: %i \n', ite+1, M(ite)-d(ite-1));
            fprintf('Maximum degree of f_{%i}: %i \n', ite+1, M(ite)-1);
        catch
            lower_lim = 1;
            upper_lim = M(ite)-1;
        end

            
        % (fx_n,gx_n,dx, ux, vx, alpha, theta, t , lambda,mu)
        [fx{ite},~,fx{ite+1}, ux{ite} ,vx{ite},~,vTheta(ite+1),M(ite+1),~,~] ...
            = o_gcd_mymethod( fx{ite} , Differentiate(fx{ite}) , [lower_lim,upper_lim]);
        
        % Get number of distinct roots of f(ite)
        d(ite) = M(ite) - M(ite+1);
        
        
       
        fprintf('The computed deg(GCD(f_{%i},f_{%i}) is : %i \n',ite,ite,M(ite+1))
        
        fprintf('Number of distinct roots in f_{%i} : %i \n',ite,d(ite))
        
        fprintf('Degree of f_{%i} : %i \n',ite + 1, M(ite+1))
        

        % increment iteration number.
        ite = ite+1;
        
        
    elseif M(ite) == 1
        
        % if m=1, then n = 0, GCD has maximum degree 0.
        dx = 1;
        
        %theta_vec(ite_num+1) = 1;
        M(ite+1) = 0;
        
        fx{ite+1} = 1;
        ux{ite} = fx{ite}
        ite = ite+1;
        
        break;
        
        
    end
end


% Get the degree structure of the polynomials h_{i}
deg_struct_h = diff([M]);

% Get the degree structure of the polynomials w_{i}
deg_struct_w = diff([deg_struct_h 0]);

vMultiplicities = find(deg_struct_w~=0);


% Deconvolve the first set of polynomials q_{i}(x), which are the outputs
% of the series of GCD computations, to obtain the set of polynomaials h_{x}

% We can either obtain h(x) from the series of deconvolutions on f(x){i} or
% we can use the already calculated values ux{i}
method = 'By Deconvolution';
switch method
    case 'By Deconvolution'
        hx = Deconvolve(fx);
        
    
        hx_new = Deconvolve_Batch_Constrained(fx,vMultiplicities);
        
        arr_wx_new = Deconvolve(hx_new);
        display(arr_wx_new)
        
        % set the w1{max} = h1{max}
        arr_wx_new{length(arr_wx_new)+1} = hx{length(hx)};
        
    case 'From ux'
        hx = ux;
    otherwise
        error('err')
        
end

% for each w_{i}
for i = 1:1:length(arr_wx_new)
    wx_new = arr_wx_new{i}
    wx_new = wx_new./wx_new(2)

end


% Get the number of polynomials in h_{x}
[~,nCols_hx] = size(hx);

roots_wrt_x = [];


% If only one entry in h_{x}
% If only one entry in h_{x}
if nCols_hx == 1
    
    
    % if number of cols in h1 is only 1, then do not perform second
    % deconvolution, since only one entry in h1.
    % Note - this is a rare exception.
    vRoots = [];

    factor_x = hx{1};
    
    
    % Normalise the polynomial coefficients by the leading coefficient x^m
    factor_x = factor_x./factor_x(end);
    
    % Get the root from the factor ( x - r);
    rt = - factor_x(1);
    
    % get the roots with respect to y, and their multiplicities all set
    % to one.
    roots_wrt_x = [rt];
        
    % add the roots to the array of roots
    vRoots = [vRoots ; roots_wrt_x];
else
    % perform deconvolutions
    
    % Deconvolve the second set of polynomials
    wx = Deconvolve(hx);
    
    % w1 yields the simple, double, triple roots of input polynomial f.
    % w1{i} yields the roots of multiplicity i.
    
    % set the w1{max} = h1{max}
    wx{ite-1} = hx{ite-1};
    
    % get number of entries in w1
    [~,ncols_wx] = size(wx);
    
    % initialise an empty set
    vRoots = [];
    
    % for each multiplicity in w1.
    for i = 1:1:ncols_wx
        
        % Get the polynomial w_{i}
        poly_wi = wx{i};
        
        % Get the degree of w_{i}
        [nRows_wi,~] = size(wx{i});
        deg_wi = nRows_wi -1;
        
        % If w_{i} is of degree one, then is of the form (ax+b)
        % and has only one root. Add it to the list of roots.
        if deg_wi == 1;
                        
                 
            % Normalise the coefficients of w_{i} by dividing by the
            % leading coefficient. Since coefficients are orders in
            % asending powers, LC is the second coefficient.
            poly_wi = poly_wi./poly_wi(2);
            
            % Get the root
            rt = - poly_wi(1);
                      
            % Add the root to the [root, mult] matrix
            vRoots = [vRoots ; rt];
            
        elseif (deg_wi > 1)
            % The given multiplicity contains more than one root, such that
            % number of coefficients in greater than 2, use MATLAB roots
            % function to find roots.
                  
            % Get the roots in terms of w_{i} using Matlab Roots function.
            roots_wrt_x = roots(flipud(poly_wi));    
            
            % Add the computed roots to the array of roots.
            vRoots = [vRoots ; roots_wrt_x];
        end
    end
end

%% Get multiplicities of the roots
% Obtaining multiplicities of the calculated roots



% create a matrix where the first column contains the multiplicities, and
% the second column contains the number of roots of that multiplicity
nPolys_wi = length(deg_struct_w);
mat = [(1:1:nPolys_wi)' deg_struct_w'];

count = 1;
root_mult_array = [];
for i = 1:1:size(mat,1)
    % Get the number of roots of multiplicity i
    nRoots_of_Multi = mat(i,2);
    for j = 1:1:nRoots_of_Multi
        root_mult_array = [root_mult_array ; vRoots(count,1) i];
        count= count +1;
    end
end

%% Print the calculated roots and the corresponding multiplicities.
PrintRoots(root_mult_array,'MY METHOD');

end