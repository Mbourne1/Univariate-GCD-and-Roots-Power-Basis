function [] = o_roots_mymethod(fx)
% Given the polynomial f(x) calculate its real roots by square free 
% decomposition.

% Get derivative of f(x)
gx = Differentiate(fx);

% Initialise an iteration counter
ite_num = 1;

% Print current iteration number
fprintf('GCD Calculation Loop iteration = %i \n\n', ite_num);

% Perform an intial GCD Calculation, from this we obtain new forms of f and
% g and a calculated 'd'.
[fx,~,dx, ~ ,~,~,theta] = o1(fx,gx);

% Initialise an array 'q' which stores the gcd outputs from each gcd
% calculation
q{1} = fx;
q{2} = dx;

% let degree_vector store the degrees corresponding to the array of
% GCDs stored in q.
degree_vec(1) = length(fx)-1;
degree_vec(2) = length(dx)-1;

% Let theta_vec store all theta values used in each iteration.
theta_vec(1) = 1;
theta_vec(2) = theta;

% increment the iteration number
ite_num = ite_num+1;


% Whilst the most recently calculated GCD has a degree greater than
% zero. ie is not a constant, perform a gcd calculation on it and its
% derivative.
while length(q{ite_num})-1 > 0
    
    % set polynomial f to be the most recently calculated GCD.
    fx = q{ite_num};
    
    % set polynomial g to be the derivative of the GCD.
    gx = Differentiate(q{ite_num});
    
    % get degrees m and n of polynomials f and g respectively.
    m = size(fx,1) - 1;
    
    % if degree of f is greater than one
    if m > 1
        
        fprintf('GCD Calculation Loop iteration = %i \n\n',ite_num );
        [fx,~,dx, ~ ,~,~,theta] = o1(fx,gx);
        
        % add the value of theta used in this GCD calculation to the theta
        % vector
        theta_vec(ite_num+1) = theta;
        
        % add the degree of the calculated GCD to the degree vector
        degree_vec(ite_num+1) = length(dx)-1;
        
        % replace input f with updated f
        q{ite_num} = fx;
        
        % add vector dx to the array of gcds 'q'
        q{ite_num+1} = dx;
        
        % increment iteration number.
        ite_num = ite_num+1;
        
        
    elseif m == 1
        
        % if m=1, then n = 0, GCD has maximum degree 0.
        fprintf('Only one subresultant exists \n')
        dx = 1;
        
        %theta_vec(ite_num+1) = 1;
        degree_vec(ite_num+1) = length(dx)-1;
        q{ite_num+1} = dx;
        ite_num = ite_num+1;
        break;
        
        
    end
end


 fprintf('All Polynomials from GCD Calculations \n.')
 for i = 1:1:length(q)
     fprintf('\n')
     disp(q{i})
 end

% Deconvolve the first set of polynomials.
h1 = Deconvolve(q);

[~,c] = size(h1);

if c == 1
    % if number of cols in h1 is only 1, then do not perform second
    % deconvolution, since only one entry in h1.
    % Note - this is a rare exception.
    root_arr = [];

    ax = h1{1};
    
    % Normalise the polynomial coefficients
    ax = ax./ax(1);
       
    
    % get the roots in terms of z^{i} where z^{i} =(\frac{y}{(1-y)})^{i}.
    ax_roots = roots(flipud(ax));
        
    % Initialise a vector of ones
    one = ones(length(ax_roots),1);
    
    % get the roots with respect to y, and their multiplicities all set
    % to one.
    roots_wrt_z = [ax_roots one];
        
    % add the roots to the array of roots
    root_arr = [root_arr ; roots_wrt_z];
else
    
    % perform deconvolutions
    
    % Deconvolve the second set of polynomials
    w1 = Deconvolve(h1);
    
    
    % w1 yields the simple, double, triple roots of input polynomial f.
    % w1{i} yields the roots of multiplicity i.
    
    % set the w1{max} = h1{max}
    w1{ite_num-1} = h1{ite_num-1};
    
    % get number of entries in w1
    [~,c] = size(w1);
    
    % initialise an empty set
    root_arr = [];
    
    % for each multiplicity in w1.
    for i = 1:1:c
        
        % if the polynomial of said multiplicity is of length 2, degree 1, then
        % only one root exists for this multiplicity. Add it to the list wp1.
        if (length(w1{i}) == 2)
            
            
            % get the polynomial, whose roots have multiplicty i, in bernstein
            % form, where coefficients are in terms of (1-y)^{m-i}y^{i}.
            ax = w1{i};
            
            % Normalise the polynomial coefficients
            ax = ax./ax(1);
            
            % Convert to power form, so that coefficients are in terms of y^{i}
            % rather than (1-y)^{m-i}y^{i}.
            a_rt = [ax(1,:) ; ax(2,:)];
            
            % add the root to the [root, mult] matrix
            root_arr = [root_arr ; a_rt];
            
        elseif (length(w1{i}) > 2)
            % The given multiplicity contains more than one root, such that
            % number of coefficients in greater than 2, use MATLAB roots
            % function to find roots.
            % display('Multiplicity contains more than one root. ')
            
            % get the polynomial a(w), whose roots have multiplicity i, in bernstein
            % form.
            ax = w1{i};
            
            % Normalise the polynomial coefficients
            ax = ax./ax(1);
            
            % get the degree of a(w)
            n = length(w1{i}) - 1;
                       
            % get the roots in terms of z^{i} where z^{i} =(\frac{y}{(1-y)})^{i}.
            roots_wrt_z = roots(flipud(ax));    
            
            % Initialise a vector of ones
            one = ones(length(roots_wrt_z),1);
            
            % get the roots with respect to y, and their multiplicities all set
            % to one.
            roots_wrt_z = [roots_wrt_z one];
            
            % add the roots to the array of roots
            root_arr = [root_arr ; roots_wrt_z];
        end
    end
end


% Obtaining multiplicities of the calculated roots
roots_multiplicty = [];
while sum(degree_vec) ~= 0
    % Get index of first zero.;
    index = min(find(degree_vec==0)) - 1;
    minus_vector = zeros(1,length(degree_vec));
    minus_vector(1,1:index) = index:-1:1;
    degree_vec = degree_vec - minus_vector;
    roots_multiplicty = [roots_multiplicty ; index];
end



roots_calc = [root_arr(:,1) flipud(roots_multiplicty)];

%% Print the calculated roots and the corresponding multiplicities.
PrintRoots(roots_calc,'MY METHOD');

end