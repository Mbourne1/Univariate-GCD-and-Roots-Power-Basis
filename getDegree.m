function t = getDegree(fw,gw)
% Get the degree of the gcd of fx and gx, by sylvester matrix method.
%
%   Inputs
%
%   fw :
%
%   gw :

global PLOT_GRAPHS

% Get Degree of polynomial f(x)
[r,~] = size(fw);
m = r - 1;

% Get degree of polynomial g(x)
[r,~] = size(gw);
n = r - 1;

min_sing_val = zeros(1,min(m,n));
min_distance = zeros(1,min(m,n));

matrix = [];

% Build the Sylvester Matrix
C_f = BuildC1(fw,n,1);
C_g = BuildC1(gw,m,1);

Sk = [C_f C_g];
% Get QR Decomposition of Sk
[Q,R] = qr(Sk);

for k = 1:1:min(m,n)
    
    if k > 1
        
        % update C_f
        C_f = C_f(1:m+n-k+1,1:n-k+1);
        C_g = C_g(1:m+n-k+1,1:m-k+1);
        Sk = [C_f C_g];
        
        % Perform QR Decomposition of Sk, by QR delete.
        % Remove the last column
        [Q,R] = qrdelete(Q,R,m+n+2-((2*k)-2),'col');
        % Remove last column of C_{1}(f)
        [Q,R] = qrdelete(Q,R,n+2-k,'col');
        % Remove last row
        [Q,R] = qrdelete(Q,R,m+n+2-k,'row');
        
        %[Q,R] = qr(Sk)
    end
    
    
    % Get the diagonal values in the R Matrix
    diags = diag(R);
    num_diags = size(diags,1);
    % Save the diagonal values
    matrix = [matrix ;k.*ones(num_diags,1) diags];
    
    
    min_sing_val(k) = min(svd(Sk));
    
    min_distance(k) = getMinDistance(Sk);
    
end


% Check to see if more than one subresultant exsits
if min(m,n) == 1
    % Only one subresultant
    fprintf('Only one subresultant exists')
    switch PLOT_GRAPHS
        case 'y'
            
            figure('name','getDegree - Singular values of Sk')
            hold on
            title('Singular values of S_{1}')
            plot(log10(svd(Sk)))
            hold off
        case 'n'
    end
    
    % Get maximum difference between the singular values
    [val,index] = max(abs(diff(log10(svd(Sk)))));
    THRESHOLD = 1.5;
    if val < THRESHOLD
        % The subresultant is of full rank, in which case t = 0
        t = 0;
        fprintf('The only Subresultant S_{1} appears to be of full rank.');
        fprintf('Polynomials f(x) and g(x) are coprime');
        return
    else % val > threshold
        % The Subresultant is rank deficient, in which case t = 1
        t = 1;
        fprintf('The only Subresultant S_{1} appears to be rank deficient');
        fprintf('Polynomials f(x) and g(x) have a GCD of degree 1');
        return
    end
    
end

switch PLOT_GRAPHS
    case 'y'
        figure('name','QR Scatter data')
        hold on
        scatter(matrix(:,1),log10(matrix(:,2)))
        hold off
        
        % plot the minimum singular values
        figure('name','Plot Min Sing val')
        hold on
        plot(log10(min_sing_val),'-o')
        hold off
        
        figure('name','Plot Min Distance')
        hold on
        plot(log10(min_distance) ,'-s')
        hold off
        
    case 'n'
    otherwise
        error('error: plot_graphs either y or n')
        
end



differences = diff(log10(min_sing_val));
differences(isinf(differences)) = 0;

[~, index] = max(differences);

% if the maximum difference is less than a threshold then all subresultants
% are of full rank or all subresultants are rank deficient.
threshold = 1.9;
criterion = (abs(max(differences)));

if  criterion < threshold
    fprintf('All subresultants are either full rank, or rank deficient\n')
    
    if norm(min_sing_val) < 0
        % if all singular values are close to zero, then rank deficient, degree of
        % gcd is min(m,n)
        t = min(m,n)
    else
        % if all singular values are not close to zero, then full rank, degree
        % of gcd is 0
        
        t = 0;
    end
else
    fprintf('Significant change in min singular values \n')
    t = index;
end

fprintf('The degree of the GCD is %i \n',t)














end