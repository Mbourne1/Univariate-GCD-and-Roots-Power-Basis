function t = GetDegree(fw,gw)
% Get the degree of the gcd of fx and gx, by sylvester matrix method.
%
%   Inputs
%
%   fw :
%
%   gw :

global PLOT_GRAPHS

% Set the threshold for determining whether all subresultants are rank
% defficient or all subresultants are full rank, OR some subresultants are
% rank deficient, and others are of full rank, meaning degree of GCD can be
% identified.
THRESHOLD = 1;


% Get Degree of polynomial f(x)
[r,~] = size(fw);
m = r - 1;

% Get degree of polynomial g(x)
[r,~] = size(gw);
n = r - 1;

% Initialise a vector to store the minimum singular values for each
% S_{k}.
vMin_sing_val = zeros(1,min(m,n));

% Initialise a vector to store the minimum distnaces for each S_{k}.
vMin_distance = zeros(1,min(m,n));

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
    
    % Add to the vector of minimum Singular values.
    vMin_sing_val(k) = min(svd(Sk));
    
    % Add to the vector of minimum distances
    vMin_distance(k) = GetMinDistance(Sk);
    
end


% Check to see if more than one subresultant exsits
if min(m,n) == 1
    % Only one subresultant
    fprintf('Only one subresultant exists. \n')
    switch PLOT_GRAPHS
        case 'y'
            
            figure('name','getDegree - Singular values of Sk')
            hold on
            title('Singular values of S_{1}')
            plot(log10(svd(Sk)))
            hold off
        case 'n'
    end
    
    % Get the differences between minimum singular value S_{k} and S_{k+1}
    differences = diff(log10(vMin_sing_val));
    
    
    % Get the index of the largest change in minimum singular values
    [delta_min_sing_val, index] = max(differences);
    
    if delta_min_sing_val < THRESHOLD
        % The subresultant is of full rank, in which case t = 0
        t = 0;
        fprintf('The only Subresultant S_{1} appears to be of full rank. \n');
        fprintf('Polynomials f(x) and g(x) are coprime \n');
        return
    else % val > threshold
        % The Subresultant is rank deficient, in which case t = 1
        t = 1;
        fprintf('The only Subresultant S_{1} appears to be rank deficient \n');
        fprintf('Polynomials f(x) and g(x) have a GCD of degree 1 \n');
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
        plot(log10(vMin_sing_val),'-o')
        hold off
        
        figure('name','Plot Min Distance')
        hold on
        plot(log10(vMin_distance) ,'-s')
        hold off
        
    case 'n'
    otherwise
        error('error: plot_graphs either y or n')
        
end


% Get the differences between minimum singular value S_{k} and S_{k+1}
differences = abs(diff(log10(vMin_sing_val)));
differences(isnan(differences)) = 0 ;
differences(differences==inf)=0;
% Get the index of the largest change in minimum singular values
[delta_min_sing_val, index] = max(differences);



% if the maximum difference is less than a threshold then all subresultants
% are of full rank or all subresultants are rank deficient.


criterion = abs(delta_min_sing_val);

% Check to see if the largest change in minimum singular value from S_{k}
% to S_{k+1} is significant. If not significant, must assume that all
% subresultants are either ALL RANK DEFICIENT or ALL FULL RANK

if  criterion < THRESHOLD
    fprintf('All subresultants are either full rank, or rank deficient\n')
    
    avgMinSingularValue = log10(norm(vMin_sing_val))
    
    if  avgMinSingularValue < -12
        % if all singular values are close to zero, then rank deficient, degree of
        % gcd is min(m,n)
        t = min(m,n);
    else
        % if all singular values are not close to zero, then full rank, degree
        % of gcd is 0
        
        t = 0;
    end
else
    %fprintf('The minimal singular values of S_{k} are : \n')
    %display(vMin_sing_val)
    fprintf('Significant change in min singular values \n')
    t = index;
end

fprintf('The degree of the GCD is %i \n',t)














end