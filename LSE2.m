function y = LSE2(E,f,C,g)
% This function uses the QR decomposition to solve the LSE problem
% minimise ||Ey-f||  subject to  Cy=g
% The output is the vector y.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % a.) Compute the QR Decomposition of CT
    [m1,~]=size(C);
    [Q,R] = qr(C');
    
    % b.) Set w_{1} = R_{1}^{-T} * q
    w1 = R1'\g;
    
    % Partition EQ as EQ = [E1 E2]
    
    
    %Q1 = Q(:,1:m1);
    Q2 = Q(:,m1+1:end);

    w2 = Q2'*f;

    y = Q*[w1;w2];



end
