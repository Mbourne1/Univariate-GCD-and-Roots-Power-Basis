function [res,index] = GetMinDistance(Sk)
% Given a sylvester matrix get the minimum distance for the removal of a
% column

[~,c] = size(Sk);

residual_vec = zeros(c,1);

for i = 1:1:c
    
    % Get Ak, Sk with ck removed
    Ak = Sk;
    Ak(:,i) = [];
    
    % Get the column ck
    ck = Sk(:,i);
    
    inversion_method = 'QR';
    %inversion_method = 'SVD'
    switch inversion_method
        case 'QR'
            [~,n2] = size(Ak);
            [Q2,R] = qr(Ak);
            R1 = R(1:n2,:);
            cd = Q2'*ck;
            c1 = cd(1:n2,:);
            x_ls = R1\c1;
        case 'SVD'
            x_ls = pinv(Ak) * ck;
    end
    
    % Get the residual
    residual_vec(i) = norm(ck - (Ak*x_ls));
    
end

% Get the minimal residual and the index of the corresponding column ck
[res,index] = min(residual_vec);

end