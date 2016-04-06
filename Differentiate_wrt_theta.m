function f_wrt_theta = Differentiate_wrt_theta(fw,theta)
% Given the polynomial f(\omega,\theta) differentiate with respect to theta

m = GetDegree(fw);

mat = diag(0:1:m) ./ theta;

f_wrt_theta = mat * fw;


end