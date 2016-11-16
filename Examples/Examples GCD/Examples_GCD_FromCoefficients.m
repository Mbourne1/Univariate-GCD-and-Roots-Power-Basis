
function [fx,gx,dx,ux,vx] = Examples_GCD_FromCoefficients(ex_num)

addpath('../Examples')
[f,g,d,u,v] = Univariate_GCD_Examples(ex_num);


fx = GetCoefficientsFromSymbolicRoots(f);
gx = GetCoefficientsFromSymbolicRoots(g);
dx = GetCoefficientsFromSymbolicRoots(d);
ux = GetCoefficientsFromSymbolicRoots(u);
vx = GetCoefficientsFromSymbolicRoots(v);
end


