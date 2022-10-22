function [xEst_new, etaEst_new] = algo_kalman_diag(svdA, y, gmaV, svdF, beta, xEst_old, etaEst_old, rW, gmaW)
Nx = length(xEst_old);

%% Prediction
F = svdF.matrix;    sF = svdF.s;

rPrd = beta * F * xEst_old + sqrt(1-beta^2) * rW;
gmaPrd = 1 / (beta^2*(1/etaEst_old)*(sF')*sF/Nx + (1/gmaW));

%% Update
UA = svdA.U;        VA = svdA.V;        sA = svdA.s;
rankA = length(sA);
yTld = diag(1./sA)*UA'*y;

d = gmaV*abs(sA).^2 ./ (gmaV*abs(sA).^2 + gmaPrd);
xEst_new = rPrd + VA*(d.* (yTld - VA'*rPrd));
etaEst_new = gmaPrd / (1 - (rankA/Nx)*mean(d));
end