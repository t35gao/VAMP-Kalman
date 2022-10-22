function [etaEst_new, mseX_new] = algo_kalman_diag_SE(svdA, gmaV, gmaV0, svdF, beta, etaEst_old, mseX_old, gmaW, tauW)
A = svdA.matrix;
Nx = size(A, 2);

%% Prediction
sF = svdF.s;

gmaPrd = 1 / (beta^2*(1/etaEst_old)*(sF')*sF/Nx + (1/gmaW));

tauPrd = beta^2*mseX_old*(sF')*sF/Nx + (1-beta^2)*tauW;

%% Update
sA = svdA.s;    sA(end+1:Nx) = 0;

etaEst_new = 1 / mean(1./(gmaV*abs(sA).^2 + gmaPrd), 'all');

mseX_new = mean((gmaV^2*abs(sA).^2/gmaV0 + gmaPrd^2*tauPrd)./ (gmaV*abs(sA).^2 + gmaPrd).^2, 'all');

end