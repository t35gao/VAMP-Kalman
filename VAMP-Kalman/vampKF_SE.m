function [etaX, mseX] = vampKF_SE(svdA, gmaV, gmaV0, svdF, beta, gmaSI, mseSI, svdAF, priorW, K_iter)
svdH = svdA;
svdH.matrix = sqrt(1-beta^2) * svdA.matrix;
svdH.s = sqrt(1-beta^2) * svdA.s;

sAF = svdAF.s;

Ny = size(svdA.matrix, 1);

gmaB = 1 / (beta^2*(1/gmaSI)*(sAF')*sAF/Ny + (1/gmaV));
gmaB0 = 1 / (beta^2*(1/gmaSI)*(sAF')*sAF/Ny + (1/gmaV0));

%% Simulation
% VAMP
[gma2, tau2, mseW] = algo_vamp_SE(svdH, gmaB, gmaB0, priorW, K_iter);

% Kalman
[etaX, mseX] = algo_kalman_diag_SE(svdA, gmaV, gmaV0, svdF, beta, gmaSI, mseSI, gma2, tau2);

end

