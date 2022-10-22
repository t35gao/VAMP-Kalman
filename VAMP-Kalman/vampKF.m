function [xHat, etaX] = vampKF(svdA, y, gmaV, svdF, beta, rSI, gmaSI, svdAF, priorW, w, K_iter, ratioDamp)
svdH = svdA;
svdH.matrix = sqrt(1-beta^2) * svdA.matrix;
svdH.s = sqrt(1-beta^2) * svdA.s;

AF = svdAF.matrix;      sAF = svdAF.s;

Ny = length(y);

z = y - beta * AF * rSI;
gmaB = 1 / (beta^2*(1/gmaSI)*(sAF')*sAF/Ny + (1/gmaV));

%% Simulation
% VAMP
[wHat, etaW, r2, gma2] = algo_vamp_SVD(svdH, z, gmaB, priorW, w, K_iter, ratioDamp);

% Kalman
[xHat, etaX] = algo_kalman_diag(svdA, y, gmaV, svdF, beta, rSI, gmaSI, r2, gma2);

end

