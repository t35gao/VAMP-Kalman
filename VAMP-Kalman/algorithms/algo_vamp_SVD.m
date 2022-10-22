function [xHat, eta, rWvy, gmaWvy] = algo_vamp_SVD(svdA, y, gmaV, priorX, x, K_iter, ratioDamp)
A = svdA.matrix;
U = svdA.U;     V = svdA.V;     s = svdA.s;

Nx = size(V, 1);
R = length(s);
yTld = diag(1./ s) * U' * y;
sSqr = abs(s).^2;

r = randn(Nx, 1) + 1i*randn(Nx, 1);
gma = 1;
xHat = 0;

gmaMin = 1e-11;
gmaMax = 1e11;

nrmseX = zeros(K_iter, 1);
nrmseAx = zeros(K_iter, 1);

% Function definition
computeNRMSE = @(X1, X0) sqrt( mean(abs(X1-X0).^2, 'all') / mean(abs(X0).^2, 'all') );

for k = 1:K_iter
    %% Estimation
    [g1, g1Prime] = priorX.gfuncAndPrime(r, gma);
    
    xHat = ratioDamp * g1 + (1-ratioDamp) * xHat;
    
    alp = mean(g1Prime, 'all');
    eta = gma/alp;
    
    rWvy = (xHat - alp*r) / (1 - alp);
    
    gmaWvy = gma * (1 - alp) / alp;
    gmaWvy = max(min(gmaWvy, gmaMax), gmaMin);
    
    d = (gmaV*sSqr)./ (gmaV*sSqr + gmaWvy);
    
    gma = ratioDamp * (gmaWvy * mean(d) / (Nx/R - mean(d))) + (1-ratioDamp) * gma;
    gma = max(min(gma, gmaMax), gmaMin);
    
    r = rWvy + (Nx/R) * V * ((d/mean(d)).* (yTld - V'*rWvy));
    
    %% Error
    % nrmseX(k) = computeNRMSE(xHat, x);
    % nrmseAx(k) = computeNRMSE(A*xHat, A*x);
end
end

