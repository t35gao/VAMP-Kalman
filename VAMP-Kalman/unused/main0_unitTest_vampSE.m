clear all;
clear all variables;
close all;
clc;

randnCplx = @(numRow, numCol) randn(numRow, numCol) + 1i*randn(numRow, numCol);
normalizeByPower = @(X, numRow) sqrt(numRow) * X / norm(X, 'fro');
computePower = @(X) mean(abs(X).^2, 'all');
computeVar = @(X) mean(abs(X-mean(X, 'all')).^2, 'all');
computeRMSE = @(X1, X0) sqrt( mean(abs(X1-X0).^2, 'all'));
computeNRMSE = @(X1, X0) sqrt( mean(abs(X1-X0).^2, 'all') / mean(abs(X0).^2, 'all') );

Nx = 4;
Ny = 6;

snr_options = 5:5:30;
K_iter_VAMP = 500;
K_iter_SE = 100;
N_MonteCarlo = 1;
N_innovation = 10;

gmaX = 1;
rhoBG = 0.5;
gmaBG = (1-rhoBG) * gmaX;
priorX = Prior_BGCplx(rhoBG, gmaBG);

eX = zeros(N_innovation, N_MonteCarlo, length(snr_options));
eX_SE = zeros(N_innovation, N_MonteCarlo, length(snr_options));

for iSNR = 1:length(snr_options)
    snr = snr_options(iSNR);
    
    for n_mc = 1:N_MonteCarlo
        A = normalizeByPower(randnCplx(Ny, Nx), Ny);
        svdA = func_svdByRank(A);
        
        
        n_inov_tame = 1;
        n_inov_wild = 1;
        while n_inov_wild <= N_innovation
            x = priorX.generateRandForced(Nx);
            [gmaV, v, y] = func_observe(A, x, snr);
            %gmaV0 = 1/computeVar(v);
            gmaV0 = gmaV;
            
            % VAMP
            [xHat, etaX, ~, ~, nrmseX, nrmseAx] = algo_vamp_SVD(svdA, y, gmaV, priorX, x, K_iter_VAMP, 0.9);
            
            % VAMP SE
            %[~, ~, mseX_SE] = algo_vamp_SE(svdA, gmaV, gmaV0, priorX, K_iter_SE);
            mseX_SE = 0;
            
            if isnan(sqrt(mseX_SE(end))) || sqrt(mseX_SE(end)) > 1
                N_innovation = N_innovation + 1;
            else
                eX(n_inov_tame, n_mc, iSNR) = computeRMSE(xHat, x);
                eX_SE(n_inov_tame, n_mc, iSNR) = sqrt(mseX_SE(end));
                n_inov_tame = n_inov_tame + 1;
            end
            
            fprintf(1, '[snr = %02d][MC = %03d][inov = %03d]: ', snr, n_mc, n_inov_wild);
            fprintf(1, 'gmaV/V0 = %f, ', gmaV/gmaV0);
            fprintf(1, 'nrmseX = %f, \n', computeRMSE(xHat, x));
            %fprintf(1, 'nrmseX_SE = %f\n', sqrt(mseX_SE(end)));
            
            n_inov_wild = n_inov_wild + 1;
        end
        N_innovation = 10;
    end
end

plot(snr_options, squeeze(mean(eX, [1 2])), 'b', 'DisplayName', 'VAMP')
hold on
plot(snr_options, squeeze(mean(eX_SE, [1 2])), 'k', 'DisplayName', 'VAMP SE')
hold off
legend show