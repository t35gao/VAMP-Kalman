clear all;
clear all variables;
close all;
clc;

randnCplx = @(numRow, numCol) randn(numRow, numCol) + 1i*randn(numRow, numCol);
normalizeByPower = @(X, numRow) sqrt(numRow) * X / norm(X, 'fro');
computePower = @(X) mean(abs(X).^2, 'all');
computeVar = @(X) mean(abs(X-mean(X, 'all')).^2, 'all');
computeNRMSE = @(X1, X0) sqrt( mean(abs(X1-X0).^2, 'all') / mean(abs(X0).^2, 'all') );

addpath('functions')
addpath('algorithms')
addpath('bunchmark')

prior_type = 1;

%% Initial Parameters
Nx = 10;
Ny_init = Nx*2;
Ny = 12;

snr_options = 5:5:30;
N_MonteCarlo = 10;
N_innovation = 10;
K_iter = 500;

% innovation factor
beta = 0.5;

% process matrix
F = normalizeByPower(randnCplx(Nx, Nx), Nx);    
svdF = func_svdByRank(F);

% process noise
gmaW = 1;
switch prior_type
    case 1
        muW = 0;
        priorW = Prior_GaussianCplx(muW, gmaW);        
    case 2
        rhoBG = 0.9;
        gmaBG = (1-rhoBG) * gmaW;
        priorW = Prior_BGCplx(rhoBG, gmaBG);
end


% error storage matrices
size_storage = [N_innovation, N_MonteCarlo, length(snr_options)];
eX_KF           = zeros(size_storage);
eX_GSF          = zeros(size_storage);
eX_EnKF         = zeros(size_storage);
eX_CF           = zeros(size_storage);
eX_MCCCF        = zeros(size_storage);
eX_vampKF       = zeros(size_storage);
eX_vampKF_SE    = zeros(size_storage);


%% Monte Carlo Simulations
for n_snr = 1:length(snr_options)
    snr = snr_options(n_snr);
    
    for n_mc = 1:N_MonteCarlo
        %% Initial condition
        A_init = normalizeByPower(randnCplx(Ny_init, Nx), Ny_init);
        svdA_init = func_svdByRank(A_init);

        priorX_init = priorW;
        x_init = priorX_init.generateRand(Nx);
        
        [gmaV, v, y] = func_observe(A_init, x_init, snr);
        
        [xHat_init, etaX_init, ~, ~] = algo_vamp_SVD(... 
            svdA_init, y, gmaV, priorX_init, x_init, K_iter, 1);
        
        x = x_init;
        
        % Update KF
        xHat_KF_old = xHat_init;        etaX_KF_old = etaX_init;
        
        % Update GSF
        xHat_GSF_old = xHat_init;       P_GSF_old = (1/etaX_init)*eye(Nx);
        
        % Update EnKF
        xHat_EnKF_old = xHat_init;      P_EnKF_old = (1/etaX_init)*eye(Nx);
        N_EnKF = 1e4;
        no = sqrt(P_EnKF_old) * repmat(randnCplx(Nx, 1), 1, N_EnKF);
        X_EnKF_old = repmat(xHat_EnKF_old, 1, N_EnKF) + no;
        
        % Update CF
        xHat_CF_old = xHat_init;
        
        % Update MCCKF
        xHat_MCCCF_old = xHat_init;     P_MCCCF_old = (1/etaX_init)*eye(Nx);
        
        % Update VAMP-KF
        xHat_vampKF_old = xHat_init;       etaX_vampKF_old = etaX_init;
               
        %gmaSI_vampKF_SE = gmaSI;
        %mseSI_vampKF_SE = 1/gmaSI_vampKF;        
        
        %% Temporal innovations
        A = normalizeByPower(randnCplx(Ny, Nx), Ny);
        svdA = func_svdByRank(A);

        AF = A*F;
        svdAF = func_svdByRank(AF);                
                
        for n_inov = 1:N_innovation
            w = priorW.generateRand(Nx);
            x = beta*F*x + sqrt(1-beta^2)*w;
            [gmaV, v, y] = func_observe(A, x, snr);
            gmaV0 = gmaV;
            
            Q = (1-beta^2)*(1/gmaW)*eye(Nx);    % process noise
            R = (1/gmaV)*eye(Ny);               % measurement noise
            
            %% Kalman Filter (KF)
            [xHat_KF, etaX_KF] = algo_kalman_diag(svdA, y, gmaV, svdF, beta, xHat_KF_old, etaX_KF_old, 0, gmaW);
            xHat_KF_old = xHat_KF;                  etaX_KF_old = etaX_KF;
            
            %% Adaptive Gaussian Sum Filter (GSF)
            [xHat_GSF, P_GSF] = GSF(xHat_GSF_old, P_GSF_old, y, beta*F, A, Q, R, Nx, Ny);
            xHat_GSF_old = xHat_GSF;                P_GSF_old = P_GSF; 
            
            %% Ensemble Kalman Filter (EnKF)
            [xHat_EnKF, X_EnKF] = EnKF(X_EnKF_old, y, beta*F, A, Q, R, N_EnKF);
            X_EnKF_old = X_EnKF;
            
            %% C Filter (CF)
            [xHat_CF] = CF(xHat_CF_old, y, beta*F, A);
            xHat_CF_old = xHat_CF;
            
            %% MCC-CF
            [xHat_MCCCF, P_MCCCF] = MCCCF(xHat_MCCCF_old, P_MCCCF_old, y, beta*F, A, Q, R, Nx);
            xHat_MCCCF_old = xHat_MCCCF;            P_MCCCF_old = P_MCCCF; 
            
            %% VAMP-Kalman
            [xHat_vampKF, etaX_vampKF] = vampKF(svdA, y, gmaV, svdF, beta, xHat_vampKF_old, etaX_vampKF_old, svdAF, priorW, w, K_iter, 1);            
            xHat_vampKF_old = xHat_vampKF;          etaX_vampKF_old = etaX_vampKF;                                   
                       
            eX_KF(n_inov, n_mc, n_snr) = computeNRMSE(xHat_KF, x);   
            eX_GSF(n_inov, n_mc, n_snr) = computeNRMSE(xHat_GSF, x); 
            eX_EnKF(n_inov, n_mc, n_snr) = computeNRMSE(xHat_EnKF, x);
            eX_CF(n_inov, n_mc, n_snr) = computeNRMSE(xHat_CF, x);
            eX_MCCCF(n_inov, n_mc, n_snr) = computeNRMSE(xHat_MCCCF, x);
            eX_vampKF(n_inov, n_mc, n_snr) = computeNRMSE(xHat_vampKF, x);
            
            % VAMP-Kalman SE
            %{
            [etaX_vampKF_SE, mseX_vampKF_SE] = cntr_gkalman_SE(svdA, gmaV, gmaV0, svdF, beta, gmaSI_vampKF_SE, mseSI_vampKF_SE, svdAF, priorW, K_iter);
            mseSI_vampKF_SE = mseX_vampKF_SE;
            eX_vampKF_SE(n_inov, n_mc, n_snr) = sqrt(mseX_vampKF_SE(end) / computePower(x));
            %}
            
            fprintf(1, '[snr = %02d][MC = %03d][inov = %03d]: ', snr, n_mc, n_inov);
            fprintf(1, 'eX_KF = %f', eX_KF(n_inov, n_mc, n_snr));
            fprintf(1, ', eX_GSF = %f', eX_GSF(n_inov, n_mc, n_snr));
            fprintf(1, ', eX_EnKF = %f', eX_EnKF(n_inov, n_mc, n_snr));
            fprintf(1, ', eX_CF = %f', eX_CF(n_inov, n_mc, n_snr));
            fprintf(1, ', eX_MCCCF = %f', eX_MCCCF(n_inov, n_mc, n_snr));
            fprintf(1, ', eX_vampKF = %f', eX_vampKF(n_inov, n_mc, n_snr));
            %fprintf(1, ', eX_SE = %f', eX_vampKF_SE(n_inov, n_mc, n_snr));
            fprintf(1, '\n');
        end
    end
end


%% Plotting
plot(snr_options, squeeze(mean(eX_KF, [1 2])), 'DisplayName', 'KF');                    hold on
plot(snr_options, squeeze(mean(eX_GSF, [1 2])), 'DisplayName', 'GSF');                 	hold on
plot(snr_options, squeeze(mean(eX_EnKF, [1 2])), 'DisplayName', 'EnKF');             	hold on
plot(snr_options, squeeze(mean(eX_CF, [1 2])), 'DisplayName', 'CF');                    hold on
plot(snr_options, squeeze(mean(eX_MCCCF, [1 2])), 'DisplayName', 'MCCCF');              hold on
plot(snr_options, squeeze(mean(eX_vampKF, [1 2])), 'DisplayName', 'vampKF');            hold on
%plot(snr_options, squeeze(mean(eX_vampKF_SE, [1 2])), 'k', 'DisplayName', 'eX-SE');     hold on
hold off
title('NRMSE vs SNR')
xlabel('NRMSE')
ylabel('SNR')
legend show



