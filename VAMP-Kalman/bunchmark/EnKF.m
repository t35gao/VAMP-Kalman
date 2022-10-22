function [xEst, X_enkf] = EnKF(X_enkf, z, A, B, Q, R, N_Enkf)
    X_enkf = A * X_enkf;
    no     = sqrt(Q) * randn(size(X_enkf));
    X_enkf = X_enkf + no;
    xEst   = mean(X_enkf,2);
    PPred  = (1/(N_Enkf-1)) * ((X_enkf - repmat(xEst,1,N_Enkf))*(X_enkf - repmat(xEst,1,N_Enkf)).');
    K      = PPred * B' / (B * PPred * B' + R);
    X_enkf = X_enkf + K * (repmat(z,1,N_Enkf)- B * X_enkf);
    xEst   = mean(X_enkf,2);
    %PEst   = (1/(N_Enkf-1)) * ((X_enkf - repmat(xEst,1,N_Enkf)) * (X_enkf - repmat(xEst,1,N_Enkf)).');
end