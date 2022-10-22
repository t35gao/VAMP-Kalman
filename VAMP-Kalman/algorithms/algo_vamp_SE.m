function [gma2, tau2, mse] = algo_vamp_SE(svdA, gmaV, gmaV0, priorX, K_iter)
A = svdA.matrix;
Nx = size(A, 2);
s = svdA.s;
s(end+1:Nx) = 0;

N_sample = 5e3;
x0 = priorX.generateRandForced(N_sample);

tau1 = 1;
gma1 = 1;

mse = zeros(K_iter, 1);

for k = 1:K_iter
    %% MMSE block
    r1 = x0 + sqrt(tau1/2)*(randn(N_sample, 1) + 1i*randn(N_sample, 1));
    [g1, g1Prime] = priorX.gfuncAndPrime(r1, gma1);
    
    alp1 = mean(g1Prime, 'all');
    eta1 = gma1/alp1;           
    
    gma2 = eta1-gma1;
    
    err1 = mean(abs(g1-x0).^2, 'all');
    
    tau2 = (err1 + alp1^2*tau1 - alp1^2*(1/gma1)) / (1-alp1)^2;
    
    %% Compute MSE
    mse(k) = err1;
    
    %% LMMSE block
    alp2 = mean(gma2./(gmaV*abs(s).^2 + gma2), 'all');
    eta2 = gma2/alp2;           
    
    gma1 = eta2-gma2;
    
    err2 = mean((gmaV^2*abs(s).^2/gmaV0 + gma2^2*tau2)./ (gmaV*abs(s).^2 + gma2).^2, 'all');
    
    tau1 = (err2 - alp2^2*tau2) / (1-alp2)^2;
end
end

