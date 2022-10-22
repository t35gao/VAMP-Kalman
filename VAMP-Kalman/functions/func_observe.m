function [gmaV, v, y] = func_observe(A, x, snr)

snr2gmaV = @(sig, snr) 1 / (mean(abs(sig).^2, 'all') * 10^(-snr/10));

Ny = size(A, 1);

y0 = A * x;
gmaV = snr2gmaV(y0, snr);
v = sqrt(1/(2*gmaV)) * (randn(Ny, 1) + 1i*randn(Ny, 1));
y = y0 + v;
end

