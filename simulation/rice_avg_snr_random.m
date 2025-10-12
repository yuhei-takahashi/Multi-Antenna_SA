function snr_samples = rice_avg_snr_random(K, avg_snr_dB, num_samples)
% K: Rician K-factor
% avg_snr_dB: Average received SNR in dB
% num_samples: Number of random samples to generate

% Convert average SNR from dB to linear scale
avg_snr_linear = 10^(avg_snr_dB / 10);

% Compute sigma^2 based on K-factor
sigma = sqrt(avg_snr_linear / (2 * (K + 1)));

% Compute nu (the LOS component)
nu = sqrt(2 * K * sigma^2);

% Generate random samples following Rician distribution
fading_samples = rice_rnd(nu, sigma, num_samples);

% Compute SNR as the squared magnitude of the fading samples
snr_samples = (abs(fading_samples).^2);
end

function r = rice_rnd(nu, sigma, num_samples)
% nu: Noncentrality parameter (mean of LOS component)
% sigma: Standard deviation of scattered component
% num_samples: Number of random samples to generate

% Generate in-phase (I) and quadrature (Q) components from Gaussian distributions
X1 = nu + sigma * randn(1, num_samples);  % In-phase component
X2 = sigma * randn(1, num_samples);       % Quadrature component

% Compute amplitude (magnitude of complex value)
r = sqrt(X1.^2 + X2.^2);
end
