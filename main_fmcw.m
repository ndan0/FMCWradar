clear; clc; close all;

%% =========================================================
%  FMCW Radar Simulation up to Fused 2D CFAR Detection
%
%  "An Improved Scheme for High-Resolution Point Cloud Map
%   Generation Based on FMCW Radar"
%
%  stages:
%   1) IF signal simulation per receiver channel
%   2) Range FFT
%   3) Doppler FFT
%   4) Fused Range Doppler Map RDM Z(kv,kr)
%   5) Improved 2D CFAR detection

%
%  No clustering
%  No AOA
%% =========================================================

%% =========================
% 1) Radar Parameters
% ==========================
c  = 3e8;                  % speed of light (m/s)
fc = 77e9;                 % carrier frequency (Hz)
lambda = c / fc;           % wavelength (m)

B  = 2.5e9;                % bandwidth (Hz)
S  = 50.018e12;            % chirp slope (Hz/s)
Tchirp = 130e-6;           % chirp repetition period (s)

N   = 256;                 % fast-time samples per chirp
M   = 128;                 % number of chirps
Nrx = 8;                   % simulated receiver channels

Ts = Tchirp / N;           % fast-time sampling interval
fs = 1 / Ts;               % fast-time sampling rate
d  = lambda / 2;           % antenna spacing

%% =========================
% 2) Ground-Truth Targets
% ==========================
R_targets         = [2.0, 2.8, 1.4];     % meters
v_targets         = [1.0, -0.8, 0.0];    % m/s
theta_targets_deg = [20, -15, 35];       % degrees
A_targets         = [1.0, 0.85, 0.65];   % amplitudes

SNR_dB = 10;
numTargets = numel(R_targets);

fprintf('\n=== Ground Truth Targets ===\n');
for k = 1:numTargets
    fprintf('Target %d: R = %.3f m, v = %.3f m/s, theta = %.3f deg, A = %.2f\n', ...
        k, R_targets(k), v_targets(k), theta_targets_deg(k), A_targets(k));
end
fprintf('\n');

%% =========================
% 3) Derived Metrics
% ==========================
range_resolution    = c / (2 * B);
Rmax                = c * fs / (4 * S);
velocity_resolution = lambda / (2 * M * Tchirp);
vmax                = lambda / (4 * Tchirp);

fprintf('=== Radar Metrics ===\n');
fprintf('Range resolution          = %.4f m\n', range_resolution);
fprintf('Max unambiguous range     = %.4f m\n', Rmax);
fprintf('Velocity resolution       = %.4f m/s\n', velocity_resolution);
fprintf('Max unambiguous velocity  = %.4f m/s\n\n', vmax);

%% =========================
% 4) Indices
% ==========================
n = (0:N-1).';             % fast-time index
m = 0:M-1;                 % slow-time index

%% =========================
% 5) Simulate IF Signal Cube xIF(n,m,rx)
% ==========================
x = zeros(N, M, Nrx);

for k = 1:numTargets
    Rk = R_targets(k);
    vk = v_targets(k);
    thetak = deg2rad(theta_targets_deg(k));
    Ak = A_targets(k);

    fr_k     = 2 * S * Rk / c;               % range frequency
    fd_k     = 2 * fc * vk / c;              % Doppler frequency
    fTheta_k = d * sin(thetak) / lambda;     % spatial frequency

    for rx = 1:Nrx
        phase_space = 2*pi*fTheta_k*(rx-1);

        sig_nm = Ak * exp(1j * ( ...
            2*pi*fr_k*n*Ts + ...
            2*pi*fd_k*m*Tchirp + ...
            phase_space ));

        x(:,:,rx) = x(:,:,rx) + sig_nm;
    end
end

%% =========================
% 6) Add Complex AWGN
% ==========================
signalPower = mean(abs(x(:)).^2);
noisePower  = signalPower / (10^(SNR_dB/10));
noise = sqrt(noisePower/2) * (randn(size(x)) + 1j*randn(size(x)));
x_noisy = x + noise;

%% =========================
% 7) Range FFT
% ==========================
wr = hann(N);
rangeFFT = zeros(N, M, Nrx);

for rx = 1:Nrx
    xw = x_noisy(:,:,rx) .* wr;
    rangeFFT(:,:,rx) = fft(xw, [], 1);
end

%% =========================
% 8) Doppler FFT
% ==========================
wd = hann(M).';
RD = zeros(N, M, Nrx);

for rx = 1:Nrx
    temp = rangeFFT(:,:,rx) .* wd;
    RD(:,:,rx) = fftshift(fft(temp, [], 2), 2);
end

%% =========================
% 9) Axes
% ==========================
f_range_axis = (0:N-1) * (fs / N);
R_axis = c * f_range_axis / (2 * S);

fd_axis = (-M/2:M/2-1) * (1 / (M * Tchirp));
v_axis  = (c * fd_axis) / (2 * fc);

Nr = floor(N/2);  % positive range bins only

%% =========================
% 10) Single-Channel RDM (for reference only)
% ==========================
RD_single = abs(RD(1:Nr,:,1)).^2;

figure;
imagesc(v_axis, R_axis(1:Nr), 10*log10(RD_single / max(RD_single(:)) + eps));
xlabel('Velocity (m/s)');
ylabel('Range (m)');
title('Single-Channel Range-Doppler Map (Rx 1)');
colorbar;
axis xy;
ylim([0, Rmax]);

%% =========================
% 11) Fused RDM Z(kv,kr)
% ==========================
Z = mean(abs(RD(1:Nr,:,:)).^2, 3);

figure;
imagesc(v_axis, R_axis(1:Nr), 10*log10(Z / max(Z(:)) + eps));
xlabel('Velocity (m/s)');
ylabel('Range (m)');
title('Fused Range-Doppler Map Z(k_v,k_r)');
colorbar;
axis xy;
ylim([0, Rmax]);

%% =========================
% 12) Improved 2D CFAR
% ==========================
window = [11, 11];   % odd-sized rectangular window [rows, cols]
Tscale = 3.3;        % threshold scale factor

[peakBins, detMask, Umap, Smap] = cfar_detection_paper(Z, Tscale, window);

fprintf('Raw CFAR detections = %d\n', size(peakBins,1));

%% =========================
% 13) Plot CFAR detections on fused RDM
% ==========================
figure;
imagesc(v_axis, R_axis(1:Nr), 10*log10(Z / max(Z(:)) + eps));
xlabel('Velocity (m/s)');
ylabel('Range (m)');
title('Fused RDM with Improved 2D CFAR Detections');
colorbar;
axis xy;
ylim([0, Rmax]);
hold on;
plot(v_axis(peakBins(:,2)), R_axis(peakBins(:,1)), 'ro', ...
    'MarkerSize', 6, 'LineWidth', 1.2);
hold off;

%% =========================
% 14) Display detections
% ==========================
fprintf('\n=== Paper-style CFAR detections ===\n');
for k = 1:size(peakBins,1)
    rbin = peakBins(k,1);
    dbin = peakBins(k,2);
    fprintf('%2d: R = %.3f m, v = %.3f m/s, power = %.6f\n', ...
        k, R_axis(rbin), v_axis(dbin), Z(rbin, dbin));
end
fprintf('\n');
