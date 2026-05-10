clear; clc; close all;

%% =========================================================
%  FMCW Radar Simulation with Paper-Style Improved 2D CFAR
%  and Switchable TDM-MIMO Virtual Array
%
%  Supported configurations:
%   1) "9x16" : 9 Tx x 16 Rx, tx_step = 9, 88 virtual elements
%   2) "4x8"  : 4 Tx x 8 Rx,  tx_step = 4, 20 virtual elements
%
%  Flow:
%   1) Define targets and radar system
%   2) Simulate TDM-MIMO FMCW radar returns
%   3) Build overlapping virtual array
%   4) Generate range-Doppler response
%   5) Fuse virtual channels into Z(kv,kr)
%   6) Apply improved 2D CFAR
%   7) Apply post-CFAR grouping / NMS
%   8) Estimate AOA using DFT on virtual array
%   9) Plot PCM: Truth vs Detected
%% =========================================================

%% User parameters

target_dist  = [41   43   31   60];                    % meters
target_speed = [100.0 100.0 109.84 92.0190]*1000/3600;  % m/s
target_az    = [-25.00 -23.70 13.203 16.640];          % degrees
target_rcs   = [100 100 20 20];                        % m^2

radar_speed = 100*1000/3600;                           % m/s

debug_plot = 0;

% Choose antenna configuration:
%   "9x16" -> paper-matching 88 virtual elements
%   "4x8"  -> reduced 20 virtual elements
array_config = "9x16";
% array_config = "4x8";

%% Radar waveform parameters

fc     = 77e9;
c      = physconst('LightSpeed');
lambda = c / fc;

range_max = 200;
T         = 5 * range2time(range_max, c);
range_res = 1;
bw        = rangeres2bw(range_res, c);
S         = bw / T;

fr_max = range2beat(range_max, S, c);
v_max  = 230 * 1000/3600;
fd_max = speed2dop(2*v_max, lambda);
fb_max = fr_max + fd_max;
fs     = max(2*fb_max, bw);

%% Define FMCW pulse

waveform = phased.FMCWWaveform( ...
    'SweepTime', T, ...
    'SweepBandwidth', bw, ...
    'SampleRate', fs);

sig = waveform();

if debug_plot
    figure;
    subplot(211);
    plot(0:1/fs:T-1/fs, real(sig));
    xlabel('Time (s)');
    ylabel('Amplitude (V)');
    title('FMCW Signal');
    axis tight;

    subplot(212);
    spectrogram(sig, 32, 16, 32, fs, 'yaxis');
    title('FMCW Signal Spectrogram');
end

%% Define targets

target_pos = [target_dist.*cosd(target_az);
              target_dist.*sind(target_az);
              0.5*ones(1, length(target_dist))];

target = phased.RadarTarget( ...
    'MeanRCS', target_rcs, ...
    'PropagationSpeed', c, ...
    'OperatingFrequency', fc);

target_motion = phased.Platform( ...
    'InitialPosition', target_pos, ...
    'Velocity', [target_speed;
                 zeros(1,length(target_dist));
                 zeros(1,length(target_dist))]);

channel = phased.FreeSpace( ...
    'PropagationSpeed', c, ...
    'OperatingFrequency', fc, ...
    'SampleRate', fs, ...
    'TwoWayPropagation', true);

for i = 1:length(target_dist)
    fprintf(1, "Target %d at range %.1f m, abs vel %.2f km/h, az %.3f deg, rcs %.0f\n", ...
        i, target_dist(i), target_speed(i)*3.6, target_az(i), target_rcs(i));
end

%% Setup radar system

ant_aperture = 6.06e-4;
ant_gain     = aperture2gain(ant_aperture, lambda);

tx_ppower = db2pow(25) * 1e-3;   % 316 mW
tx_gain   = 9  + ant_gain;
rx_gain   = 15 + ant_gain;
rx_nf     = 4.5;

%% =========================================================
%  Select TDM-MIMO antenna configuration
%% =========================================================

switch array_config

    case "9x16"
        % Paper-matching configuration:
        % 9 Tx x 16 Rx
        % Tx positions: 0,9,18,27,36,45,54,63,72
        % Rx positions: 0,1,2,...,15
        % Nvirt = (9 - 1)*9 + 16 = 88
        Nt = 9;
        Nr = 16;
        tx_step = 9;

    case "4x8"
        % Reduced overlapping configuration:
        % 4 Tx x 8 Rx
        % Tx positions: 0,4,8,12
        % Rx positions: 0,1,2,...,7
        % Nvirt = (4 - 1)*4 + 8 = 20
        Nt = 4;
        Nr = 8;
        tx_step = 4;

    otherwise
        error('Unknown array_config. Use "9x16" or "4x8".');
end

dr = lambda / 2;          % Rx spacing
dt = tx_step * dr;        % Tx spacing
Nrx_virt = (Nt - 1)*tx_step + Nr;

txarray = phased.ULA(Nt, dt);
rxarray = phased.ULA(Nr, dr);

transmitter = phased.Transmitter( ...
    'PeakPower', tx_ppower, ...
    'Gain', tx_gain);

receiver = phased.ReceiverPreamp( ...
    'Gain', rx_gain, ...
    'NoiseFigure', rx_nf, ...
    'SampleRate', fs);

txradiator = phased.Radiator( ...
    'Sensor', txarray, ...
    'OperatingFrequency', fc, ...
    'PropagationSpeed', c, ...
    'WeightsInputPort', true);

rxcollector = phased.Collector( ...
    'Sensor', rxarray, ...
    'OperatingFrequency', fc, ...
    'PropagationSpeed', c);

radar_speed_kmh = radar_speed * 3.6;

radarmotion = phased.Platform( ...
    'InitialPosition', [0;0;0.5], ...
    'Velocity', [radar_speed;0;0]);

fprintf(1, "\nSelected array configuration: %s\n", array_config);
fprintf(1, "Radar: %d Tx, %d Rx, %d Virtual (unique), speed %.1f km/h\n", ...
        Nt, Nr, Nrx_virt, radar_speed_kmh);
fprintf(1, "Tx spacing: %d Rx spacings\n", tx_step);
fprintf(1, "Angle resolution: %.4f deg\n\n", ...
        rad2deg(lambda / (Nrx_virt * dr)));

truth_targets = [target_dist.' , (radar_speed - target_speed).' , target_az.'];

%% Simulation loop

if debug_plot
    specanalyzer = spectrumAnalyzer( ...
        'SampleRate', fs, ...
        'Method', 'welch', ...
        'AveragingMethod', 'running', ...
        'PlotAsTwoSidedSpectrum', true, ...
        'FrequencyResolutionMethod', 'rbw', ...
        'Title', 'Spectrum for Received and Dechirped Signal', ...
        'ShowLegend', true);
end

rng(11);

Dn            = 2;
fs            = fs / Dn;
Nsweep_per_Tx = 32;
Nsweep        = Nt * Nsweep_per_Tx;

xr = complex(zeros(floor(fs*waveform.SweepTime), Nr, Nsweep));

w0     = zeros(Nt, 1);
w0(Nt) = 1;

for m = 1:Nsweep

    [radar_pos, radar_vel]   = radarmotion(waveform.SweepTime);
    [target_pos, target_vel] = target_motion(waveform.SweepTime);

    [~, target_ang] = rangeangle(target_pos, radar_pos);

    sig   = waveform();
    txsig = transmitter(sig);

    % TDM-MIMO Tx switching
    w0    = [w0(Nt); w0(1:Nt-1)];
    txsig = txradiator(txsig, target_ang, w0);

    txsig = channel(txsig, radar_pos, target_pos, radar_vel, target_vel);
    txsig = target(txsig);

    rxsig      = rxcollector(txsig, target_ang);
    rxsig      = receiver(rxsig);
    dechirpsig = dechirp(rxsig, sig);

    for n = size(xr,2):-1:1
        xr(:,n,m) = decimate(dechirpsig(:,n), Dn, 'FIR');
    end

    if debug_plot
        specanalyzer([txsig dechirpsig]);
    end
end

%% Virtual array assembly using overlapping TDM-MIMO demux

xrv_sum   = complex(zeros(size(xr,1), Nrx_virt, Nsweep_per_Tx));
xrv_count = zeros(1, Nrx_virt);

for tx = 1:Nt

    tx_data = xr(:, :, tx:Nt:end);   % [fast_time, Nr, Nsweep_per_Tx]

    for rx = 1:Nr

        virt_idx = (tx-1)*tx_step + rx;

        xrv_sum(:, virt_idx, :) = xrv_sum(:, virt_idx, :) + tx_data(:, rx, :);
        xrv_count(virt_idx) = xrv_count(virt_idx) + 1;

    end
end

xrv = xrv_sum;

for v = 1:Nrx_virt
    if xrv_count(v) > 0
        xrv(:, v, :) = xrv(:, v, :) ./ xrv_count(v);
    end
end

fprintf(1, "Virtual array positions used: 1 to %d\n", Nrx_virt);
fprintf(1, "Virtual overlap counts per position:\n");
disp(xrv_count);

%% Range-Doppler response

nfft_r = 2^nextpow2(size(xrv,1));
nfft_d = 2^nextpow2(size(xrv,3));

rngdopresp = phased.RangeDopplerResponse( ...
    'PropagationSpeed',       c, ...
    'DopplerOutput',          'Speed', ...
    'OperatingFrequency',     fc, ...
    'SampleRate',             fs, ...
    'RangeMethod',            'FFT', ...
    'PRFSource',              'Property', ...
    'RangeWindow',            'Hann', ...
    'PRF',                    1/(Nt*waveform.SweepTime), ...
    'SweepSlope',             S, ...
    'RangeFFTLengthSource',   'Property', ...
    'RangeFFTLength',         nfft_r, ...
    'DopplerFFTLengthSource', 'Property', ...
    'DopplerFFTLength',       nfft_d, ...
    'DopplerWindow',          'Hann');

[resp, r, sp] = rngdopresp(xrv);

if debug_plot
    figure;
    plotResponse(rngdopresp, squeeze(xrv(:,1,:)));
end

%% Fused RDM Z(kv,kr)

mag_resp = abs(resp(nfft_r/2+1:end, :, :)).^2;

Z = reshape(sum(mag_resp, 2), [nfft_r/2, nfft_d]);

r = r(nfft_r/2+1:end);

%% Improved 2D CFAR on fused RDM Z

window = [7, 5];      % [range bins, Doppler bins], must be odd
Tscale = 6.0;         % threshold scaling factor

[detects_raw, detMask, Umap, Smap] = cfar_detection_paper(Z, Tscale, window);

fprintf(1, "Raw CFAR detections: %d\n", size(detects_raw, 1));

if debug_plot
    figure;
    imagesc(sp, r, 10*log10(Z ./ max(Z(:)) + eps));
    axis xy;
    colorbar;
    xlabel('Relative Velocity (m/s)');
    ylabel('Range (m)');
    title('Fused RDM with Raw CFAR Detections');
    hold on;
    if ~isempty(detects_raw)
        plot(sp(detects_raw(:,2)), r(detects_raw(:,1)), ...
             'ro', 'MarkerSize', 7, 'LineWidth', 1.2);
    end
    hold off;
end

%% Post-CFAR peak grouping / Non-Maximum Suppression

if isempty(detects_raw)

    detects = zeros(0,2);
    num_dets = 0;

    fprintf(1, "After peak grouping: 0 detections\n\n");

else

    range_bin_size   = abs(mean(diff(r)));
    doppler_bin_size = abs(mean(diff(sp)));

    group_r_m   = 1.5 * range_bin_size;
    group_d_mps = 1.5 * doppler_bin_size;

    fprintf(1, "Range bin: %.3f m  Doppler bin: %.3f m/s\n", ...
            range_bin_size, doppler_bin_size);
    fprintf(1, "Grouping thresholds: range < %.2f m AND Doppler < %.2f m/s\n", ...
            group_r_m, group_d_mps);

    powers = Z(sub2ind(size(Z), detects_raw(:,1), detects_raw(:,2)));
    [~, sort_idx] = sort(powers, 'descend');

    detects_sorted = detects_raw(sort_idx, :);

    kept = true(size(detects_sorted, 1), 1);

    for i = 1:size(detects_sorted, 1)

        if ~kept(i)
            continue;
        end

        for j = i+1:size(detects_sorted, 1)

            if ~kept(j)
                continue;
            end

            dr_m   = abs(r(detects_sorted(i,1))  - r(detects_sorted(j,1)));
            dd_mps = abs(sp(detects_sorted(i,2)) - sp(detects_sorted(j,2)));

            if dr_m < group_r_m && dd_mps < group_d_mps
                kept(j) = false;
            end
        end
    end

    detects = detects_sorted(kept, :);
    num_dets = size(detects, 1);

    fprintf(1, "After peak grouping: %d detections\n\n", num_dets);
end

if debug_plot
    figure;
    imagesc(sp, r, 10*log10(Z ./ max(Z(:)) + eps));
    axis xy;
    colorbar;
    xlabel('Relative Velocity (m/s)');
    ylabel('Range (m)');
    title('Fused RDM with Final Grouped Detections');
    hold on;
    if ~isempty(detects)
        plot(sp(detects(:,2)), r(detects(:,1)), ...
             'rd', 'MarkerSize', 9, 'LineWidth', 1.5);
    end
    hold off;
end

%% AOA estimation — DFT on virtual array

Nrx = Nrx_virt;

detected_targets = zeros(num_dets, 3);   % [range(m), rel_vel(m/s), angle(deg)]

fprintf(1, "%-8s %-10s %-12s %-14s %-14s\n", ...
        "Detect", "Range(m)", "RelVel(m/s)", "DFT_az(deg)", "Truth_az(deg)");

for i = 1:num_dets

    ri = detects(i,1);
    di = detects(i,2);

    Y = resp(nfft_r/2 + ri, :, di);   % 1 x Nrx_virt

    Pi = fft(Y, Nrx);

    [~, max_idx] = max(abs(Pi));

    if max_idx > Nrx/2
        k_norm = (max_idx - 1 - Nrx) / Nrx;
    else
        k_norm = (max_idx - 1) / Nrx;
    end

    k_norm = -k_norm;

    val = max(min(k_norm / (dr/lambda), 1), -1);
    angle_deg = rad2deg(asin(val));

    v_rel = sp(di);

    [~, ti] = min(abs(target_dist - r(ri)));
    truth_az = target_az(ti);

    fprintf(1, "%-8d %-10.2f %-12.2f %-14.2f %-14.2f\n", ...
            i, r(ri), v_rel, angle_deg, truth_az);

    detected_targets(i,1) = r(ri);
    detected_targets(i,2) = v_rel;
    detected_targets(i,3) = angle_deg;
end

%% PCM visualization

num_truths = size(truth_targets, 1);

pcm_dots = zeros(num_truths + num_dets, 3);

for i = 1:num_truths
    pcm_dots(i,1) = target_dist(i);
    pcm_dots(i,2) = target_dist(i) * sind(target_az(i));
    pcm_dots(i,3) = truth_targets(i,2);
end

for i = 1:num_dets
    pcm_dots(i+num_truths,1) = detected_targets(i,1);
    pcm_dots(i+num_truths,2) = detected_targets(i,1) * sind(detected_targets(i,3));
    pcm_dots(i+num_truths,3) = detected_targets(i,2);
end

% Swap to [cross-range, range, velocity]
pcm_dots(:,[1,2]) = pcm_dots(:,[2,1]);

fig = figure('Name','PCM Plot','NumberTitle','off','Color','w');
ax  = axes(fig);

hold(ax, 'on');

truthColor = [0 0.4470 0.7410];
detColor   = [0.85 0 0];

truth_x = pcm_dots(1:num_truths,1);
truth_y = pcm_dots(1:num_truths,2);

if num_dets > 0
    detIdx = (num_truths+1):(num_truths+num_dets);
    det_x = pcm_dots(detIdx,1);
    det_y = pcm_dots(detIdx,2);
else
    det_x = [];
    det_y = [];
end

if num_dets > 0
    hDet = plot(ax, det_x, det_y, 'd', ...
        'MarkerSize', 10, ...
        'MarkerFaceColor', detColor, ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 1.2, ...
        'LineStyle', 'none');
else
    hDet = [];
end

hTruth = plot(ax, truth_x, truth_y, 'o', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', truthColor, ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.0, ...
    'LineStyle', 'none');

%% Label offsets that work for any number of detections

base_truth_offsets = [
     0.55, -0.10;
     0.55,  0.25;
     0.45,  0.35;
     0.45,  0.35
];

base_det_offsets = [
    -0.95, -0.85;
    -0.95, -0.85;
    -0.85, -1.00;
    -0.85, -0.90;
     0.65, -1.10;
     0.65,  0.75;
    -1.10,  0.70;
     0.75,  1.00
];

% Truth labels
for i = 1:num_truths

    idx = mod(i-1, size(base_truth_offsets,1)) + 1;

    text(ax, truth_x(i)+base_truth_offsets(idx,1), ...
             truth_y(i)+base_truth_offsets(idx,2), ...
         sprintf('T%d', i), ...
         'FontSize', 9, ...
         'FontWeight', 'bold', ...
         'Color', truthColor);
end

% Detected labels
for i = 1:num_dets

    idx = mod(i-1, size(base_det_offsets,1)) + 1;

    text(ax, det_x(i)+base_det_offsets(idx,1), ...
             det_y(i)+base_det_offsets(idx,2), ...
         sprintf('D%d', i), ...
         'FontSize', 9, ...
         'FontWeight', 'bold', ...
         'Color', detColor);
end

xlabel(ax, 'Cross-range (m)', 'FontWeight','bold');
ylabel(ax, 'Range (m)', 'FontWeight','bold');
title(ax, 'Point Cloud Map — Truth vs Detected', 'FontWeight','bold');

grid(ax, 'on');
axis(ax, 'equal');

set(ax, 'FontSize', 10, 'LineWidth', 0.9);
ax.GridAlpha = 0.22;

if num_dets > 0
    legend(ax, [hTruth, hDet], {'Truth','Detected'}, ...
        'Location','northwest', ...
        'FontSize', 9);
else
    legend(ax, hTruth, {'Truth'}, ...
        'Location','northwest', ...
        'FontSize', 9);
end

x_all = [truth_x(:); det_x(:)];
y_all = [truth_y(:); det_y(:)];

xlim(ax, [min(x_all)-5, max(x_all)+5]);
ylim(ax, [min(y_all)-4, max(y_all)+4]);

hold(ax, 'off');

%% =========================================================
%  Local function: Paper-style improved 2D CFAR
%% =========================================================

function [detects, detMask, Umap, Smap] = cfar_detection_paper(Z, Tscale, window)

    if numel(window) ~= 2
        error('window must be [rows, cols].');
    end

    if mod(window(1),2) == 0 || mod(window(2),2) == 0
        error('window size must be odd.');
    end

    r_offset = (window(1)-1)/2;
    c_offset = (window(2)-1)/2;

    filt = ones(window);
    filt(r_offset+1, c_offset+1) = 0;
    filt = filt / sum(filt(:));

    padded_Z = padarray(Z, [r_offset, c_offset], 'replicate');

    detMask = false(size(Z));
    Umap    = zeros(size(Z));
    Smap    = zeros(size(Z));
    detects = [];

    for r = 1:size(Z,1)
        for c = 1:size(Z,2)

            local_win = padded_Z(r:r+2*r_offset, c:c+2*c_offset);

            U = sum(local_win .* filt, "all");
            S = Tscale * U;

            Umap(r,c) = U;
            Smap(r,c) = S;

            cutVal = padded_Z(r+r_offset, c+c_offset);

            if cutVal > S

                detMask(r,c) = true;
                detects = [detects; r, c]; %#ok<AGROW>

                padded_Z(r+r_offset, c+c_offset) = U;
            end
        end
    end
end