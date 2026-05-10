% This script will borrow the info from https://www.mathworks.com/help/radar/ug/automotive-adaptive-cruise-control-using-fmcw-technology.html
% to create a phase array.
debug_plot = 0;

fc = 77e9;
c = physconst('LightSpeed'); % 299792458 m/s
lambda = c/fc; % meters, wavelength of signal

% TODO: what should our max range and res be?
range_max = 200; % meters, depends on situation. 200m is for cruise control driving situation.
T = 5*range2time(range_max,c); % sec, sweep time. scale should be between 5-6x
% T = 5.5*range_max*c/2;

range_res = 1; % meters, desired range resolution. 1 for cruise control driving situation.
bw = rangeres2bw(range_res, c); % Hz, bandwidth of pulse
%bw = c/(range_res*2);
S = bw/T; % Hz/sec, sweep rate or chirp rate

fr_max = range2beat(range_max, S, c);
%fr_max = 2*S*range_max/c;

v_max = 230*1000/3600; % meters/sec, expected maximum speed of objects. 230km/h for cruise control driving situation
fd_max = speed2dop(2*v_max,lambda); % Hz/sec, maximum Dopper shift
fb_max = fr_max+fd_max; % Hz/sec, maximum beat frequency

fs = max(2*fb_max, bw); % Hz, sampling rate

%% Define FMCW pulse
waveform = phased.FMCWWaveform('SweepTime', T, 'SweepBandwidth', bw, 'SampleRate', fs);
sig = waveform();
if (debug_plot)
    subplot(211); plot(0:1/fs:T-1/fs,real(sig));
    xlabel('Time (s)'); ylabel('Amplitude (v)');
    title('FMCW signal'); axis tight;
    subplot(212); spectrogram(sig,32,16,32,fs,'yaxis');
    title('FMCW signal spectrogram');
end

%% Setup target
% target_dist = [20 50]; % meters         % TODO: what should be targets be?
% target_speed = [99 96]*1000/3600; % meters/sec (note, simulation also sets radar to have a speed)
% target_az = [-15.029342 63.2709384]; % Azimuth angle in degrees
% target_rcs = [20 40]; % radar cross section

%% Final target setup?
% 4 targets. 1 pair to check range resolution, 1 pair to check angle resolution
% How difficult is this?
target_dist = [44 42 31 60];
target_speed = [105.4724 90.24 109.84 92.0190]*1000/3600;
target_az = [-25.35 -26.12 13.203 16.640];
target_rcs = [20 20 20 20];

target_pos = [target_dist.*cosd(target_az); target_dist.*sind(target_az); 0.5*ones(1,length(target_dist))];
target = phased.RadarTarget('MeanRCS', target_rcs, 'PropagationSpeed', c, 'OperatingFrequency', fc);
target_motion = phased.Platform('InitialPosition', target_pos, 'Velocity', [target_speed; zeros(1, length(target_dist)); zeros(1, length(target_dist))]);

channel = phased.FreeSpace('PropagationSpeed', c, 'OperatingFrequency', fc, 'SampleRate', fs, 'TwoWayPropagation', true);

for i = 1:length(target_dist)
    fprintf(1, "Target %d at range %f m, abs vel %f m/s, relative az %f degrees, rcs %f\n", ...
        i, target_dist(i), target_speed(i)*3.6, target_az(i), target_rcs(i));
end
%% Setup Radar System
ant_aperture = 6.06e-4; % square meters
ant_gain = aperture2gain(ant_aperture, lambda); % dB
tx_ppower = db2pow(5)*1e-3; % watts
tx_gain = 9+ant_gain; % dB
rx_gain = 15+ant_gain; % dB
rx_nf = 4.5; % dB

Nt = 4; % num tx        % TODO: should we implement the 9x16 array in the paper?
Nr = 8; % num rx

dt = Nr*lambda/2; % meters, tx spacing (half wavelength)
dr = lambda/2; % meters, rx spacing

txarray = phased.ULA(Nt, dt);
rxarray = phased.ULA(Nr, dr);
varray = phased.ULA(Nt*Nr, dr); % The physical representation of the virtual array
% viewArray(txarray,'ShowTaper',true,'ShowIndex','All');
% viewArray(rxarray,'ShowTaper',true,'ShowIndex','All');
% viewArray(varray, 'ShowTaper',true,'ShowIndex','All');

transmitter = phased.Transmitter('PeakPower', tx_ppower, 'Gain', tx_gain);
receiver = phased.ReceiverPreamp('Gain', rx_gain, 'NoiseFigure', rx_nf, 'SampleRate', fs);

txradiator = phased.Radiator('Sensor', txarray, 'OperatingFrequency', fc, 'PropagationSpeed', c, 'WeightsInputPort', true);
rxcollector = phased.Collector('Sensor', rxarray, 'OperatingFrequency', fc, 'PropagationSpeed', c);

radar_speed = 100*1000/3600; % meters/sec, assuming radar is also in motion
radar_speed_msec = radar_speed*3.6;
radarmotion = phased.Platform('InitialPosition', [0;0;0.5], 'Velocity', [radar_speed;0;0]);

fprintf(1, "Radar with %d Tx, %d Rx, %d Virtual, abs vel %f m/s\n", Nt, Nr, Nt*Nr, radar_speed_msec);
% Generate "Truth" target matrix
truth_targets = [target_dist.' (radar_speed - target_speed).' target_az.'];

%% Simulation
if (debug_plot)
    specanalyzer = spectrumAnalyzer('SampleRate',fs, ...
    'Method','welch','AveragingMethod','running', ...
    'PlotAsTwoSidedSpectrum',true, 'FrequencyResolutionMethod','rbw', ...
    'Title','Spectrum for received and dechirped signal', ...
    'ShowLegend',true);
end
rng(11);
Dn = 2; % Decimation factor
fs = fs/Dn;
Nsweep = 64/2*Nt;
xr = complex(zeros(fs*waveform.SweepTime,Nr, Nsweep));

w0 = zeros(Nt,1); % weights to enable/disable radiating elements
w0(Nt) = 1;

for m = 1:Nsweep
    % Update radar and target positions
    [radar_pos, radar_vel] = radarmotion(waveform.SweepTime);
    [target_pos, target_vel] = target_motion(waveform.SweepTime);
    [~, target_ang] = rangeangle(target_pos, radar_pos);

    % Transmit FMCW waveform
    sig = waveform();
    txsig = transmitter(sig);

    % Toggle transmit element
    w0 = [w0(Nt); w0(1:Nt-1)];
    txsig = txradiator(txsig, target_ang, w0);

    % Propagate the signal and reflect off the target
    txsig = channel(txsig, radar_pos, target_pos, radar_vel, target_vel);
    txsig = target(txsig);

    % Dechirp the received radar return
    rxsig = rxcollector(txsig, target_ang);
    rxsig = receiver(rxsig);
    dechirpsig = dechirp(rxsig, sig);

    % Decimate the return to reduce computation requirements
    for n = size(xr,2):-1:1
        xr(:,n,m) = decimate(dechirpsig(:,n),Dn, 'FIR');
    end

    % Visualize the spectrum
    if (debug_plot)
        specanalyzer([txsig dechirpsig]);
    end

end

%% Virtual Array processing
xr1 = xr(:,:,1:Nt:end); % taking every other page to recover the measurements corresponding to the two transmit antenna elements
xr2 = xr(:,:,2:Nt:end); % When Nt > 2, we still toggle 1 receiver per sweep
xr3 = xr(:,:,3:Nt:end); % Note: Need to go up to xr<Nt>
xr4 = xr(:,:,4:Nt:end);
% xr5 = xr(:,:,5:Nt:end);
% xr6 = xr(:,:,6:Nt:end);
% xr7 = xr(:,:,7:Nt:end);
% xr8 = xr(:,:,8:Nt:end);
% xr9 = xr(:,:,9:Nt:end);

%xrv = cat(2,xr1, xr2);
xrv = cat(2,xr1, xr2, xr3, xr4); % Xrv size is [num range bins (positive only), num virtual rx, num vel bins];
%xrv = cat(2,xr1, xr2, xr3, xr4, xr5, xr6); %xr7, xr8, xr9);
%% Range and Doppler Estimation
nfft_r = 2^nextpow2(size(xrv,1));
nfft_d = 2^nextpow2(size(xrv,3));

rngdopresp = phased.RangeDopplerResponse('PropagationSpeed',c,...
    'DopplerOutput','Speed','OperatingFrequency',fc,'SampleRate',fs,...
    'RangeMethod','FFT','PRFSource', 'Property', ...
    'RangeWindow', 'Hann', 'PRF', 1/(Nt*waveform.SweepTime), ...
    'SweepSlope', S,...
    'RangeFFTLengthSource','Property','RangeFFTLength',nfft_r,...
    'DopplerFFTLengthSource','Property','DopplerFFTLength',nfft_d, ...
    'DopplerWindow', 'Hann');

[resp, r, sp] = rngdopresp(xrv);    % resp
                                    % r: range values of resp bins dim 1
                                    % sp: velocity values of resp bins dim 3

if (debug_plot)
    figure();
    plotResponse(rngdopresp, squeeze(xrv(:,1,:)));
end

%% Sum the magnitude of response
mag_resp = abs(resp(nfft_r/2+1:end, :,:)).^2; % Dump the negative range bins, sum the virtual receivers along dim 2
Z = reshape(sum(mag_resp, 2), [nfft_r/2, nfft_d]);
r = r(nfft_r/2+1:end); % Dumped negative range bins, is now positive range values of bins of Z

if (debug_plot)
imagesc(sp, r, Umap);
colorbar;
title('BNE Umap');
xlabel('velocity (m/s)');
ylabel('range (m)');
end

%% Detection
%Pfa = 0.001;
window = [5, 3]; % CHECKME: what size? the resolution of resp bins is not the same as our input range_res from the start. (it's closer to half of range_res)
                 % If we window 5 range bins, do we actually have 1 m resolution?
%detects = cfar_detection(Z, Pfa, window); % Deprecated cfar detector function
Tscale = 4.0;
[detects, detMask, Umap, Smap] = cfar_detection_paper(Z, Tscale, window);


%% Determine AOA of detects
Nrx = Nt * Nr; % Total number of virtual channels
num_dets = size(detects, 1);
detected_targets = zeros(num_dets,3);
fprintf(1, "Detected %d range/vel combinations\n", num_dets);
% FIXME: I assumed the Angle range was -PI to PI, but the actual range might be smaller. How do we know?
ktheta = [(0:Nt*Nr/2-1) -Nt*Nr/2:-1]/(Nt*Nr)*2*pi; % Angle of each index in Pi

musicEstimator = phased.MUSICEstimator('SensorArray', varray, ...
                    OperatingFrequency=fc,...
                    DOAOutputPort=true,NumSignalsSource="Property",...
                    NumSignals=num_dets); % NOTE: If we trust our CFAR to only find unique targets, we could use detected_targets here.
[spectrum, doa] = musicEstimator(xrv);
doa = broadside2az(sort(doa)) % zero elevation
% this function seems to have the opposite sign for az as the radar equations, so negate the doas
doa = -1*doa;
plotSpectrum(musicEstimator,NormalizeResponse=true)

for i = 1:num_dets
    this_range_idx = detects(i,1);
    this_vel_idx = detects(i,2);

    % Extract sequence Y
    Y = resp(nfft_r/2 + this_range_idx,:,this_vel_idx); % Range-doppler map values for this detect at each virtual receiver

    % Use DFT to determine angle information
    Pi = fft(Y, Nrx);
    [~, max_angle_idx] = max(abs(Pi));

    if max_angle_idx > Nrx/2
        k_norm = (max_angle_idx - 1 - Nrx) / Nrx;
    else
        k_norm = (max_angle_idx - 1) / Nrx;
    end

    %Flip the sign
    k_norm = -k_norm;

    % Calculate Angle using the Inverse Sine
    % Formula: theta = asin( f_theta * lambda / dr )
    % Since k_norm represents (d/lambda * sin(theta)), we divide by (d/lambda)
    val_to_asin = k_norm / (dr / lambda);

    val_to_asin = max(min(val_to_asin, 1), -1); 

    angle_rad = asin(val_to_asin);
    angle_deg = rad2deg(angle_rad);

    % Use the results of MUSIC Estimator to find closest
    [~, doa_index] = min(abs(angle_deg - doa));
    
    v_rel = sp(this_vel_idx);
    v_target_abs = v_rel + radar_speed_msec;
    fprintf(1, "Detect %d at range %f and rel vel %f m/s.", i, r(this_range_idx), v_rel);
    

    fprintf(1, " Found angle of %f deg. MUSIC found %f deg.\n", angle_deg, doa(doa_index));
    detected_targets(i,1) = r(this_range_idx);
    detected_targets(i,2) = v_rel;
    detected_targets(i,3) = angle_deg;
end

% --- PCM Visualization ---
num_truths = size(truth_targets, 1);
% [Range, CrossRange, Rel Vel]
pcm_dots = zeros(num_truths + num_dets, 2);

for i = 1:num_truths
    pcm_dots(i, 1) = target_dist(i);
    pcm_dots(i, 2) = target_dist(i) * sin(deg2rad(target_az(i)));
    pcm_dots(i, 3) = target_speed(i)*3.6 - radar_speed_msec;
end

for i = 1:num_dets
    pcm_dots(i + num_truths, 1) = detected_targets(i, 1);
    pcm_dots(i + num_truths, 2) = detected_targets(i, 1) * sin(deg2rad(detected_targets(i,3)));
    pcm_dots(i + num_truths, 3) = detected_targets(i, 2);
end

%Swap X and Y
% [Cross Range, Range, Rel Vel]
pcm_dots(:, [1, 2]) = pcm_dots(:, [2, 1]);

% Prepare figure
fig = figure('Name','PCM Plot','NumberTitle','off');
ax = axes(fig); 
axis(ax, 'equal');
hold(ax, 'on');

% Colors and marker styles
truthColor = [0 0.4470 0.7410]; 
detColor = [.85 0 0];

% Initialize legend containers
legHandles = []; % Using an empty array for easier concatenation
legEntries = {};

% --- TRUTH DOTS & ARROWS ---
numTruthsToPlot = min(num_truths, size(pcm_dots,1));
detIdx = 1:numTruthsToPlot; % Use numTruthsToPlot for safety

dots_x = pcm_dots(detIdx, 1);
dots_y = pcm_dots(detIdx, 2);
v_rels = pcm_dots(detIdx, 3);

% 1. Calculate Truth Arrows (pointing toward/away from origin)
u_unit = dots_x;
v_unit = dots_y;
magnitudes = sqrt(u_unit.^2 + v_unit.^2);
magnitudes(magnitudes == 0) = 1; 

u_unit = u_unit ./ magnitudes;
v_unit = v_unit ./ magnitudes;
scale_factor = 0.5;
u_vel = u_unit .* v_rels * scale_factor;
v_vel = v_unit .* v_rels * scale_factor;

hArrow2 = quiver(ax, dots_x, dots_y, u_vel, v_vel, 0, ...
       'MaxHeadSize', 0.5, 'Color', truthColor, 'LineWidth', 1.2);


h = plot(ax, dots_x, dots_y, 'Marker', 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', truthColor, 'LineStyle', 'none', 'Color', truthColor);

% Add Text Labels for Truth
for k = 1:length(dots_x)
    tip_x = dots_x(k) + u_vel(k);
    tip_y = dots_y(k) + v_vel(k);
    label = sprintf('%.1f m/s', abs(v_rels(k)));
    text(ax, tip_x, tip_y, label, 'FontSize', 9, 'FontWeight', 'bold');
end

legHandles(end+1) = h; 
legEntries{end+1} = 'Truth'; 

% --- DETECTED DOTS & ARROWS ---
if num_truths < size(pcm_dots,1)
    detIdx = (num_truths+1):size(pcm_dots,1);
    
    dots_x = pcm_dots(detIdx, 1);
    dots_y = pcm_dots(detIdx, 2);
    v_rels = pcm_dots(detIdx, 3);
    
    u_unit = dots_x;
    v_unit = dots_y;
    magnitudes = sqrt(u_unit.^2 + v_unit.^2);
    magnitudes(magnitudes == 0) = 1; 
    
    u_unit = u_unit ./ magnitudes;
    v_unit = v_unit ./ magnitudes;
    u_vel = u_unit .* -v_rels * scale_factor;
    v_vel = v_unit .* -v_rels * scale_factor;
    
    % Draw Detected arrows and capture handle for legend
    hArrow = quiver(ax, dots_x, dots_y, u_vel, v_vel, 0, ...
           'MaxHeadSize', 0.5, 'Color', [1 0 0], 'LineWidth', 1.5);
    
    % Draw Detected dots
    hDet = plot(ax, dots_x, dots_y, 'd', 'MarkerSize', 10, ...
        'MarkerFaceColor', detColor, 'Color', 'k', 'LineWidth', 1.5);
    
    legHandles(end+1) = hDet;
    legEntries{end+1} = 'Detected';
    
    % Add Text Labels for Detected
    for k = 1:length(dots_x)
        tip_x = dots_x(k) + u_vel(k);
        tip_y = dots_y(k) + v_vel(k);
        label = sprintf('%.1f m/s', abs(v_rels(k)));
        text(ax, tip_x, tip_y, label, 'FontSize', 9, 'FontWeight', 'bold');
    end
end

% Finalize plot
xlabel(ax,'Cross-range (meters)');
ylabel(ax,'Range (meters)');
grid(ax, 'on');
title(ax,'Truth and Detected Dots');

valid = isgraphics(legHandles);
if any(valid)
    legend(ax, legHandles(valid), legEntries(valid), 'Location','best');
end
hold(ax, 'off');