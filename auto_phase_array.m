% This script will borrow the info from https://www.mathworks.com/help/radar/ug/automotive-adaptive-cruise-control-using-fmcw-technology.html
% to create a phase array.
debug_plot = 0;

fc = 77e9;
c = physconst('LightSpeed'); % 299792458 m/s
lambda = c/fc; % meters, wavelength of signal

% TODO: what should our max range and res be?
range_max = 200; % meters, depends on situation. 200m is for cruise control driving situation.
T = 5.5*range2time(range_max,c); % sec, sweep time. scale should be between 5-6x
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
target_dist = [43 50]; % meters         % TODO: what should be targets be?
target_speed = [-80 96]*1000/3600; % meters/sec (note, simulation also sets radar to have a speed)
target_az = [-10 10]; % Azimuth angle in degrees
target_rcs = [20 40]; % radar cross section
target_pos = [target_dist.*cosd(target_az); target_dist.*sind(target_az); 0.5 0.5];
target = phased.RadarTarget('MeanRCS', target_rcs, 'PropagationSpeed', c, 'OperatingFrequency', fc);
target_motion = phased.Platform('InitialPosition', target_pos, 'Velocity', [target_speed; 0 0; 0 0]);

channel = phased.FreeSpace('PropagationSpeed', c, 'OperatingFrequency', fc, 'SampleRate', fs, 'TwoWayPropagation', true);

for i = 1:length(target_dist)
    fprintf(1, "Target %d at range %f, abs vel %f, relative az %f degrees, rcs %f\n", i, target_dist(i), target_speed(i), target_az(i), target_rcs(i));
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
dr = lambda/Nt; % meters, rx spacing

txarray = phased.ULA(Nt, dt);
rxarray = phased.ULA(Nr, dr);
%varray = phased.ULTA(Nt*Nr, dr); % The physical representation of the virtual array

transmitter = phased.Transmitter('PeakPower', tx_ppower, 'Gain', tx_gain);
receiver = phased.ReceiverPreamp('Gain', rx_gain, 'NoiseFigure', rx_nf, 'SampleRate', fs);

txradiator = phased.Radiator('Sensor', txarray, 'OperatingFrequency', fc, 'PropagationSpeed', c, 'WeightsInputPort', true);
rxcollector = phased.Collector('Sensor', rxarray, 'OperatingFrequency', fc, 'PropagationSpeed', c);

radar_speed = 100*1000/3600; % meters/sec, assuming radar is also in motion
radarmotion = phased.Platform('InitialPosition', [0;0;0.5], 'Velocity', [radar_speed;0;0]);

fprintf(1, "Radar with %d Tx, %d Rx, %d Virtual, abs vel %f\n", Nt, Nr, Nt*Nr, radar_speed);
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

    % Toggle transmit element with TDMA
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

xrv = cat(2,xr1, xr2, xr3, xr4); % Xrv size is [num range bins (positive only), num virtual rx, num vel bins];

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
% FIXME: the maximum velocity is too low, and target1's velocity is aliasing from 50m/s to -16m/s
if (debug_plot)
    figure();
    plotResponse(rngdopresp, squeeze(xrv(:,1,:)));
end

%% Sum the magnitude of response
mag_resp = abs(resp(nfft_r/2+1:end, :,:)).^2; % Dump the negative range bins, sum the virtual receivers along dim 2
Z = reshape(sum(mag_resp, 2), [nfft_r/2, nfft_d]);
r = r(nfft_r/2+1:end); % Dumped negative range bins, is now positive range values of bins of Z

% colormap('winter');
% figure(), imagesc(sp, r, Z);
% colorbar;
% title('Magnitude of response');
% xlabel('velocity (m/s)');
% ylabel('range (m)');

%% Detection
%Pfa = 0.001;
window = [5, 3]; % CHECKME: what size? the resolution of resp bins is not the same as our input range_res from the start. (it's closer to half of range_res)
                 % If we window 5 range bins, do we actually have 1 m resolution?
%detects = cfar_detection(Z, Pfa, window); % Deprecated cfar detector function
Tscale = 3.5;
[detects, detMask, Umap, Smap] = cfar_detection_paper(Z, Tscale, window);


%% Determine AOA of detects
num_dets = size(detects, 1);
detected_targets = zeros(num_dets,3);
fprintf(1, "Detected %d range/vel combinations\n", num_dets);
ktheta = [(0:Nt*Nr/2-1) -Nt*Nr/2:-1]/(Nt*Nr)*2*pi; % Angle of each index in Pi
for i = 1:num_dets
    this_range_idx = detects(i,1);
    this_vel_idx = detects(i,2);
    fprintf(1, "Detect %d at range %f and relative vel %f\n", i, r(this_range_idx), sp(this_vel_idx));
    
    % Extract sequence Y
    Y = resp(nfft_r/2 + this_range_idx,:,this_vel_idx); % Range-doppler map values for this detect at each virtual receiver
    
    % Use DFT to determine angle information
    Pi = fft(Y);
    % How to determine the angle?
    % TODO: Make this a peak finder, what if we have two targets at the same range but different az?
    [~, max_angle_idx] = max(Pi);
    fprintf(1, "Found angle of %f\n", 180/pi*ktheta(max_angle_idx))

    detected_targets(i,1) = r(this_range_idx);
    detected_targets(i,2) = sp(this_vel_idx);
    detected_targets(i,3) = 180/pi*ktheta(max_angle_idx);
end

%% Generate PCM using detected_targets and truth_targets
% my_pcm(truth_targets, detected_targets);