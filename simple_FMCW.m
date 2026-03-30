% Using variable names from the fundamentals of mmwave ... paper
c = physconst('LightSpeed'); % 299792458 m/s
noise_std = 3;

hann_L = 256;
hann = 1/hann_L*cos(pi*([0:hann_L-1] - hann_L/2)/hann_L).^2;

Fs = 10e3; % MHz
t = 0:1/Fs:40-1/Fs; % usec

%% Synthesize FMCW truth pulse
fc = 1e3; % MHz, start frequency of synthesized pulse
S = 100; % MHz/usec, chirp rate of pulse (slope)
Tc = 20; % us, width of pulse
B = Tc*S; % MHz, bandwidth of pulse
SNR_truth = 15; % dB, SNR of cosine is 10*log10(A^2/(2*noise_std^2));
A = sqrt(10^(0.1*SNR_truth)*2*noise_std^2); % amplitude of pulse
phi0 = pi/10; % Initial phase of pulse

synth = zeros(size(t));
synth(1:floor(Tc*Fs)) = A*cos(2*pi*(fc*t(1:floor(Tc*Fs)) + S/2*t(1:floor(Tc*Fs)).^2) + phi0);
%spectrogram(synth, hann);

max_d = Tc*1e-6*c/2; % meters, range of truth pulse
d_res = c/(2*B*1e6); % meters, range resolution

%% Simulate echo from object 1
d1 = 1400; % meters, Distance to object 1 (must be < Tc e-6*c/2)
if (d1 >= max_d)
    warning('Distance [%f] to object 1 is too large > [%f]', d1, max_d);
end
tau1 = 2*d1/c*1e6; % usec, time offset from tx start to rx start (must be < Tc)
SNR_1 = 15; % dB
Ar1 = sqrt(10^(0.1*SNR_1)*2*noise_std^2); % echo signal amplitude
t_ind_1 = ceil(tau1*Fs):floor(Fs*(tau1+Tc));
echo1 = zeros(size(t));
echo1(t_ind_1) = Ar1*cos(2*pi*(fc*(t(t_ind_1)-tau1) + S/2*(t(t_ind_1)-tau1).^2) + phi0);

w = noise_std*randn(size(t));

rx = echo1+w;

mix = synth.*rx; % mixed IF from echo
%mix = conv(hann, mix);
%mix(1:hann_L/2) = [];
%mix = mix(1:end-hann_L/2+1);

%% Predict object 1 IF
lambda = c/(fc*1e6); % meters, wavelength of synth signal
phif_1 = 4*pi*d1/lambda; % radians, expected initial phase of IF signal
f0_1 = S*1e6*2*d1/c; % MHz, predicted intermediate frequency
predicted_IF = zeros(size(t));
predicted_IF(ceil(tau1*Fs):floor(Tc*Fs)) = A*Ar1*sin(2*pi*f0_1*t(ceil(tau1*Fs):floor(Tc*Fs)) + phif_1 + phi0);