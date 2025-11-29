%% Cleanup
clc;
clear all;
close all;

%   DDC PARAMETERS 
Fs = 60e6;
Fc = 15e6;
F0 = Fc;                 % Desired output frequency (Hz)
Deltaf = 0.05;           % Desired frequency resolution (Hz)
SFDR = 90;               % Spurious free dynamic range (dB)
Ts = 1/Fs;               % Sample period (s)
phOffd = pi/2;           % Desired phase offset (rad)

%% NCO Accumulator Size Calculation
N = ceil(log2(1/(Ts * Deltaf)));
Q = ceil((SFDR - 12)/6);
phIncr = round(F0 * 2^N * Ts);
phOff = (2^N * phOffd)/(2*pi);
ditherBits = N - Q;

%%   FIR LOW-PASS FILTER 
L = 15;  
Wn = Fc/(Fs/2);     % normalized cutoff (0.5)
win = blackman(L+1);
h = fir1(L, Wn, "low", win, "scale");
hf = fi(h, 1, 12, 11);      % fixed point representation

figure;
stem(h);
title("FIR Low-Pass Filter Coefficients");

figure;
freqz(h);
title("Frequency Response of Low-Pass Filter");


%%   RUN SIMULINK MODEL
% sim('YourDDC_SimulinkModelName');  

%   EXTRACT FFT DATA FROM SIMULINK
FFT = 512;     % FFT size used in Simulink

%% Find the first index where valid == 1
k = find(out.SciInfo_Valid, 1, 'first');

% Extract magnitude starting at that valid index
fftData = out.SciInfo(k:end);

%%   PLOT FIRST (INSTANTANEOUS) FFT FRAME
a = (1:FFT)';

figure;
subplot(2,1,1);
plot(a, fftData(1:FFT), 'k', 'LineWidth', 1);
title("Instantaneous Spectrum");
xlabel("FFT Bin");
ylabel("Magnitude");
xlim([0 FFT]);

%%   SPECTRAL AVERAGING
specAvg = floor(length(fftData)/FFT);

% Keep only full frames
fftData_snr1 = fftData(1:specAvg * FFT);

% Reshape into a matrix of size 512 x (numAvg)
fftData_snr2 = reshape(fftData_snr1, FFT, specAvg);

% Sum along columns
fftData_snr3 = sum(fftData_snr2, 2);

subplot(2,1,2)
plot(a, fftData_snr3, 'r');
xlim([0 FFT])
