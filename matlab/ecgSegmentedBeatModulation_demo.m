%%
% ECGSEGMENTEDBEATMODULATION_DEMO - Demo for the function
% ecgSegmentedBeatModulation.
%
% Created July 16, 2018.
% Arturo Moncada-Torres
%   arturomoncadatorres@gmail.com
%   http://www.arturomoncadatorres.com


%% Preliminaries.
clc;
clear variables;
close all;


%% Define execution parameters.

% Obtained from the corresponding .info file.
fs = 500;           % Sampling rate [Hz].
Ts = 1/fs;          % Sampling period [s].
gain = 100;         % Gain [u.a.]

deltaT = 40e-3;             % Suggested value is 40e-3 (p. 50) [s]
deltaT_N = deltaT * fs;     % [samples]


%% Load data.

% Data is contained in "ecg.mat".
% The original data was obtained from the MIT-BIH Polysomnographic
% Database (http://www.physionet.org/cgi-bin/atm/ATM).
% - Database: Motion Artifact Contaminated ECG Database (macegdb)
% - Record: test05_45s
% - Signals: all
% - Length: to end
% - Time format: samples
% - Data format: standard
tmp = load(fullfile('..','data','ecg'));
ecg = tmp.val;

% Each row corresponds to different signals.
% For now, we will work with ECG3.
ecg_raw = ecg(3,:);


%% Preprocess data.

% Convert to mV.
ecg_mV = ecg_raw / gain;    % [mV]

% Filter ECG (as suggested in the paper).

% Filter parameters.
n = 8;      % Desired order.
fc1 = 0.5;  % Cut-off frequency 1 (suggested value is 0.5) [Hz]
fc2 = 45;   % Cut-off frequency 2 (suggested value is 45) [Hz]

% We divide the desired filter order by 4 because creating a bandpass
% filter doubles the order, as well as the FILTFILT function.
[num, den] = butter(round(n/4),[fc1, fc2]/fs);

% Actual filtering.
ecg_pre  = filtfilt(num, den, ecg_mV);  % [mV]


%% Find R peaks.
% Here, we will use a simple thresholding operation.
% However, there are better (i.e., more robust) ways to do so.

minPeakHeight = 0.7;
[~, RPeaks] = findpeaks(ecg_pre,'MinPeakHeight',minPeakHeight, ...
                'MinPeakDistance',deltaT_N);     % [samples]


%% Segmented beat modulation processing.
[ecg_clean, mCC, modulatedCC] = ecgSegmentedBeatModulation(ecg_pre, fs, RPeaks, deltaT);


%% Plots.
blue = [30,150,255]/255;
gray = [0.66,0.66,0.66];

% ECG signal.
t = 0:Ts:length(ecg_mV)*Ts-Ts;

figure('Name','Pre-processed ECG Signal');
hold('on');
plot(t, ecg_pre, 'LineWidth',2, 'Color',blue);
plot(t, ones(size(ecg_pre))*minPeakHeight, 'LineStyle','--', 'Color',gray);
for ii = 1:length(RPeaks)
    plot([RPeaks(ii),RPeaks(ii)].*Ts,[-100,100], 'LineStyle','--', 'Color',gray);
end
hold('off');
ylim([-0.75,1.25]);
title('Pre-processed ECG Signal');
xlabel('Time [s]');
ylabel('Amplitude [mV]')


% mCC Computation.
% Equivalent of Fig. 3 of the paper.
t_mCC = 0:Ts:length(mCC)*Ts-Ts;

figure('Name','mCC Computation');
hold('on');
plot(t_mCC,modulatedCC', 'Color',gray);
plot(t_mCC,mCC, 'LineWidth',2, 'Color',blue);
hold('off');
ylim([-0.75,1.25]);
title('mCC Computation');
xlabel('Time [s]');
ylabel('Amplitude [mV]')


% Signal and motion artifact.
% Equivalent of Fig. 5 of the paper.

% Trim signal.
ecg_pre2 = ecg_pre(RPeaks(1)-deltaT_N:RPeaks(1)-deltaT_N+length(ecg_clean)-1);
artifact = ecg_pre2 - ecg_clean;

t = 0:Ts:length(artifact)*Ts-Ts;

figure('Name','ECG and Artifact');

subplot(2,1,1);
plot(t,ecg_pre2, 'LineWidth',2, 'Color',blue);
ylim([-0.75,1.25]);
title('Pre-processed ECG');
ylabel('Amplitude [mV]')
subplot(2,1,2);
hold('on');
plot(t,artifact, 'LineWidth',2, 'Color',gray);
plot(t,ecg_clean, 'LineWidth',2, 'Color',blue);
hold('off');
ylim([-0.75,1.25]);
title('Clean ECG and Artifact');
xlabel('Time [s]');
ylabel('Amplitude [mV]')
legend({'Motion Artifact','Clean ECG'});
legend boxoff;


% Signal comparison.
t = 0:Ts:length(ecg_mV)*Ts-Ts;

% To make a better comparison, we align the clean ECG to the
% original ECG with NaN.
% We could, of course, with further processing align the clean ECG with
% actual filtered ECG data.
ecg_clean2 = nan(size(ecg_mV));
ecg_clean2(RPeaks(1)-deltaT_N:RPeaks(1)-deltaT_N+length(ecg_clean)-1) = ecg_clean;

figure('Name','Signal Comparison');
subplot(3,1,1);
plot(t,ecg_mV,'Color',blue);
ylim([-0.75,1.25]);
title('Original ECG');
ylabel('Amplitude [mV]')
subplot(3,1,2);
plot(t,ecg_pre,'Color',blue);
ylim([-0.75,1.25]);
title('Preprocessed ECG');
ylabel('Amplitude [mV]')
subplot(3,1,3);
plot(t,ecg_clean2,'Color',blue);
ylim([-0.75,1.25]);
title('Segmented-Beat Modulation Filtered ECG');
xlabel('Time [s]');
ylabel('Amplitude [mV]')
