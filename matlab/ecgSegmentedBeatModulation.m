function [ecg_clean, varargout] = ecgSegmentedBeatModulation(ecg, fs, RPeaks, varargin)
%%
% ECGSEGMENTEDBEATMODULATION - ECG artifact removal using
%   segmented-beatmodulation.
%
% ECG_CLEAN = ECGSEGMENTEDBEATMODULATION(ECG, FS, RPEAKS);
%   Remove artifacts of the ECG signal with sampling rate FS and 
%   RPEAKS using the Segmented-beat Modulation algorithm,
%   yielding ECG_CLEAN.
%
% [...] = ECGSEGMENTEDBEATMODULATION(..., DELTAT);
%   Use a different DELTAT value than the one originally proposed
%   (40e-3 s = 40 ms)
%
% [..., MCC, MODULATEDCC] = ECGSEGMENTEDBEATMODULATION(...);
%   Output the median of the cardiac cyle MCC and the modulated
%   cardiac cycles (MODULATEDCC).
%
%   Inputs
%       ECG             Original ECG signal [uV].
%       FS              Sampling rate [Hz]
%       RPEAKS          Location of R peaks [samples]
%       DELTAT          (Optional) Half duration of the QRS segment.
%                       Default value is 40e-3 [s]
%
%   Outputs
%       ECG_CLEAN       Clean ECG signal [uV]
%       MCC             Median of the cardiac cycle [uV]
%       MODULATEDCC     Modulated cardiac cycles [uV]
%
%   This is the MATLAB implementation of the algorithm proposed in
%   the paper
%       Agostinelli, Angela, Corrado Giuliani, and Laura Burattini.
%       "Extracting a clean ECG from a noisy recording: a new method
%       based on segmented-beat modulation." Computing in Cardiology
%       Conference (CinC), 2014. IEEE, 2014.
%       https://ieeexplore.ieee.org/abstract/document/7042976
%   Full credit goes to the authors of the paper.
%
% Created July 16, 2018.
% Arturo Moncada-Torres
%   arturomoncadatorres@gmail.com
%   http://www.arturomoncadatorres.com


%% Manage optional inputs.
nVarArgs = length(varargin);

% Set defaults.
optArgs = {40e-3};  % deltaT. Preferred value is 40e-3 (p. 50) [s]

% Overwrite arguments specified in varargin.
optArgs(1:nVarArgs) = varargin;

% Empty values in memorable variables.
[deltaT] = optArgs{:};


%% Important parameters.
deltaT_N = deltaT * fs;     % [samples]


%% Median cardiac cycle duration (mCCd) computation
% Left side of the block diagram of Fig. 2.
nBeats = length(RPeaks);
            
% Compute mCCd.
% mCCd is the median duration of the cardiac cycles.
% It is computed from the R-R intervals.
RRIntervals = diff(RPeaks);
mCCd = round(median(RRIntervals));     % [samples]


%% Median cardiac cycle (mCC) computation
% Right side of the block diagram of Fig. 2.

% (Beginning of each) cardiac cycle identification.
cc = RPeaks - deltaT_N;

% Cardiac cycle segmentation.
% Notice that we won't take into account the last cardiac cycle, since it
% might very well be incomplete.
% With further processing, the last cycle could be used, filling the
% missing data with NaN.
% The more cardiac cycles that are included, the smoother the ECG
% signal will result.
for ii = 1:nBeats-1
    CCSegmentation{ii} = ecg(cc(ii):cc(ii+1)-1); %#ok
end

% QRS and TUP segments.
for ii = 1:nBeats-1
    QRSSegments{ii} = CCSegmentation{ii}(1:2*deltaT_N); %#ok
    TUPSegments{ii} = CCSegmentation{ii}((2*deltaT_N)+1:end); %#ok
end

% Modulate TUP segments.
mTUPd = mCCd - 2*deltaT_N;  % Median length of TUP segment [samples]
for ii = 1:nBeats-1
    TUPSegmentsModulated{ii} = resample(TUPSegments{ii},mTUPd,length(TUPSegments{ii})); %#ok
end

% Cardiac cycle construction.
modulatedCC = zeros(nBeats-1,mCCd);
for ii = 1:nBeats-1
    modulatedCC(ii,:) = [QRSSegments{ii} TUPSegmentsModulated{ii}];
end

% Compute mCC (median of the cardiac cycle).
% The original paper uses the mean.
% If wished, this could be replaced by the mean.
mCC = median(modulatedCC,1);


%% Clean ECG extraction.
% Fig. 4.

% mCC segmentation.
mQRSSegment = mCC(1:2*deltaT_N);
mTUPSegment = mCC((2*deltaT_N)+1:end);

% Demodulate mTUP segment.
for ii = 1:nBeats-1
    TUPSegmentsDemodulated{ii} = resample(mTUPSegment,length(TUPSegments{ii}),mTUPd); %#ok
end

% Cardiac cycles construction.
for ii = 1:nBeats-1
    demodulatedCC{ii} = [mQRSSegment TUPSegmentsDemodulated{ii}]; %#ok
end

% Cardiac cycle concatenation.
ecg_clean = [];
for ii = 1:nBeats-1
    ecg_clean = [ecg_clean demodulatedCC{ii}]; %#ok
end


%% Pack (optionl) outputs.
varargout{1} = mCC;
varargout{2} = modulatedCC;