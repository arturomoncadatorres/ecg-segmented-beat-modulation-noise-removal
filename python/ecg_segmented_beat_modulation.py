# -*- coding: utf-8 -*-
"""
ecg_segmented_beat_modulation.py
Implementation of the ECG noise removal algorithm using segmented-beat 
modulation proposed by Agostinelli et al., (2014)

Created on Sun Jun  2 16:23:05 2019
@author: Arturo Moncada-Torres
arturomoncadatorres@gmail.com
"""


#%% Preliminaries
import pathlib
import numpy as np
from scipy import signal, io
import matplotlib as mpl
from matplotlib import pyplot as plt


# Define paths.
PATH_DATA = pathlib.Path(r'../data/')
PATH_IMAGES = pathlib.Path(r'../images/')

# Make sure directories exist
if not PATH_DATA.exists():
    raise Exception("Data directory does not exist.")
if not PATH_IMAGES.exists():
    PATH_IMAGES.mkdir()


#%%
def ecg_segmented_beat_modulation_noise_removal(ecg, fs, r_peaks, delta_t=40e-3):
    """
    ECG artifact removal using segmented-beat modulation.
    Full credit goes to the authors of the paper.
    
    Parameters
    ----------
    ecg: type
        Original (pre-processed) ECG signal [uV]
        
    fs:
        Sampling rate [Hz]
        
    r_peaks:
        Location of R peaks [samples]
        
    delta_t:
        Half duration of the QRS segment [s]
        Default value is 40e-3
    
    Returns
    -------
    ecg_clean:
        Clean ECG signal [uV]
        
    mCC:
        Median of the cardiac cycle [uV]
        
    cc_modulated:
        Modulated cardiac cycles [uV]
        
    References
    ----------
    Agostinelli, Angela, Corrado Giuliani, and Laura Burattini. "Extracting a 
    clean ECG from a noisy recording: a new method based on segmented-beat 
    modulation." Computing in Cardiology Conference (CinC), 2014. IEEE, 2014.
    https://ieeexplore.ieee.org/abstract/document/7042976
    """

    # Important parameters.
    delta_t_n = delta_t * fs
    
    
    #%% Median cardiac cycle duration (mCCd) computation.
    
    # Left side of the block diagram of Fig. 2.
    n_beats = len(r_peaks)
                
    # Compute mCCd.
    # mCCd is the median duration of the cardiac cycles.
    # It is computed from the R-R intervals.
    rr_intervals = np.diff(r_peaks)
    mCCd = round(np.median(rr_intervals)) # [samples]
    
    # Compute the median duration of TUP cycles (mTUPd).
    mTUPd = mCCd - 2*delta_t_n # Median length of TUP segment [samples]

    
    #%% Median cardiac cycle (mCC) computation
    # Right side of the block diagram of Fig. 2.
    # (Beginning of each) cardiac cycle identification.
    cc = r_peaks - delta_t_n
    
    # Cast to int.
    cc = cc.astype(int)
    
    # Cardiac cycle segmentation.
    # Notice that we won't take into account the last cardiac cycle, since it
    # might very well be incomplete.
    # With further processing, the last cycle could be used, filling the
    # missing data with NaN (for instance).
    # The more cardiac cycles that are included, the smoother the ECG
    # signal will result.
    cc_segmentation = list()
    qrs_segments = list()
    tup_segments = list()
    tup_segments_modulated = list()
    cc_modulated = list()
    for ii in range(0, n_beats-1):
        cc_segmentation.append(ecg[cc[ii]:cc[ii+1]])
        qrs_segments.append(cc_segmentation[ii][0:int(delta_t_n*2)])
        tup_segments.append(cc_segmentation[ii][int(delta_t_n*2):])

        # Modulate TUP segments.
        tup_segments_modulated.append(signal.resample(tup_segments[ii], int(mTUPd)))
        
        # Cardiac cycle reconstruction.
        cc_modulated.append(np.concatenate((qrs_segments[ii], tup_segments_modulated[ii])))
        

    # Compute mCC (median of the cardiac cycle).
    # The original paper uses the median.
    # If wished, this could be replaced by the mean.
    mCC = np.median(cc_modulated, 0)
     
    #%% Clean ECG extraction (Fig. 4)
    
    # mCC segmentation.
    qrs_segment_median = mCC[0:int(delta_t_n*2)]
    tup_segment_median = mCC[int(delta_t_n*2):]
    
    # Demodulate tup_segment_median.
    tup_segments_demodulated = list()
    cc_demodulated = list()
    for ii in range(0, n_beats-1):
        tup_segments_demodulated.append(signal.resample(tup_segment_median, len(tup_segments[ii])))
        cc_demodulated.append(np.concatenate((qrs_segment_median, tup_segments_demodulated[ii])))


    # Cardiac cycle concatenation.
    ecg_clean = [item for sublist in cc_demodulated for item in sublist]

    return (ecg_clean, mCC, cc_modulated)


#%%
def main():
    """ 
    Main to demo the ecg_segmented_beat_modulation_noise_removal function.
    """
    
    # Set (default) plotting parameters.
    mpl.rcParams['font.sans-serif'] = "Calibri"
    mpl.rcParams['font.family'] = "sans-serif"
    plt.rc('axes.spines', top=False, right=False)

    blue = np.array([30,150,255])/255
    gray = [0.66,0.66,0.66]
    gray_dark = [0.33, 0.33, 0.33]
    

    #%% Load data.
    # Data is contained in "ecg.mat".
    # The original data was obtained from the MIT-BIH Polysomnographic
    # Database (http://www.physionet.org/cgi-bin/atm/ATM).
    # - Database: Motion Artifact Contaminated ECG Database (macegdb)
    # - Record: test05_45s
    # - Signals: all
    # - Length: to end
    # - Time format: samples
    # - Data format: standard
    ecg = io.loadmat(PATH_DATA/'ecg.mat')['val']

    # Each row corresponds to different signals.
    # For now, we will work with ECG3.
    ecg_raw = ecg[2, :]
    
    # Data parameters (from the .info file).
    fs = 500    # Sampling rate [Hz].
    Ts = 1/fs   # Sampling period [s].
    gain = 100  # Gain [u.a.]    
    
    f_nyquist = fs / 2
    
    # Time vector.
    t = np.arange(0, len(ecg_raw)*Ts, Ts)


    #%% Define execution parameters.
    delta_t = 40e-3 # Suggested value is 40e-3 (p. 50) [s]
    delta_t_n = delta_t * fs # [samples]
    
    
    #%% Pre-process data.
    # Convert to mV.
    ecg_mV = ecg_raw / gain     # [mV]
    
    
    #%% Filter ECG (as suggested in the paper).
    
    # Filter parameters.
    n = 8 # Desired (effective) order.
    f_low = 0.5 # Cut-off frequency 1 (suggested value is 0.5) [Hz]
    f_high = 45 # Cut-off frequency 2 (suggested value is 45) [Hz]
    
    # When creating the filter, we will divide the order by 2.
    # This is because we will be doing forward and backward filtering later on.
    # Thus, we wish to keep the original filter specification.
    b, a = signal.butter(round(n/2), [f_low/f_nyquist, f_high/f_nyquist], btype='band')

    # Actual filtering.
    ecg_pre = signal.filtfilt(b, a, ecg_mV) # [mV]
    
    # Plot
    fig, ax = plt.subplots(1, 1, figsize=[10, 5])
    plt.plot(t, ecg_mV, color=gray)
    plt.plot(t, ecg_pre, color=blue)
    ax.legend(["Raw", "Pre-processed"], frameon=False)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Amplitude [mV]")
    
    
    #%% Find R peaks.
    # Here, we will use a simple thresholding operation.
    # However, there are better (i.e., more robust) ways to do so.
    peak_threshold = 0.7;
    r_peaks, r_peaks_prop = signal.find_peaks(ecg_pre, height=peak_threshold, distance=delta_t_n)
                
    
    # Plot.
    fig, ax = plt.subplots(1, 1, figsize=[10, 5])
    plt.plot(t, ecg_pre, color=blue)
    plt.plot(r_peaks*Ts, r_peaks_prop['peak_heights'], 'or')
    ax.legend(["ECG", "Peaks"], frameon=False)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Amplitude [mV]")
    
                
    #%% Segmented beat modulation processing.
    (ecg_clean, mCC, cc_modulated) = ecg_segmented_beat_modulation_noise_removal(ecg_pre, fs, r_peaks, delta_t);

    # Plots.
    
    # mCC Computation.
    # Equivalent of Fig. 3 of the paper.
    t_mCC = np.arange(0, len(mCC)*Ts, Ts)
    
    fig, ax = plt.subplots(1, 1, figsize=[5, 2.5])
    for ii, cc_modulated_ in enumerate(cc_modulated):
        if ii == 1:
            line_legend = "Noisy cycle"
        else:
            line_legend = None
        plt.plot(t_mCC, cc_modulated_, color=gray, label=line_legend)
    plt.plot(t_mCC, mCC, color=blue, label="Median cycle")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Amplitude [mV]")
    ax.set_ylim([-0.75, 1.25])
    ax.legend(frameon=False)
    fig.savefig(PATH_IMAGES/'ecg_template.png', dpi=1000, bbox_inches='tight')
    
    
    # Signal and motion artifact.
    # Equivalent of Fig. 5 of the paper.
    # Trim signal.
    ecg_pre2 = ecg_pre[int(r_peaks[0]-delta_t_n):int(r_peaks[0]-delta_t_n+len(ecg_clean))]
    artifact = ecg_pre2 - ecg_clean
    t_artifact = np.arange(0, len(artifact)*Ts, Ts)
    
    fig, ax = plt.subplots(1, 1, figsize=[10, 5])
    plt.plot(t_artifact, artifact, color=gray)
    plt.plot(t_artifact, ecg_clean, color=blue)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Amplitude [mV]")
    ax.set_ylim([-0.75, 1.25])
    ax.legend(["Motion artifact", "Clean ECG"], frameon=False)
    

    # Signal comparison.
    # To make a better comparison, we align the clean ECG to the
    # original ECG with NaN.
    # We could, of course, with further processing align the clean ECG with
    # actual filtered ECG data.
    ecg_clean2 = np.full_like(ecg_mV, np.nan)
    ecg_clean2[int(r_peaks[0]-delta_t_n):int(r_peaks[0]-delta_t_n+len(ecg_clean))] = ecg_clean

    fig, ax = plt.subplots(1, 1, figsize=[10, 5])
    plt.plot(t, ecg_mV, color=gray_dark)
    plt.plot(t, ecg_pre, color=gray)
    plt.plot(t, ecg_clean2, color=blue)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Amplitude [mV]")
    ax.set_ylim([-0.75, 1.25])
    ax.legend(["Raw", "Pre-processed", "Clean"], frameon=False)
    fig.savefig(PATH_IMAGES/'ecg_comparison.png', dpi=1000, bbox_inches='tight')
    
    
#%%
if __name__ == '__main__':
    main()