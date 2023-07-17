% -----------------------------------------
% Preprocessing script for BASIS fNIRS ----
% Borja Blanco - bb579@cam.ac.uk
% -----------------------------------------

% Initialize environment
clear; close all; clc

% Add required toolboxes/functions to path
addpath(genpath('/Users/borjablanco/GitHub/BASIS'))
addpath(genpath('/Users/borjablanco/Documents/MATLAB/homer2'))

% Paths to raw data
path_cont = '/Users/borjablanco/Library/CloudStorage/OneDrive-UniversityofCambridge/BB_documents/BASIS/STAARS_backup/STAARS/BB_analysis/analysis_072021/control';
path_highRisk = '/Users/borjablanco/Library/CloudStorage/OneDrive-UniversityofCambridge/BB_documents/BASIS/STAARS_backup/STAARS/BB_analysis/analysis_072021/highRisk';

% Paths for storing preprocessed data and figures
path_prep_control = '/Users/borjablanco/Library/CloudStorage/OneDrive-UniversityofCambridge/BB_documents/BASIS/STAARS_backup/STAARS/BB_analysis/analysis_072021/control/preprocessed';
path_prep_highRisk = '/Users/borjablanco/Library/CloudStorage/OneDrive-UniversityofCambridge/BB_documents/BASIS/STAARS_backup/STAARS/BB_analysis/analysis_072021/highRisk/preprocessed';
path_figures = '/Users/borjablanco/Library/CloudStorage/OneDrive-UniversityofCambridge/BB_documents/BASIS/STAARS_backup/STAARS/BB_analysis/analysis_072021/figures';

% List participants 
cd(path_cont)
sub_cont = dir('*.mat');

cd(path_highRisk)
sub_highRisk = dir('*.mat');

% Parameters for motion detection
% Homer2 - hmrMotionArtifact
tMotion = 1;
tMask = 1;
STDEVthresh = 15;
AMPthresh = 0.5;

% Parameters for motion correction
% Homer2 - Spline and Wavelet denoising
p = 0.99;
iqr = 0.8;

% Time range for block averaging
tRange = [-4 20];

% Time range for stimulus rejection
tRange_tr = [-4 16];

% Bandpass filtering
hpf = 0;
lpf = 0.6;

% Parameters for channel rejection (Pollonini et al., 2016)
sci_th = 0.7;
power_th = 0.05;

% Reorder channels according to 
% 1_Programs for NIRS analysis_newsoftware.doc
nch = 26;
chlab = [6 8 7 9 16 10 17 19 18 20 25 1 11 2 3 ...
    12 4 13 22 5 14 23 15 24 21 26]; 
chlab = [chlab, chlab + nch];
% ----------------------------------------

% CHANGE THIS PART DEPENDING ON THE GROUP BEING ANALYZED
data_path = path_highRisk;
prep_path = path_prep_highRisk;
sub = sub_highRisk;
% ----------------------------------------

cd(data_path)

for nsub = 1:length(sub)
    
    % Load participant data
    cd(data_path)
    data = load(sub(nsub).name);
    data = data.data;
    clc
    disp(['Preprocessing participant: ' data.name '  '])
    
    % Add info to data structure
    data.name = sub(nsub).name;
    data.DPF = [5.29 4.25]; % From Scholkmann 2013
    data.sf = 10; % Sampling frequency
    data.time = (0:size(data.d,1)-1)/data.sf;
    data.nchannels = size(data.d,2);
    data.SD.MeasListAct = ones(size(data.d,2),1);
    
    % Reorder channel according to BASIS layout (for plotting)
    data.d = data.d (:, chlab);
    data.SD.MeasList = data.SD.MeasList(chlab,:);
         
    % Display number of included trials per condition
    % At this point trials have been excluded due to looking times
    %disp(['Included trials condition S: ' num2str(data.lookTrials(1).trials) '  ']);
    %disp(['Included trials condition N: ' num2str(data.lookTrials(2).trials) '  ']);
    %disp(['Included trials condition V: ' num2str(data.lookTrials(3).trials) '  ']);
    
    % Exclude channels (Pollonini et al., 2016)
    [data.bad_links, data.bad_windows] = sciLUCA_Ch_BB(data, path_figures, sci_th, power_th);
    data.SD.MeasListAct(data.bad_links) = 0;    
    disp([num2str(data.bad_links') ' excluded channels']);

    % Convert raw intensity to optical density
    data.OD = hmrIntensity2OD(data.d);
    
    % Motion detection (by channel)
    [~, tIncCh] = hmrMotionArtifactByChannel(data.OD, data.sf, data.SD, [], tMotion, tMask, STDEVthresh, AMPthresh);
    
    % Motion correction (spline)
    data.spline = hmrMotionCorrectSpline (data.OD, data.t, data.SD, tIncCh,  p);
    
    % Motion correction (wavelets)
    data.wav = hmrMotionCorrectWavelet(data.spline, data.SD, iqr);
    
    % Motion detection (for trial rejection)
    [tInc, ~] = hmrMotionArtifactByChannel(data.wav, data.sf, data.SD, [], tMotion, tMask, STDEVthresh, AMPthresh);
    
    % Exclude trials based on looking times and update data.s
    % data.s colums  = C S N V (4)
    % data.lookTrials =  S N V (3) that's why indexing differs
    idx_S = find(data.s(:,2));
    data.s(idx_S, 2) = -1;
    data.s(idx_S(data.lookTrials(1).trials), 2) = 1;
    
    idx_N = find(data.s(:,3));
    data.s(idx_N, 3) = -1;
    data.s(idx_N(data.lookTrials(2).trials), 3) = 1;
    
    idx_V = find(data.s(:,4));
    data.s(idx_V, 4) = -1;
    data.s(idx_V(data.lookTrials(3).trials), 4) = 1;
    
    % Trial rejection
    [data.s, ~] = enStimRejection (data.t, data.s, tInc, [], tRange_tr);
     
    % Optical density to concentration
    data.conc = hmrOD2Conc (data.wav, data.SD, data.DPF);
    
    % Low pass filtering (hpf=0)
    data.filt = hmrBandpassFilt (data.conc, data.sf, hpf, lpf);

    % Block averaging (detrending + normalization (mean))
    [hrf_results, ystd, tHRF, nTrials, ysum2, yTrials] = bb_hmrBlockAvg(data.filt, data.s, data.t, tRange);
    data.ystd = ystd;
    data.tHRF = tHRF;
    data.nTrials = nTrials;
    data.ysum2 = ysum2;
    data.yTrials = yTrials;
    
    % Store results (HRF, Hb, channels, conditions)
    data.hbo_S = squeeze(hrf_results(:,1,:,2)); % hbo visual social
    data.hbr_S = squeeze(hrf_results(:,2,:,2)); % hbr visual social
    data.hbo_N = squeeze(hrf_results(:,1,:,3)); % hbo non-vocal
    data.hbr_N = squeeze(hrf_results(:,2,:,3)); % hbr non-vocal
    data.hbo_V = squeeze(hrf_results(:,1,:,4)); % hbo vocal
    data.hbr_V = squeeze(hrf_results(:,2,:,4)); % hbr vocal

    % Save preprocessed data
    cd(prep_path)
    save([data.name(1:end-5), '_prepro.mat'], 'data')
       
end



