% clear; clc; 
% close all;
format compact;
% set(0,'DefaultFigureWindowStyle','normal')

addpath(genpath('main'));
addpath(genpath('chronux'));
addpath(genpath('Kilosort'));
addpath(genpath('npy-matlab'));
addpath(genpath('spikes-master'));
addpath(genpath('PreProcessing'));
addpath(genpath('Plotting'));
addpath(genpath('Analysis'));
addpath(genpath('Dependancies'));
rmpath(genpath('Dependancies/MVGC1'));
%%  PreProcessing
load GridsLowDenNeedle_chanmap.mat;  % load the channel map for the IntanConcatenate function
parameters.rows = 8;  % Number of rows of electrodes on the Grid
parameters.cols = 4;  % Number of colums of electrodes on the Grid
parameters.Fs = 1000;
parameters.ts = 1/parameters.Fs;
parameters.windowBeforePull = 1.5; % in seconds
parameters.windowAfterPull = 1.5; % in seconds
parameters.windowBeforeCue = 1.5; % in seconds 
parameters.windowAfterCue = 1.5; % in seconds 
parameters.windowBeforeMI = 1.5; % in seconds 
parameters.windowAfterMI = 1.5; % in seconds 
parameters.experiment = 'cue'; % self - internally generated, cue - cue initiated 
parameters.opto = 0; % 1 - opto ON , 0 - opto OFF
parameters.xspacing = 0.06; % Grid spacing in mm between columns 
parameters.yspacing = 0.06; % Grid spacing in mm between rows
parameters.shank = 0; % 1 - if UCLA 64Ch single shank data is present

if parameters.shank == 1
    load UCLASingle64Ch_chanmap.mat; % load the channel map for the the shank data
    UCLAProbeMap = UCLAProbeMap + size(electrode_map,1);
    finalElectrodeMap = [electrode_map;UCLAProbeMap];
    parameters.nShank = size(UCLAProbeMap,1);
else
    finalElectrodeMap = electrode_map;
end

IntanConcatenate
fpath = Intan.path; % where on disk do you want the analysis? ideally and SSD...

% Generating time series from Intan data
Ts = 1/Intan.offsetSample;
Intan.Tmax = Ts * size(Intan.allIntan,2);
Intan.t = 0:Ts:Intan.Tmax-Ts;

%% Removing bad channels from impedance values
[Z,Intan.goodChMap,Intan.badChMap] = readImp(electrode_map,10e6);
figure('Name','Impedance Test at 1kHz');boxchart(Z); xlabel('n = ' + string(size(Z,1)));ylabel('Impedance (in \Omega)');set(gca,'xticklabel',{[]})
% Intan.badChMap =[]; [1;2;7];
%Intan = removeBadCh(Intan,Intan.badCh);
%% LFP
set(0,'DefaultFigureWindowStyle','normal')
LFP = fastpreprocess_filtering(Intan.allIntan,5000);
% LFP = bestLFP(LFP);
% LFP = bandFilter(LFP,'depth'); % Extract LFPs based on 'depth' or 'single'
% LFPplot(LFP);
LFP = createDataCube(LFP,parameters.rows,parameters.cols,Intan.badChMap); % Creating datacube

%% Loading Lever Data 
plotOption = 0;
[Behaviour] = readLever(parameters,LFP.times,plotOption);

%% Reading behaviour data from Intan traces 
plotOption = 1;
IntanBehaviour = readLeverIntan(parameters,LFP.times,Intan.analog_adc_data,Intan.dig_in_data,Behaviour,plotOption);

%% Generalized Phase 
LFP.xf = bandpass_filter(LFP.LFPdatacube,5,40,4,1000);
[LFP.xgp, LFP.wt] = generalized_phase(LFP.xf,1000,0);
% LFP.xfbeta = bandpass_filter(LFP.LFPdatacube,10,30,4,1000);
% [LFP.xgpbeta, LFP.wtbeta] = generalized_phase(LFP.xfbeta,1000,0);
% LFP.xftheta = bandpass_filter(LFP.LFPdatacube,4,10,4,1000);
% [LFP.xgptheta, LFP.wttheta]  = generalized_phase(LFP.xftheta,1000,0);
% LFP.xfgamma = bandpass_filter(LFP.LFPdatacube,30,40,4,1000);
% [LFP.xgpgamma, LFP.wtgamma]  = generalized_phase(LFP.xfgamma,1000,0);
% [parameters.X,parameters.Y] = meshgrid( 1:parameters.cols, 1:parameters.rows );
% LFP.xfwide = bandpass_filter(LFP.LFPdatacube,5,90,4,1000);
% LFP.xfbetanarrow = bandpass_filter(LFP.LFPdatacube,6,9,4,1000);
% [LFP.xgpbetanarrow, LFP.wtbetanarrow] = generalized_phase(LFP.xfbetanarrow,1000,0);
% % GP for spatial mean LFP 
% [LFP.xgpbetamean, ~] = generalized_phase(mean(LFP.xfbetanarrow,[1,2]),1000,0);

% Add trial segmented data to IntanBehaviour Variable
IntanBehaviour = addLFPToBehaviour(IntanBehaviour,LFP,parameters);
%% Saving paramters, path, IntanBehaviour to bin file 
savepath = uigetdir(path);
sessionName = [savepath,'/','M1WavesCooling.mat'];
% save(sessionName,"IntanBehaviour","fpath","parameters","-v7.3");
save(sessionName,"IntanBehaviour","fpath","parameters","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",

%% Combining  multiple Intanbehaviour structs from multiple sessions
combIntanBehaviour = horzcat(IntanBehaviour1, IntanBehaviour2);
IntanBehaviour.cueHitTrace = horzcat(combIntanBehaviour(1:end).cueHitTrace);
IntanBehaviour.cueMissTrace = horzcat(combIntanBehaviour(1:end).cueMissTrace);
IntanBehaviour.hitTrace = horzcat(combIntanBehaviour(1:end).hitTrace);
IntanBehaviour.missTrace = horzcat(combIntanBehaviour(1:end).missTrace);
IntanBehaviour.reactionTime = horzcat(combIntanBehaviour(1:end).reactionTime);
IntanBehaviour.MIHitTrace = horzcat(combIntanBehaviour(1:end).MIHitTrace);
IntanBehaviour.MIFATrace = horzcat(combIntanBehaviour(1:end).MIFATrace);

WavesComb = horzcat(Waves1,Waves2);
Waves.wavesHit = horzcat(WavesComb(1:end).wavesHit);
Waves.wavesMiss = horzcat(WavesComb(1:end).wavesMiss);
Waves.wavesHitReward = horzcat(WavesComb(1:end).wavesHitReward);
Waves.wavesFA = horzcat(WavesComb(1:end).wavesFA);
Waves.wavesMIHit = horzcat(WavesComb(1:end).wavesMIHit);
Waves.wavesMIFA = horzcat(WavesComb(1:end).wavesMIFA);

clear combIntanBehaviour IntanBehaviour1 IntanBehaviour2 WavesComb Waves1 Waves2;
%% Power Spectrum during task across channels 
[PSD.ChHit , PSD.f] = getAvgPSD(IntanBehaviour.cueHitTrace,parameters);
[PSD.ChMiss , PSD.f] = getAvgPSD(IntanBehaviour.cueMissTrace,parameters);
PSD.trialPSDHit = squeeze(10*log10(mean(PSD.ChHit,2,"omitnan")));
PSD.trialPSDMiss = squeeze(10*log10(mean(PSD.ChMiss,2,"omitnan")));

[PSD.ChHitReward , PSD.f] = getAvgPSD(IntanBehaviour.hitTrace,parameters);
[PSD.ChFA , PSD.f] = getAvgPSD(IntanBehaviour.missTrace,parameters);
PSD.trialPSDHitReward = squeeze(10*log10(mean(PSD.ChHitReward,2,"omitnan")));
PSD.trialPSDFA = squeeze(10*log10(mean(PSD.ChFA,2,"omitnan")));

[PSD.ChHitMI , PSD.f] = getAvgPSD(IntanBehaviour.MIHitTrace,parameters);
[PSD.ChFAMI , PSD.f] = getAvgPSD(IntanBehaviour.MIFATrace,parameters);
PSD.trialPSDHitMI = squeeze(10*log10(mean(PSD.ChHitMI,2,"omitnan")));
PSD.trialPSDFAMI = squeeze(10*log10(mean(PSD.ChFAMI,2,"omitnan")));

remove_artifact = 1;
if remove_artifact == 1
    % Removing trials with artifacts
    PSD.rejectThres = 39; % in db
    PSD.rejectFreq = 3; %  not Hz but index in PSD.f
    PSD.artifactTrialIndex = find(PSD.trialPSDHit(:,PSD.rejectFreq)>PSD.rejectThres);
    disp(['Number of Hit trials rejected ', num2str(size(PSD.artifactTrialIndex,1))]);
    PSD.trialPSDHit(PSD.artifactTrialIndex,:) = [];
    IntanBehaviour.cueHitTrace(PSD.artifactTrialIndex)=[];
    PSD.ChHit(PSD.artifactTrialIndex,:,:) = [];
%     IntanBehaviour.reactionTime(PSD.artifactTrialIndex,:,:) = [];
    
    PSD.artifactTrialIndex = find(PSD.trialPSDMiss(:,PSD.rejectFreq)>PSD.rejectThres);
    disp(['Number of Miss trials rejected ', num2str(size(PSD.artifactTrialIndex,1))]);
    PSD.trialPSDMiss(PSD.artifactTrialIndex,:) = [];
    IntanBehaviour.cueMissTrace(PSD.artifactTrialIndex)=[];
    PSD.ChMiss(PSD.artifactTrialIndex,:,:) = [];
    
    PSD.artifactTrialIndex = find(PSD.trialPSDHitReward(:,PSD.rejectFreq)>PSD.rejectThres);
    disp(['Number of Hit trials rejected ', num2str(size(PSD.artifactTrialIndex,1))]);
    PSD.trialPSDHitReward(PSD.artifactTrialIndex,:) = [];
    IntanBehaviour.hitTrace(PSD.artifactTrialIndex)=[];
    PSD.ChHitReward(PSD.artifactTrialIndex,:,:) = [];
    
    PSD.artifactTrialIndex = find(PSD.trialPSDFA(:,PSD.rejectFreq)>PSD.rejectThres);
    disp(['Number of FA trials rejected ', num2str(size(PSD.artifactTrialIndex,1))]);
    PSD.trialPSDFA(PSD.artifactTrialIndex,:) = [];
    IntanBehaviour.missTrace(PSD.artifactTrialIndex)=[];
    PSD.ChFA(PSD.artifactTrialIndex,:,:) = [];

    PSD.artifactTrialIndex = find(PSD.trialPSDHitMI(:,PSD.rejectFreq)>PSD.rejectThres);
    disp(['Number of Hit MI trials rejected ', num2str(size(PSD.artifactTrialIndex,1))]);
    PSD.trialPSDHitMI(PSD.artifactTrialIndex,:) = [];
    IntanBehaviour.MIHitTrace(PSD.artifactTrialIndex)=[];
    PSD.ChHitMI(PSD.artifactTrialIndex,:,:) = [];

    PSD.artifactTrialIndex = find(PSD.trialPSDFAMI(:,PSD.rejectFreq)>PSD.rejectThres);
    disp(['Number of FA MI trials rejected ', num2str(size(PSD.artifactTrialIndex,1))]);
    PSD.trialPSDFAMI(PSD.artifactTrialIndex,:) = [];
    IntanBehaviour.MIFATrace(PSD.artifactTrialIndex)=[];
    PSD.ChFAMI(PSD.artifactTrialIndex,:,:) = [];
end

PSD.avgPSDHit = squeeze(10*log10(mean(PSD.ChHit,[1 2],"omitnan")));
PSD.avgPSDMiss = squeeze(10*log10(mean(PSD.ChMiss,[1 2],"omitnan")));
PSD.avgPSDHitReward = squeeze(10*log10(mean(PSD.ChHitReward,[1 2],"omitnan")));
PSD.avgPSDFA = squeeze(10*log10(mean(PSD.ChFA,[1 2],"omitnan")));
PSD.avgPSDHitMI = squeeze(10*log10(mean(PSD.ChHitMI,[1 2],"omitnan")));
PSD.avgPSDFAMI = squeeze(10*log10(mean(PSD.ChFAMI,[1 2],"omitnan")));

figure();
subplot(1,2,1);
plot(PSD.f(1:81),PSD.trialPSDHit(:,1:81),'Color', [0 0 1 0.1]);
hold on;
plot(PSD.f(1:81),PSD.avgPSDHit(1:81),'Color', [0 0 1 1],'LineWidth',1.5);
ylim([0 50]);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');
title('Average Power Spectral Density for Hit Trials');
box off;

subplot(1,2,2);
plot(PSD.f(1:81),PSD.trialPSDMiss(:,1:81),'Color', [1 0 0 0.1]);
hold on;
plot(PSD.f(1:81),PSD.avgPSDMiss(1:81),'Color', [1 0 0 1],'LineWidth',1.5);
ylim([0 50]);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');
title('Average Power Spectral Density for Miss Trials');
box off;

figure();
subplot(1,2,1);
plot(PSD.f(1:81),PSD.trialPSDHitReward(:,1:81),'Color', [0 0 1 0.1]);
hold on;
plot(PSD.f(1:81),PSD.avgPSDHitReward(1:81),'Color', [0 0 1 1],'LineWidth',1.5);
ylim([0 50]);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');
title('Average Power Spectral Density for Hit Trials');
box off;

subplot(1,2,2);
plot(PSD.f(1:81),PSD.trialPSDFA(:,1:81),'Color', [1 0 0 0.1]);
hold on;
plot(PSD.f(1:81),PSD.avgPSDFA(1:81),'Color', [1 0 0 1],'LineWidth',1.5);
ylim([0 50]);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');
title('Average Power Spectral Density for FA Trials');
box off;

figure();
subplot(1,2,1);
plot(PSD.f(1:81),PSD.trialPSDHitMI(:,1:81),'Color', [0 0 1 0.1]);
hold on;
plot(PSD.f(1:81),PSD.avgPSDHitMI(1:81),'Color', [0 0 1 1],'LineWidth',1.5);
ylim([0 50]);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');
title('Average Power Spectral Density for Hit MI Trials');
box off;

subplot(1,2,2);
plot(PSD.f(1:81),PSD.trialPSDFAMI(:,1:81),'Color', [1 0 0 0.1]);
hold on;
plot(PSD.f(1:81),PSD.avgPSDFAMI(1:81),'Color', [1 0 0 1],'LineWidth',1.5);
ylim([0 50]);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');
title('Average Power Spectral Density for FA MI Trials');
box off;

%% Wavelet spectrogram
[Spectrogram.hitAvgSpectrogram, Spectrogram.hitSpectrogramCWT,Spectrogram.AvgHitTrace ,Spectrogram.fwt] = getAvgSpectogram(IntanBehaviour.cueHitTrace,parameters,[5 80]);
[Spectrogram.missAvgSpectrogram, Spectrogram.missSpectrogramCWT,Spectrogram.AvgMissTrace,Spectrogram.fwt] = getAvgSpectogram(IntanBehaviour.cueMissTrace,parameters,[5 80]);
[Spectrogram.hitRewardAvgSpectrogram, Spectrogram.hitRewardSpectrogramCWT,Spectrogram.AvgHitRewardTrace ,Spectrogram.fwt] = getAvgSpectogram(IntanBehaviour.hitTrace,parameters,[5 80]);
[Spectrogram.FAAvgSpectrogram, Spectrogram.FASpectrogramCWT,Spectrogram.AvgFATrace,Spectrogram.fwt] = getAvgSpectogram(IntanBehaviour.missTrace,parameters,[5 80]);
[Spectrogram.hitMIAvgSpectrogram, Spectrogram.hitMISpectrogramCWT,Spectrogram.AvgHitMITrace ,Spectrogram.fwt] = getAvgSpectogram(IntanBehaviour.MIHitTrace,parameters,[5 80]);
[Spectrogram.FAMIAvgSpectrogram, Spectrogram.FAMISpectrogramCWT,Spectrogram.AvgFAMITrace,Spectrogram.fwt] = getAvgSpectogram(IntanBehaviour.MIFATrace,parameters,[5 80]);
% Global average spectogram
figure('Name','Trial Averaged Wavelet Spectrogram for Hits & Misses');
subplot(1,2,1);
plotSpectrogram(10*log10((squeeze(Spectrogram.hitAvgSpectrogram))),IntanBehaviour.cueHitTrace(1).time,Spectrogram.fwt,'surf','Wavelet Based Spectrogram for Hits','Time (s)','Frequency (Hz)')
caxis([-5 12]);
hold on; yyaxis right; box off;
plot(IntanBehaviour.cueHitTrace(1).time,Spectrogram.AvgHitTrace,'-w','LineWidth',2.5);
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
ylabel('Lever deflection (mV)'); ylim([0 0.1]);
subplot(1,2,2);
plotSpectrogram(10*log10((squeeze(Spectrogram.missAvgSpectrogram))),IntanBehaviour.cueMissTrace(1).time,Spectrogram.fwt,'surf','Wavelet Based Spectrogram for Misses','Time (s)','Frequency (Hz)')
caxis([-5 12]);
hold on; yyaxis right; box off;
plot(IntanBehaviour.cueMissTrace(1).time,Spectrogram.AvgMissTrace,'-w','LineWidth',2.5);
xline(0,'--r','Cue','LabelVerticalAlignment','top');
ylabel('Lever deflection (mV)'); ylim([0 0.1]); box off;

% Global average spectogram
figure('Name','Trial Averaged Wavelet Spectrogram for Hits & FA');
subplot(1,2,1);
plotSpectrogram(10*log10((squeeze(Spectrogram.hitRewardAvgSpectrogram))),IntanBehaviour.hitTrace(1).time,Spectrogram.fwt,'surf','Wavelet Based Spectrogram for Hits','Time (s)','Frequency (Hz)')
caxis([-5 12]);
hold on; yyaxis right; box off;
plot(IntanBehaviour.hitTrace(1).time,Spectrogram.AvgHitRewardTrace,'-w','LineWidth',2.5);
xline(1,'--r','Reward','LabelVerticalAlignment','top');
% xline(-1*mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Cue Time','LabelVerticalAlignment','top');
ylabel('Lever deflection (mV)'); ylim([0 0.1]);
subplot(1,2,2);
plotSpectrogram(10*log10((squeeze(Spectrogram.FAAvgSpectrogram))),IntanBehaviour.missTrace(1).time,Spectrogram.fwt,'surf','Wavelet Based Spectrogram for Misses','Time (s)','Frequency (Hz)')
caxis([-5 12]);
hold on; yyaxis right; box off;
plot(IntanBehaviour.missTrace(1).time,Spectrogram.AvgFATrace,'-w','LineWidth',2.5);
xline(0.5,'--r','Threshold Time','LabelVerticalAlignment','top');
ylabel('Lever deflection (mV)'); ylim([0 0.1]); box off;

%% Wave detection in velocity triggered windows
% nShuffle = 1000;
% threshold = 95;
% trialno = 67;
nShuffle = 100;
threshold = 99; % zscore of 3
parameters.rhoThres = getRhoThreshold2(IntanBehaviour.cueHitTrace,IntanBehaviour.cueMissTrace,parameters,nShuffle,threshold);
% Wave detection for wide band
disp('Wave Detection for wide band ...')
xgp = arrayfun(@(s) s.xgp, IntanBehaviour.hitTrace, 'UniformOutput', false);
if isfield(IntanBehaviour,'cueHitTrace')
    xf = arrayfun(@(s) s.xf, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgp, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wt, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
    Waves.wavesHit = detectWaves(xf,xgp,wt,IntanBehaviour.cueHitTrace,parameters,parameters.rhoThres);
end
% Waves.wavesHit = detectPlanarWaves(xf,xgp,wt,IntanBehaviour.cueHitTrace,parameters,0.5);
if isfield(IntanBehaviour,'cueMissTrace')
    xf = arrayfun(@(s) s.xf, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgp, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wt, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    Waves.wavesMiss = detectWaves(xf,xgp,wt,IntanBehaviour.cueMissTrace,parameters,parameters.rhoThres);
%     Waves.wavesMiss = detectPlanarWaves(xf,xgp,wt,IntanBehaviour.cueMissTrace,parameters,0.5);
end
if isfield(IntanBehaviour,'MIHitTrace')
    xf = arrayfun(@(s) s.xf, IntanBehaviour.MIHitTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgp, IntanBehaviour.MIHitTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wt, IntanBehaviour.MIHitTrace, 'UniformOutput', false);
    Waves.wavesMIHit = detectWaves(xf,xgp,wt,IntanBehaviour.MIHitTrace,parameters,parameters.rhoThres);
end
if isfield(IntanBehaviour,'MIFATrace')
    xf = arrayfun(@(s) s.xf, IntanBehaviour.MIFATrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgp, IntanBehaviour.MIFATrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wt, IntanBehaviour.MIFATrace, 'UniformOutput', false);
    Waves.wavesMIFA = detectWaves(xf,xgp,wt,IntanBehaviour.MIFATrace,parameters,parameters.rhoThres);
end
if isfield(IntanBehaviour,'missTrace')
    xf = arrayfun(@(s) s.xf, IntanBehaviour.missTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgp, IntanBehaviour.missTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wt, IntanBehaviour.missTrace, 'UniformOutput', false);
    Waves.wavesFA = detectWaves(xf,xgp,wt,IntanBehaviour.missTrace,parameters,parameters.rhoThres);
end
if isfield(IntanBehaviour,'hitTrace')
    xf = arrayfun(@(s) s.xf, IntanBehaviour.hitTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgp, IntanBehaviour.hitTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wt, IntanBehaviour.hitTrace, 'UniformOutput', false);
    Waves.wavesHitReward = detectWaves(xf,xgp,wt,IntanBehaviour.hitTrace,parameters,parameters.rhoThres);
end

%% Wave detection for theta band
disp('Wave Detection for theta band ...')
threshold = 99;
xf = arrayfun(@(s) s.xftheta, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
xgp = arrayfun(@(s) s.xgptheta, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
wt = arrayfun(@(s) s.wttheta, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
% parameters.thetarhoThres = getRhoThreshold(xgp,IntanBehaviour.cueHitTrace,parameters,nShuffle,trialno,threshold);
thetaWaves.wavesHit = detectWaves(xf,xgp,wt,IntanBehaviour.cueHitTrace,parameters,parameters.thetarhoThres);
if isfield(IntanBehaviour,'cueMissTrace')
    xf = arrayfun(@(s) s.xftheta, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgptheta, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wttheta, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    thetaWaves.wavesMiss = detectWaves(xf,xgp,wt,IntanBehaviour.cueMissTrace,parameters,parameters.thetarhoThres);
end

% Wave detection for beta band
disp('Wave Detection for beta band ...')
threshold = 99;
xf = arrayfun(@(s) s.xfbeta, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
xgp = arrayfun(@(s) s.xgpbeta, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
wt = arrayfun(@(s) s.wtbeta, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
% parameters.betarhoThres = getRhoThreshold(xgp,IntanBehaviour.cueHitTrace,parameters,nShuffle,trialno,threshold);
betaWaves.wavesHit= detectWaves(xf,xgp,wt,IntanBehaviour.cueHitTrace,parameters,parameters.betarhoThres);
if isfield(IntanBehaviour,'cueMissTrace')
    xf = arrayfun(@(s) s.xfbeta, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgpbeta, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wtbeta, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    betaWaves.wavesMiss = detectWaves(xf,xgp,wt,IntanBehaviour.cueMissTrace,parameters,parameters.betarhoThres);
end

% Wave detection for gamma band
disp('Wave Detection for gamma band ...')
threshold = 99;
xf = arrayfun(@(s) s.xfgamma, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
xgp = arrayfun(@(s) s.xgpgamma, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
wt = arrayfun(@(s) s.wtgamma, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
% parameters.gammarhoThres = getRhoThreshold(xgp,IntanBehaviour.cueHitTrace,parameters,nShuffle,trialno,threshold);
gammaWaves.wavesHit = detectWaves(xf,xgp,wt,IntanBehaviour.cueHitTrace,parameters,parameters.gammarhoThres);
if isfield(IntanBehaviour,'cueMissTrace')
    xf = arrayfun(@(s) s.xfgamma, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgpgamma, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wtgamma, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    gammaWaves.wavesMiss = detectWaves(xf,xgp,wt,IntanBehaviour.cueMissTrace,parameters,parameters.gammarhoThres);
end

%% Wave detecion for entire time 
parameters.rhoThres= rhoThres;
allwaves.LFPIndex = (1:1:size(LFP.LFP,2))';
xf{1,1} = LFP.xf;
xgp{1,1} = LFP.xgp;
wt{1,1} = LFP.wt;
Wavesall = detectWaves(xf,xgp,wt,allwaves,parameters,0.6);
WaveStatsSingle(Wavesall,parameters,1,size(LFP.LFP,2)/20);

%% PLotting to check visually
% trialPlot = 48;
% plot_wave_examples( LFP.xf(:,:,IntanBehaviour.cueHitTrace(trialPlot).LFPIndex(1):IntanBehaviour.cueHitTrace(trialPlot).LFPIndex(end)), ...
%     options, trialPlot, Waves.wavesHit,rhoThres);

% trialPlot = 12;
% plot_wave_examples( LFP.xf(:,:,IntanBehaviour.missTrace(trialPlot).LFPIndex(1):IntanBehaviour.missTrace(trialPlot).LFPIndex(end)), ...
%     options, trialPlot, Waves.wavesMiss,rhoThres);

%% Waves accross trials 
plotOption = 1;
[WaveStats(1)] = getWaveStats(Waves.wavesHit,Waves.wavesMiss,parameters,plotOption);
[WaveStats(2)] = getWaveStats(thetaWaves.wavesHit,thetaWaves.wavesMiss,parameters,plotOption);
[WaveStats(3)] = getWaveStats(betaWaves.wavesHit,betaWaves.wavesMiss,parameters,plotOption);
[WaveStats(4)] = getWaveStats(gammaWaves.wavesHit,gammaWaves.wavesMiss,parameters,plotOption);

plotOption = 1;
[WaveStatsReward(1)] = getWaveStats(Waves.wavesHitReward,Waves.wavesFA,parameters,plotOption);

plotOption = 1;
[WaveStatsFA(1)] = getWaveStats(Waves.wavesMIHit,Waves.wavesMIFA,parameters,plotOption);

[WaveStats2(1)] = getInitRewardStats(Waves,parameters,plotOption);
[WaveStats2(2)] = getInitRewardStats(thetaWaves,parameters,plotOption);
[WaveStats2(3)] = getInitRewardStats(betaWaves,parameters,plotOption);
[WaveStats2(4)] = getInitRewardStats(gammaWaves,parameters,plotOption);

%% Mutual Information
z_score = 0;
nIterrate = 100;
MI = getMI(IntanBehaviour,z_score,nIterrate,1,parameters);

% For amplitude
xgpHit = arrayfun(@(s) abs(s.xgp), IntanBehaviour.cueHitTrace, 'UniformOutput', false);
xgpMiss = arrayfun(@(s) abs(s.xgp), IntanBehaviour.cueMissTrace, 'UniformOutput', false);
[MI.Amp] = getMutualInformation(xgpHit,xgpMiss,parameters);

figure();
title("Mututal Information across all electrodes - Amplitude")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:32,peakSort2DArray(reshape(MI.Amp,[],size(MI.Amp,3)),'descend',2)); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)"); 
h = colorbar; h.Label.String = 'Information (bits)';
xline(0,'-w','Cue','LabelVerticalAlignment','top');


%% Average PGD 
PGD.avgPGDHit = mean(vertcat(Waves.wavesHit.PGD),1);
PGD.avgPGDMiss = mean(vertcat(Waves.wavesMiss.PGD),1);
PGD.avgPGDHitReward = mean(vertcat(Waves.wavesHitReward.PGD),1);
PGD.avgPGDFA = mean(vertcat(Waves.wavesFA.PGD),1);
PGD.avgPGDMIHit = mean(vertcat(Waves.wavesMIHit.PGD),1);
PGD.avgPGDMIFA = mean(vertcat(Waves.wavesMIFA.PGD),1);

figure(); hold on;
plot(IntanBehaviour.cueHitTrace(1).time,PGD.avgPGDHit,'-r','LineWidth',1.2); hold on;
plot(IntanBehaviour.cueHitTrace(1).time,PGD.avgPGDMiss,'-k','LineWidth',1);
% plot(IntanBehaviour.cueHitTrace(1).time,PGD.avgPGDFA,'-b','LineWidth',1);
ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD)');box off;  legend('Hits','Misses');%,'False Alarms');

figure(); hold on;
plot(IntanBehaviour.hitTrace(1).time,PGD.avgPGDHitReward,'-r','LineWidth',1.2); hold on;
plot(IntanBehaviour.missTrace(1).time,PGD.avgPGDFA,'-k','LineWidth',1);
ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','Threshold','LabelVerticalAlignment','top');
xline(-mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Cue Time','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD)');box off;  legend('Hits','FAs');%,'False Alarms');

figure(); hold on;
plot(IntanBehaviour.MIHitTrace(1).time,PGD.avgPGDMIHit,'-r','LineWidth',1.2); hold on;
plot(IntanBehaviour.MIFATrace(1).time,PGD.avgPGDMIFA,'-k','LineWidth',1);
ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','MI','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD) - Motion Ininitiation');box off;  legend('Hits','FAs');%,'False Alarms');

xgp = arrayfun(@(s) s.xgp, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
[PGD.avgPGDHit,PGD.avgPGDHitNull] = getAvgPGD(xgp,Waves.wavesHit,IntanBehaviour.cueHitTrace,1,parameters);

figure(); hold on;
plot(IntanBehaviour.cueHitTrace(1).time,PGD.avgPGDHit,'-k','LineWidth',1.2); hold on;
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(PGD.avgPGDHitNull),'-r','LineWidth',1.2);
ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD)');box off;  legend('Hits','Shuffled');

%% Cross-Trial Phase Alignment
z_score = 0;
nIterrate = 2000;
PA = getPA(IntanBehaviour,z_score,nIterrate,1,parameters,0);

%% Percent Phase Locking
xgp = arrayfun(@(s) s.xgp, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
[PPLHit] = getPPL(xgp,parameters);
% xgp = arrayfun(@(s) shuffle3DMatrix(s.xgp,3), IntanBehaviour.cueHitTrace, 'UniformOutput', false);
% [PPLHitShuffle] = getPPL(xgp,parameters);
xgp = arrayfun(@(s) s.xgp, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
[PPLMiss] = getPPL(xgp,parameters);
xgp = arrayfun(@(s) s.xgp, IntanBehaviour.hitTrace, 'UniformOutput', false);
[PPLHitReward] = getPPL(xgp,parameters);
xgp = arrayfun(@(s) s.xgp, IntanBehaviour.missTrace, 'UniformOutput', false);
[PPLFA] = getPPL(xgp,parameters);
xgp = arrayfun(@(s) s.xgp, IntanBehaviour.MIHitTrace, 'UniformOutput', false);
[PPLMIHit] = getPPL(xgp,parameters);
xgp = arrayfun(@(s) s.xgp, IntanBehaviour.MIFATrace, 'UniformOutput', false);
[PPLMIFA] = getPPL(xgp,parameters);

figure();
subplot(2,1,1);
title("Percentage Phase across all electrodes - Hits")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:32,reshape(PPLHit,[],size(PPLHit,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--w','Avg. Reaction Time','LabelVerticalAlignment','top');
subplot(2,1,2);
title("Percentage Phase across all electrodes - Misses")
imagesc(IntanBehaviour.cueMissTrace(1).time,1:32,reshape(PPLMiss,[],size(PPLMiss,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');

figure();
subplot(2,1,1);
title("Percentage Phase across all electrodes - Hits")
imagesc(IntanBehaviour.hitTrace(1).time,1:32,reshape(PPLHitReward,[],size(PPLHitReward,3))); colormap(jet);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Threshold','LabelVerticalAlignment','top');
yyaxis right; box off;
plot(IntanBehaviour.hitTrace(1).time,squeeze(mean(PPLHitReward,[1 2],'omitnan')),'-r','Linewidth',0.8);
% xline(-mean(IntanBehaviour.reactionTime,'all'),'--w','Avg. Cue Time','LabelVerticalAlignment','top');
subplot(2,1,2);
title("Percentage Phase across all electrodes - FA")
imagesc(IntanBehaviour.missTrace(1).time,1:32,reshape(PPLFA,[],size(PPLFA,3))); colormap(jet);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Threshold','LabelVerticalAlignment','top');
yyaxis right; box off;
plot(IntanBehaviour.missTrace(1).time,squeeze(mean(PPLFA,[1 2],'omitnan')),'-r','Linewidth',0.8);

figure();
subplot(2,1,1);
title("Percentage Phase across all electrodes - Hits MI")
imagesc(IntanBehaviour.MIHitTrace(1).time,1:32,reshape(PPLMIHit,[],size(PPLMIHit,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','MI','LabelVerticalAlignment','top');
subplot(2,1,2);
title("Percentage Phase across all electrodes - Misses MI")
imagesc(IntanBehaviour.MIFATrace(1).time,1:32,reshape(PPLMIFA,[],size(PPLMIFA,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','MI','LabelVerticalAlignment','top');


figure();
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PPLHit,[1 2])),'-r','LineWidth',1.2); hold on;
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PPLMiss,[1 2])),'-k','LineWidth',0.2);
ylabel("Percentage Phase Locking"); xlabel("Time (s)");
xline(0,'--k','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--k','Avg. Reaction Time','LabelVerticalAlignment','top');
title('Percentage Phase Locking for hits');box off;legend('Hits','Misses');%,'Misses');

figure();
plot(IntanBehaviour.MIHitTrace(1).time,squeeze(nanmean(PPLMIHit,[1 2])),'-r','LineWidth',1.2); hold on;
plot(IntanBehaviour.MIFATrace(1).time,squeeze(nanmean(PPLMIFA,[1 2])),'-k','LineWidth',0.2);
ylabel("Percentage Phase Locking"); xlabel("Time (s)");
xline(0,'--k','MI','LabelVerticalAlignment','top');
title('Percentage Phase Locking for hits');box off;legend('Hits','FAs');%,'Misses');

% z-scoring 
nIterrate = 200;
xgpHit = arrayfun(@(s) s.xgp, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
xgpMiss = arrayfun(@(s) s.xgp, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
xgpComb = [xgpHit xgpMiss];
nHit = size(xgpHit,2);
nMiss = size(xgpMiss,2);
nTot = nHit + nMiss;
nullDistHit = zeros(parameters.rows,parameters.cols,size(PPLHit,3),nIterrate);
nullDistMiss = zeros(parameters.rows,parameters.cols,size(PPLMiss,3),nIterrate);
for j=1:nIterrate
    randIndex = randperm(nTot);
    xgpHitRand = xgpComb(randIndex(1:nHit));
    xgpMissRand = xgpComb(randIndex(nHit+1:end));
    nullDistHit(:,:,:,j) = getPPL(xgpHitRand,parameters);
    nullDistMiss(:,:,:,j) = getPPL(xgpMissRand,parameters);
    j
end
muHit = mean(nullDistHit,4); % Mean of the null distribution
sigmaHit = std(nullDistHit,0,4); % Standard deviation of null distribution
muMiss = mean(nullDistMiss,4); % Mean of the null distribution
sigmaMiss = std(nullDistMiss,0,4); % Standard deviation of null distribution

PPLHitz = (PPLHit-muHit)./sigmaHit;
PPLMissz = (PPLMiss-muMiss)./sigmaMiss;

figure();
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PPLHitz,[1 2])),'-r','LineWidth',1.2); hold on;
% plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PAMiss,[1 2])),'-k','LineWidth',1);
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PPLMissz,[1 2])),'-k','LineWidth',1);
ylabel("z-score"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
title('z-scored Phase Alignment for Hits');box off; legend('Hits','Miss');

figure();
subplot(2,1,1);
title("Percentage Phase across all electrodes - Hits")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:32,reshape(PPLHitz,[],size(PPLHitz,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');
subplot(2,1,2);
title("Percentage Phase across all electrodes - Misses")
imagesc(IntanBehaviour.cueMissTrace(1).time,1:32,reshape(PPLMissz,[],size(PPLMissz,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');

nIterrate = 200;
xgpHit = arrayfun(@(s) s.xgp, IntanBehaviour.hitTrace, 'UniformOutput', false);
xgpMiss = arrayfun(@(s) s.xgp, IntanBehaviour.missTrace, 'UniformOutput', false);
xgpComb = [xgpHit xgpMiss];
nHit = size(xgpHit,2);
nMiss = size(xgpMiss,2);
nTot = nHit + nMiss;
nullDistHit = zeros(parameters.rows,parameters.cols,size(PPLHitReward,3),nIterrate);
nullDistMiss = zeros(parameters.rows,parameters.cols,size(PPLFA,3),nIterrate);
for j=1:nIterrate
    randIndex = randperm(nTot);
    xgpHitRand = xgpComb(randIndex(1:nHit));
    xgpMissRand = xgpComb(randIndex(nHit+1:end));
    nullDistHit(:,:,:,j) = getPPL(xgpHitRand,parameters);
    nullDistMiss(:,:,:,j) = getPPL(xgpMissRand,parameters);
    j
end
muHitReward = mean(nullDistHit,4); % Mean of the null distribution
sigmaHitReward = std(nullDistHit,0,4); % Standard deviation of null distribution
muFA = mean(nullDistMiss,4); % Mean of the null distribution
sigmaFA = std(nullDistMiss,0,4); % Standard deviation of null distribution

PPLHitRewardz = (PPLHit-muHitReward)./sigmaHitReward;
PPLFAz = (PPLMiss-muFA)./sigmaFA;


%% Average LFP for Hits and Misses 
LFPHit = arrayfun(@(s) reshape(s.xf,[],size(s.xf,3)), IntanBehaviour.cueHitTrace,"UniformOutput",false);
LFPMiss = arrayfun(@(s) reshape(s.xf,[],size(s.xf,3)), IntanBehaviour.cueMissTrace,"UniformOutput",false);
LFPHitReward = arrayfun(@(s) reshape(s.xf,[],size(s.xf,3)), IntanBehaviour.hitTrace,"UniformOutput",false);
LFPFA = arrayfun(@(s) reshape(s.xf,[],size(s.xf,3)), IntanBehaviour.missTrace,"UniformOutput",false);
LFPMIHit = arrayfun(@(s) reshape(s.xf,[],size(s.xf,3)), IntanBehaviour.MIHitTrace,"UniformOutput",false);
LFPMIFA = arrayfun(@(s) reshape(s.xf,[],size(s.xf,3)), IntanBehaviour.MIFATrace,"UniformOutput",false);
avgLFPHit = zeros(parameters.rows*parameters.cols,size(LFPHit{1,1},2));
avgLFPMiss = zeros(parameters.rows*parameters.cols,size(LFPMiss{1,1},2));
avgLFPHitReward = zeros(parameters.rows*parameters.cols,size(LFPHitReward{1,1},2));
avgLFPFA = zeros(parameters.rows*parameters.cols,size(LFPFA{1,1},2));
avgLFPMIHit = zeros(parameters.rows*parameters.cols,size(LFPMIHit{1,1},2));
avgLFPMIFA = zeros(parameters.rows*parameters.cols,size(LFPMIFA{1,1},2));

% Getting trial averaged LFP for each channel for Hits
a = cell2struct(LFPHit,'lfp',1);
for i=1:(parameters.rows*parameters.cols)
    avgLFPHit(i,:) = mean(cell2mat(arrayfun(@(s) s.lfp(i,:),a, 'UniformOutput',false)),1,'omitnan');
end

% Getting trial averaged LFP for each channel for Miss
a = cell2struct(LFPMiss,'lfp',1);
for i=1:(parameters.rows*parameters.cols)
    avgLFPMiss(i,:) = mean(cell2mat(arrayfun(@(s) s.lfp(i,:),a, 'UniformOutput',false)),1,'omitnan');
end

% Getting trial averaged LFP for each channel for Hits Rewards
a = cell2struct(LFPHitReward,'lfp',1);
for i=1:(parameters.rows*parameters.cols)
    avgLFPHitReward(i,:) = mean(cell2mat(arrayfun(@(s) s.lfp(i,:),a, 'UniformOutput',false)),1,'omitnan');
end

% Getting trial averaged LFP for each channel for False alarms
a = cell2struct(LFPFA,'lfp',1);
for i=1:(parameters.rows*parameters.cols)
    avgLFPFA(i,:) = mean(cell2mat(arrayfun(@(s) s.lfp(i,:),a, 'UniformOutput',false)),1,'omitnan');
end

% Getting trial averaged LFP for each channel for Hits MI
a = cell2struct(LFPMIHit,'lfp',1);
for i=1:(parameters.rows*parameters.cols)
    avgLFPMIHit(i,:) = mean(cell2mat(arrayfun(@(s) s.lfp(i,:),a, 'UniformOutput',false)),1,'omitnan');
end

% Getting trial averaged LFP for each channel for False alarms MI
a = cell2struct(LFPMIFA,'lfp',1);
for i=1:(parameters.rows*parameters.cols)
    avgLFPMIFA(i,:) = mean(cell2mat(arrayfun(@(s) s.lfp(i,:),a, 'UniformOutput',false)),1,'omitnan');
end

figure(); % Top half is hits and bottom half is misses
title("Trial Average LFP for Hits and Misses")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:64,[avgLFPHit;avgLFPMiss]); colormap(jet);
ylabel("Electrodes");xlabel("Time (s)"); 
h = colorbar; h.Label.String = 'Amplitude (uV)';
xline(0,'-k','Cue','LabelVerticalAlignment','top');
yline(32.5,'-k');caxis([-30 30]);
% yline(64.5,'-k');caxis([-20 20]);

figure(); % Top half is hits and bottom half is FA
title("Trial Average LFP for Hits and FA")
imagesc(IntanBehaviour.hitTrace(1).time,1:64,[avgLFPHitReward;avgLFPFA]); colormap(jet);
ylabel("Electrodes");xlabel("Time (s)"); 
h = colorbar; h.Label.String = 'Amplitude (uV)';
xline(0,'-k','Reward','LabelVerticalAlignment','top');
yline(32.5,'-k');caxis([-30 30]);
% yline(64.5,'-k');caxis([-20 20]);


figure(); % Top half is hits MI and bottom half is FA MI
title("Trial Average LFP for Hits and FA- Motion Initiation")
imagesc(IntanBehaviour.MIHitTrace(1).time,1:64,[avgLFPMIHit;avgLFPMIFA]); colormap(jet);
ylabel("Electrodes");xlabel("Time (s)"); 
h = colorbar; h.Label.String = 'Amplitude (uV)';
xline(0,'-k','MI','LabelVerticalAlignment','top');
yline(32.5,'-k');caxis([-30 30]);
% yline(64.5,'-k');caxis([-20 20]);
%% Average PGD accross frequency bands
[PGDfreqHit,PGDfreqMiss] = getPGDFreqBand(LFP.LFPdatacube,IntanBehaviour,0,parameters);

figure();
PGDFreq = [5:5:100];
plot(PGDFreq,PGDfreqHit); hold on;
plot(PGDFreq,PGDfreqMiss);
xlabel('Frequency'); ylabel('Phase Gradient Directionality'); ylim([0.35 0.7]); box off;

%% Beta Burst detection using Hilbert Amplitude 
[BetaEvent] = detectBetaEvent(LFP.xgpbetamean,IntanBehaviour.cueHitTrace,parameters);

%% Beta event detection 
avgBetaband = mean(LFP.xfbeta,1);
window = Encoder.trialTime(:,3:4);
windowStop = Encoder.trialTimeStop(:,3:4);

betatrials = zeros(Encoder.nTrig,LFP.Fs*(windowAfterTrig+windowBeforeTrig));
betatrialsStop = zeros(Encoder.nTrigStop,LFP.Fs*(windowAfterTrig+windowBeforeTrig));
for i=1:Encoder.nTrig
   betatrials(i,:) = avgBetaband(window(i,1):window(i,2));
end
for i=1:Encoder.nTrigStop
   betatrialsStop(i,:) = avgBetaband(windowStop(i,1):windowStop(i,2));
end
avgbetaGroup = groupBetaBurstDetection(LFP,betatrials',window,LFP.Fs);
avgbetaGroupStop = groupBetaBurstDetection(LFP,betatrialsStop',windowStop,LFP.Fs);

% Calculating for each electrode
betatrials = zeros(Encoder.nTrig,LFP.Fs*(windowAfterTrig+windowBeforeTrig));
betatrialsStop = zeros(Encoder.nTrigStop,LFP.Fs*(windowAfterTrig+windowBeforeTrig));
for i=1:parameters.rows*parameters.cols
    for j=1:Encoder.nTrig
        betatrials(j,:) = LFP.beta_band(i,window(j,1):window(j,2));
    end
    for j=1:Encoder.nTrigStop
        betatrialsStop(j,:) = LFP.beta_band(i,windowStop(j,1):windowStop(j,2));
    end
    betaEvents(i).betagroup = groupBetaBurstDetection(LFP,betatrials',window,LFP.Fs);
    betaEventsStop(i).betagroup = groupBetaBurstDetection(LFP,betatrialsStop',windowStop,LFP.Fs);
end

% Number of beta events per trial
nbetaevents = zeros(Encoder.nTrig,1);
nbetaeventsStop = zeros(Encoder.nTrigStop,1);
for i=1:parameters.rows*parameters.cols
    nbetaevents = nbetaevents + betaEvents(i).betagroup.betaBurst.NumDetectedBeta;
    nbetaeventsStop = nbetaeventsStop + betaEventsStop(i).betagroup.betaBurst.NumDetectedBeta;
end
figure,plot(nbetaevents);hold on;plot(nbetaeventsStop);

% Number of beta events on each electrode 
nbetaperElec = zeros(parameters.cols,parameters.rows);
for i=1:parameters.rows*parameters.cols
    a = chToGrid(i,parameters);
    nbetaperElec(a(1),a(2)) = mean(betaEvents(i).betagroup.betaBurst.NumDetectedBeta,"all");
end

figure(); imagesc(nbetaperElec');%set(gca,'YDir','normal');
title('Spatial map of beta events accross all trials'); colorbar;

%% PLotting LFP 
figure();
trialno = 1;
for i=1:32
    subplot(parameters.rows,parameters.cols,i);
    if ismember(i,Intan.badChMap), continue; end
    plot(LFP.times(Encoder.trialTime(trialno,3):Encoder.trialTime(trialno,4)),squeeze(goodTrial.xf(floor((i-1)/parameters.cols)+1,mod(i-1,parameters.cols)+1,:))');xline(0,'-r');
end
