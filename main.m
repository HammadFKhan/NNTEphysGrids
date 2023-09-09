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

%% Loading previous workspace
% loadPreviousWS

%%  PreProcessing
load GridsLowDenNeedle_chanmap.mat;  % load the channel map for the IntanConcatenate function
parameters.rows = 8;  % Number of rows of electrodes on the Grid
parameters.cols = 4;  % Number of colums of electrodes on the Grid
parameters.Fs = 1000;
parameters.ts = 1/parameters.Fs;
parameters.windowBeforePull = 1; % in seconds
parameters.windowAfterPull = 1; % in seconds
parameters.windowBeforeCue = 0.5; % in seconds 
parameters.windowAfterCue = 1.5; % in seconds 
parameters.experiment = 'cue'; % self - internally generated, cue - cue initiated 
parameters.opto = 1; % 1 - opto ON , 0 - opto OFF

IntanConcatenate
fpath = Intan.path; % where on disk do you want the analysis? ideally and SSD...

% Generating time series from Intan data
Ts = 1/Intan.offsetSample;
Intan.Tmax = Ts * size(Intan.allIntan,2);
Intan.t = 0:Ts:Intan.Tmax-Ts;

%% Removing bad channels from impedance values
[Z,Intan.goodChMap,Intan.badChMap] = readImp(electrode_map,100e6);
figure('Name','Impedance Test at 1kHz');boxchart(Z); xlabel('n = ' + string(size(Z,1)));ylabel('Impedance (in \Omega)');set(gca,'xticklabel',{[]})
%Intan = removeBadCh(Intan,Intan.badCh);

%% LFP
set(0,'DefaultFigureWindowStyle','normal')
LFP = fastpreprocess_filtering(Intan.allIntan,5000);
% LFP = bestLFP(LFP);
% LFP = bandFilter(LFP,'depth'); % Extract LFPs based on 'depth' or 'single'
% LFPplot(LFP);
LFP = createDataCube(LFP,parameters.rows,parameters.cols,Intan.goodChMap); % Creating datacube

%% Loading Lever Data 
[Behaviour] = readLever(parameters,LFP.times);
figure('Name','Lever Trace');plot(Behaviour.time,Behaviour.leverTrace,'LineWidth',1.5);ylim([-5 50]);xlabel('Time (in s)');ylabel('Lever Position in mV');yline(23);
xline(squeeze(Behaviour.hit(:,2)),'-.b',cellstr(num2str((1:1:Behaviour.nHit)')),'LabelVerticalAlignment','top');
xline(squeeze(Behaviour.miss(:,2)),'-.r',cellstr(num2str((1:1:Behaviour.nMiss)')),'LabelVerticalAlignment','bottom'); xlim([410 470]); box off;

% Plotting Lever traces for Hits and Misses
figure('Name','Average Lever Traces for Hits & Misses');
subplot(1,2,1);
for i=1:Behaviour.nHit
    plot(Behaviour.hitTrace(i).time-1,Behaviour.hitTrace(i).trace,'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
end
plot(Behaviour.hitTrace(1).time-parameters.windowBeforePull,mean(horzcat(Behaviour.hitTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
yline(10,'--.b','Threshold','LabelHorizontalAlignment','left'); 
ylabel('Lever deflection (in mV)');xlabel('Time (in s)');title('Average Lever Traces for Hits');box off;

subplot(1,2,2);
for i=1:Behaviour.nMiss
    plot(Behaviour.missTrace(i).time-1,Behaviour.missTrace(i).trace,'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
end
plot(Behaviour.missTrace(1).time-parameters.windowBeforePull,mean(horzcat(Behaviour.missTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
yline(10,'--.b','Threshold','LabelHorizontalAlignment','left'); 
ylabel('Lever deflection (in mV)');xlabel('Time (in s)');title('Average Lever Traces for Misses');box off;

% Plotting Lever traces for Cue Hits and Cue miss 
figure('Name','Average Lever Traces for Cue Hits');
for i=1:Behaviour.nCueHit
    plot(Behaviour.cueHitTrace(i).time-parameters.windowBeforeCue,Behaviour.cueHitTrace(i).trace,'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
end
plot(Behaviour.cueHitTrace(1).time-parameters.windowBeforeCue,mean(horzcat(Behaviour.cueHitTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
yline(10,'--.b','Threshold','LabelHorizontalAlignment','left'); 
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(Behaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
ylabel('Lever deflection (in mV)');xlabel('Time (in s)');title('Average Lever Traces for Cue Hits');box off;

%% Reading behaviour data from Intan traces 
IntanBehaviour = readLeverIntan(parameters,LFP.times,Intan.analog_adc_data,Intan.dig_in_data,Behaviour);

% Plotting Lever traces for Hits and Misses
figure('Name','Average Lever Traces for Hits & Misses');
subplot(1,2,1);
for i=1:size(IntanBehaviour.cueHitTrace,2)
    plot(IntanBehaviour.hitTrace(i).time,IntanBehaviour.hitTrace(i).trace,'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
end
plot(IntanBehaviour.hitTrace(1).time,mean(horzcat(IntanBehaviour.hitTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
yline(IntanBehaviour.threshold,'--.b','Threshold','LabelHorizontalAlignment','left'); 
ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for Hits');box off;
subplot(1,2,2);
for i=1:IntanBehaviour.nCueMiss
    plot(IntanBehaviour.missTrace(i).time,IntanBehaviour.missTrace(i).trace,'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
end
plot(IntanBehaviour.missTrace(1).time,mean(horzcat(IntanBehaviour.missTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
yline(IntanBehaviour.threshold,'--.b','Threshold','LabelHorizontalAlignment','left'); 
ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for Misses');box off;

% Plotting Lever traces for Cue Hit 
figure('Name','Average Lever Traces for Cue Hits');
for i=1:size(IntanBehaviour.cueHitTrace,2)
    plot(IntanBehaviour.cueHitTrace(i).time,IntanBehaviour.cueHitTrace(i).trace,'Color',[0 0 0 0.1],'LineWidth',1.5);
    hold on;
end
plot(IntanBehaviour.cueHitTrace(1).time,mean(horzcat(IntanBehaviour.cueHitTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
yline(IntanBehaviour.threshold,'--.b','Threshold','LabelHorizontalAlignment','left'); 
xline(0,'--r','Cue','LabelVerticalAlignment','top');ylim([0 0.1]);
xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for Cue Hits');box off;

% Plotting Lever traces for Cue Miss 
figure('Name','Average Lever Traces for Cue Misses');
for i=1:size(IntanBehaviour.cueMissTrace,2)
    plot(IntanBehaviour.cueMissTrace(i).time,IntanBehaviour.cueMissTrace(i).trace,'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
end
plot(IntanBehaviour.cueMissTrace(1).time,mean(horzcat(IntanBehaviour.cueMissTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
yline(IntanBehaviour.threshold,'--.b','Threshold','LabelHorizontalAlignment','left'); 
xline(0,'--r','Cue','LabelVerticalAlignment','top');
ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for Cue Misses');box off;

%% Generalized Phase 
LFP.xf = bandpass_filter(LFP.LFPdatacube,5,40,4,1000);
[LFP.xgp, LFP.wt] = generalized_phase(LFP.xf,1000,0);
LFP.xfbeta = bandpass_filter(LFP.LFPdatacube,10,30,4,1000);
[LFP.xgpbeta, LFP.wtbeta] = generalized_phase(LFP.xfbeta,1000,0);
LFP.xftheta = bandpass_filter(LFP.LFPdatacube,4,10,4,1000);
[LFP.xgptheta, LFP.wttheta]  = generalized_phase(LFP.xftheta,1000,0);
LFP.xfgamma = bandpass_filter(LFP.LFPdatacube,30,80,4,1000);
[LFP.xgpgamma, LFP.wtgamma]  = generalized_phase(LFP.xfgamma,1000,0);
[parameters.X,parameters.Y] = meshgrid( 1:parameters.cols, 1:parameters.rows );
LFP.xfwide = bandpass_filter(LFP.LFPdatacube,5,90,4,1000);
LFP.xfbetanarrow = bandpass_filter(LFP.LFPdatacube,6,9,4,1000);
[LFP.xgpbetanarrow, LFP.wtbetanarrow] = generalized_phase(LFP.xfbetanarrow,1000,0);
% GP for spatial mean LFP 
[LFP.xgpbetamean, ~] = generalized_phase(mean(LFP.xfbetanarrow,[1,2]),1000,0);

%% Add trial segmented data to IntanBehaviour Variable
IntanBehaviour = addLFPToBehaviour(IntanBehaviour,LFP);
% Saving paramters, path, IntanBehaviour to bin file 
savepath = uigetdir(path);
sessionName = [savepath,'/','Day2.mat'];
% save(sessionName,"IntanBehaviour","fpath","parameters","-v7.3");
save(sessionName,"IntanBehaviour","fpath","parameters","Waves","betaWaves","thetaWaves","gammaWaves","-v7.3");

%% Combining  multiple Intanbehaviour structs from multiple sessions
combIntanBehaviour = horzcat(IntanBehaviour1, IntanBehaviour2);
IntanBehaviour.cueHitTrace = horzcat(combIntanBehaviour(1:end).cueHitTrace);
IntanBehaviour.cueMissTrace = horzcat(combIntanBehaviour(1:end).cueMissTrace);
IntanBehaviour.missTrace = horzcat(combIntanBehaviour(1:end).missTrace);
IntanBehaviour.reactionTime = horzcat(combIntanBehaviour(1:end).reactionTime);
clear combIntanBehaviour IntanBehaviour1 IntanBehaviour2;
%% Power Spectrum during task across channels 
[PSDChHit , f] = getAvgPSD(IntanBehaviour.cueHitTrace,parameters);
[PSDChMiss , f] = getAvgPSD(IntanBehaviour.cueMissTrace,parameters);

avgPSDHit = squeeze(10*log10(mean(PSDChHit,[1 2])));
trialPSDHit = squeeze(10*log10(mean(PSDChHit,2)));
avgPSDMiss = squeeze(10*log10(mean(PSDChMiss,[1 2])));
trialPSDMiss = squeeze(10*log10(mean(PSDChMiss,2)));

figure();
subplot(1,2,1);
plot(f(1:51),trialPSDHit(:,1:51),'Color', [0 0 1 0.1]);
hold on;
plot(f(1:51),avgPSDHit(1:51),'Color', [0 0 1 1],'LineWidth',1.5);
ylim([0 50]);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');
title('Average Power Spectral Density for HitTrials');
box off;

subplot(1,2,2);
plot(f(1:51),trialPSDMiss(:,1:51),'Color', [1 0 0 0.1]);
hold on;
plot(f(1:51),avgPSDMiss(1:51),'Color', [1 0 0 1],'LineWidth',1.5);
ylim([0 50]);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');
title('Average Power Spectral Density for MissTrials');
box off;

trialAvgPSD = 10*log10(squeeze(mean(PSDChHit,1)));
figure();
imagesc(trialAvgPSD(:,1:51)');
colormap("jet");set(gca,'YDir','normal');

%% Wavelet spectrogram
[hitAvgSpectrogram, hitSpectrogramCWT,AvgHitTrace ,fwt] = getAvgSpectogram(IntanBehaviour.cueHitTrace,parameters,[5 40]);
[missAvgSpectrogram, missSpectrogramCWT,AvgMissTrace,fwt] = getAvgSpectogram(IntanBehaviour.cueMissTrace,parameters,[5 40]);

% Global average spectogram
figure('Name','Trial Averaged Wavelet Spectrogram for Hits & Misses');
subplot(1,2,1);
plotSpectrogram(10*log10((squeeze(hitAvgSpectrogram))),IntanBehaviour.cueHitTrace(1).time,fwt,'surf','Wavelet Based Spectrogram for Hits','Time (s)','Frequency (Hz)')
caxis([-5 25]);hold on; yyaxis right; box off;
plot(IntanBehaviour.cueHitTrace(1).time,AvgHitTrace,'-w','LineWidth',2.5);
ylabel('Lever deflection (mV)'); ylim([0 0.1]);
subplot(1,2,2);
plotSpectrogram(10*log10((squeeze(missAvgSpectrogram))),IntanBehaviour.cueMissTrace(1).time,fwt,'surf','Wavelet Based Spectrogram for Misses','Time (s)','Frequency (Hz)')
caxis([-5 25]);hold on; yyaxis right; box off;
plot(IntanBehaviour.cueMissTrace(1).time,AvgMissTrace,'-w','LineWidth',2.5);
ylabel('Lever deflection (mV)'); ylim([0 0.1]); box off;

% Movement average spectogram 
allAvgSpectogram = mean(cat(1,hitSpectrogramCWT,missSpectrogramCWT));
allAvgLeverTrace = mean(cat(2,horzcat(IntanBehaviour.hitTrace(1:end).trace),horzcat(IntanBehaviour.missTrace(1:end).trace)),2);
figure('Name','Trial Averaged Wavelet Spectrogram for Lever Pull');
plotSpectrogram(squeeze(allAvgSpectogram),IntanBehaviour.hitTrace(1).time-parameters.windowBeforePull,fwt,'Wavelet Based Spectrogram for Lever Pull','Time (s)','Frequency (Hz)');
hold on; yyaxis right; box off;
plot(IntanBehaviour.hitTrace(1).time-parameters.windowBeforePull,allAvgLeverTrace,'-w','LineWidth',2.5);
ylabel('Lever deflection (mV)'); box off;

% Plotting for specific trial
trialno = 48;
figure('Name','Spatial Averaged Wavelet Spectrogram for Hits & Misses');
subplot(2,1,1);
plotSpectrogram(squeeze(hitSpectrogramCWT(trialno,:,:)),IntanBehaviour.hitTrace(trialno).time-parameters.windowBeforePull,fwt,'Wavelet based Average Spectogram for Hits','Time (s)','Frequency (Hz)');
hold on; yyaxis right; box off;
plot(IntanBehaviour.hitTrace(trialno).time-parameters.windowBeforePull,IntanBehaviour.hitTrace(trialno).trace,'-w','LineWidth',2.5);
ylabel('Lever deflection (mV)'); 
subplot(2,1,2);
trialno = 11;
plotSpectrogram(squeeze(missSpectrogramCWT(trialno,:,:)),IntanBehaviour.missTrace(trialno).time-parameters.windowBeforePull,fwt,'Wavelet based Average Spectogram for Misses','Time (s)','Frequency (Hz)');
hold on; yyaxis right; box off;
plot(IntanBehaviour.missTrace(trialno).time-parameters.windowBeforePull,IntanBehaviour.missTrace(trialno).trace,'-w','LineWidth',2.5);
ylabel('Lever deflection (mV)'); box off;

trialno = 1;
figure('Name','Spatial Averaged Wavelet Spectrogram for Hits & Misses');
subplot(2,1,1);
plotSpectrogram((squeeze(hitSpectrogramCWT(trialno,:,:))),IntanBehaviour.cueHitTrace(trialno).time,fwt,'Wavelet based Average Spectogram for Hits','Time (s)','Frequency (Hz)');
hold on; yyaxis right; box off;
plot(IntanBehaviour.cueHitTrace(trialno).time,squeeze(mean(LFP.xf(:,:,IntanBehaviour.cueHitTrace(trialno).LFPIndex),[1 2])),'-w','LineWidth',1.5);
plot(IntanBehaviour.cueHitTrace(trialno).time,squeeze(abs(mean(LFP.xgp(:,:,IntanBehaviour.cueHitTrace(trialno).LFPIndex),[1 2]))),'--w','LineWidth',1);
ylabel('Amplitude (\mu V)'); 
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(IntanBehaviour.reactionTime(trialno),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
subplot(2,1,2);
trialno = 11;
plotSpectrogram(squeeze(missSpectrogramCWT(trialno,:,:)),IntanBehaviour.cueMissTrace(trialno).time,fwt,'Wavelet based Average Spectogram for Misses','Time (s)','Frequency (Hz)');
hold on; yyaxis right; box off;
plot(IntanBehaviour.cueMissTrace(trialno).time,squeeze(mean(LFP.xf(:,:,IntanBehaviour.cueMissTrace(trialno).LFPIndex),[1 2])),'-w','LineWidth',1.5);
plot(IntanBehaviour.cueMissTrace(trialno).time,squeeze(abs(mean(LFP.xgp(:,:,IntanBehaviour.cueMissTrace(trialno).LFPIndex),[1 2]))),'--w','LineWidth',1);
xline(0,'--r','Cue','LabelVerticalAlignment','top');
ylabel('Amplitude (\mu V)');  box off;


%% Initializing plotting options
options.subject = 'W'; % this can be 'W' or 'T' (two marmoset subjects)
options.plot = true; % this option turns plots ON or OFF
options.plot_shuffled_examples = false; % example plots w/channels shuffled in space

%% Wave detection in velocity triggered windows
parameters.spacing = 0.1; % Grid spacing in mm
nShuffle = 1000;
threshold = 99;
trialno = 55;

% Wave detection for wide band
disp('Wave Detection for wide band ...')
xf = arrayfun(@(s) s.xf, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
xgp = arrayfun(@(s) s.xgp, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
wt = arrayfun(@(s) s.wt, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
% parameters.rhoThres = getRhoThreshold(xgp,IntanBehaviour.cueHitTrace,parameters,nShuffle,trialno,threshold);
Waves.wavesHit = detectWaves(xf,xgp,wt,IntanBehaviour.cueHitTrace,parameters,parameters.rhoThres);
Waves.wavesHit = detectPlanarWaves(xf,xgp,wt,IntanBehaviour.cueHitTrace,parameters,0.5);
if isfield(IntanBehaviour,'cueMissTrace')
    xf = arrayfun(@(s) s.xf, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgp, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wt, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    Waves.wavesMiss = detectWaves(xf,xgp,wt,IntanBehaviour.cueMissTrace,parameters,parameters.rhoThres);
    Waves.wavesMiss = detectPlanarWaves(xf,xgp,wt,IntanBehaviour.cueMissTrace,parameters,0.5);
end
if isfield(IntanBehaviour,'missTrace')
    xf = arrayfun(@(s) s.xf, IntanBehaviour.missTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgp, IntanBehaviour.missTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wt, IntanBehaviour.missTrace, 'UniformOutput', false);
    Waves.wavesFA = detectWaves(xf,xgp,wt,IntanBehaviour.missTrace,parameters,parameters.rhoThres);
end

% Wave detection for theta band
disp('Wave Detection for theta band ...')
threshold = 99;
xf = arrayfun(@(s) s.xftheta, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
xgp = arrayfun(@(s) s.xgptheta, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
wt = arrayfun(@(s) s.wttheta, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
parameters.thetarhoThres = getRhoThreshold(xgp,IntanBehaviour.cueHitTrace,parameters,nShuffle,trialno,threshold);
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
Wavesall = detectWaves(LFP.xf,LFP.xgp,LFP.wt,allwaves,parameters);
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
[WaveStats(1)] = getWaveStats(Waves,parameters,plotOption);
[WaveStats(2)] = getWaveStats(thetaWaves,parameters,plotOption);
[WaveStats(3)] = getWaveStats(betaWaves,parameters,plotOption);
[WaveStats(4)] = getWaveStats(gammaWaves,parameters,plotOption);

[WaveStats2(1)] = getInitRewardStats(Waves,parameters,plotOption);
[WaveStats2(2)] = getInitRewardStats(thetaWaves,parameters,plotOption);
[WaveStats2(3)] = getInitRewardStats(betaWaves,parameters,plotOption);
[WaveStats2(4)] = getInitRewardStats(gammaWaves,parameters,plotOption);

%% Mutual Information
% For phase
xgpHit = arrayfun(@(s) angle(s.xgp), IntanBehaviour.cueHitTrace, 'UniformOutput', false);
xgpMiss = arrayfun(@(s) angle(s.xgp), IntanBehaviour.cueMissTrace, 'UniformOutput', false);
[MIPhase] = getMutualInformation(xgpHit,xgpMiss,parameters);

figure();
title("Mututal Information across all electrodes - Phase")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:32,peakSort2DArray(reshape(MIPhase,[],size(MIPhase,3)),'descend',2)); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)"); 
h = colorbar; h.Label.String = 'Information (bits)';
xline(0,'-w','Cue','LabelVerticalAlignment','top');

% For amplitude
xgpHit = arrayfun(@(s) abs(s.xgp), IntanBehaviour.cueHitTrace, 'UniformOutput', false);
xgpMiss = arrayfun(@(s) abs(s.xgp), IntanBehaviour.cueMissTrace, 'UniformOutput', false);

[MIAmp] = getMutualInformation(xgpHit,xgpMiss,parameters);

figure();
title("Mututal Information across all electrodes - Amplitude")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:32,peakSort2DArray(reshape(MIAmp,[],size(MIAmp,3)),'descend',2)); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)"); 
h = colorbar; h.Label.String = 'Information (bits)';
xline(0,'-w','Cue','LabelVerticalAlignment','top');


%% Average PGD 
avgPGDHit = mean(vertcat(Waves.wavesHit.PGD),1);
avgPGDMiss = mean(vertcat(Waves.wavesMiss.PGD),1);
avgPGDFA = mean(vertcat(Waves.wavesFA.PGD),1);

figure(); hold on;
plot(IntanBehaviour.cueHitTrace(1).time,avgPGDHit,'-r','LineWidth',1.2); hold on;
plot(IntanBehaviour.cueHitTrace(1).time,avgPGDMiss,'-k','LineWidth',1);
plot(IntanBehaviour.cueHitTrace(1).time,avgPGDFA,'-b','LineWidth',1);
ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD)');box off;  legend('Hits','Misses','False Alarms');

xgp = arrayfun(@(s) s.xgp, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
[avgPGDHit,avgPGDHitNull] = getAvgPGD(xgp,Waves.wavesHit,IntanBehaviour.cueHitTrace,1,parameters);

figure(); hold on;
plot(IntanBehaviour.cueHitTrace(1).time,avgPGDHit,'-k','LineWidth',1.2); hold on;
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(avgPGDHitNull),'-r','LineWidth',1.2);
ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD)');box off;  legend('Hits','Shuffled');

%% Cross-Trial Phase Alignment
xgp = arrayfun(@(s) s.xgp, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
PAHit = getPhaseAlignment(xgp,parameters);
xgp = arrayfun(@(s) s.xgp, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
PAMiss = getPhaseAlignment(xgp,parameters);
xgp = arrayfun(@(s) s.xgp, IntanBehaviour.missTrace, 'UniformOutput', false);
PAFA = getPhaseAlignment(xgp,parameters);


figure();
subplot(3,1,1);
title("Phase Alignment across all electrodes - Hits")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:32,reshape(PAHit,[],size(PAHit,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');
subplot(3,1,2);
title("Phase Alignment across all electrodes - Misses")
imagesc(IntanBehaviour.cueMissTrace(1).time,1:32,reshape(PAMiss,[],size(PAMiss,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');
subplot(3,1,3);
title("Phase Alignment across all electrodes - False Alarms")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:32,reshape(PAFA,[],size(PAFA,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');

figure();
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PAHit,[1 2])),'-r','LineWidth',1.2); hold on;
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PAMiss,[1 2])),'-k','LineWidth',1);
% plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PAFA,[1 2])),'-k','LineWidth',1);
ylabel("Phase Alignment"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
title('Phase Alignment for Hits');box off;% legend('Hits','Miss');
% z-scoring 
PAHit1 = reshape(PAHit,[],size(PAHit,3));
PAMiss1 = reshape(PAMiss,[],size(PAMiss,3));
nIterrate = 5000;
xgpHit = arrayfun(@(s) s.xgp, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
xgpMiss = arrayfun(@(s) s.xgp, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
xgpComb = [xgpHit xgpMiss];
nHit = size(xgpHit,2);
nMiss = size(xgpMiss,2);
nTot = nHit + nMiss;

nullDistHit = zeros(parameters.rows,parameters.cols,size(PAHit,3),nIterrate);
nullDistMiss = zeros(parameters.rows,parameters.cols,size(PAMiss,3),nIterrate);
for j=1:nIterrate
    randIndex = randperm(nTot);
    xgpHitRand = xgpComb(randIndex(1:nHit));
    xgpMissRand = xgpComb(randIndex(nHit+1:end));
    nullDistHit(:,:,:,j) = getPhaseAlignment(xgpHitRand,parameters);
    nullDistMiss(:,:,:,j) = getPhaseAlignment(xgpMissRand,parameters);
    j
end
muHit = mean(nullDistHit,4); % Mean of the null distribution
sigmaHit = std(nullDistHit,0,4); % Standard deviation of null distribution
muMiss = mean(nullDistMiss,4); % Mean of the null distribution
sigmaMiss = std(nullDistMiss,0,4); % Standard deviation of null distribution

PAHitz = (PAHit-muHit)./sigmaHit;
PAMissz = (PAMiss-muMiss)./sigmaMiss;

figure();
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PAHitz,[1 2])),'-r','LineWidth',1.2); hold on;
% plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PAMiss,[1 2])),'-k','LineWidth',1);
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PAMissz,[1 2])),'-k','LineWidth',1);
ylabel("z-score"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
title('z-scored Phase Alignment for Hits');box off; legend('Hits','Miss');


figure();
subplot(2,1,1);
title("Phase Alignment across all electrodes - Hits")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:32,reshape(PAHitz,[],size(PAHitz,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');caxis([-4 8]);
subplot(2,1,2);
title("Phase Alignment across all electrodes - Misses")
imagesc(IntanBehaviour.cueMissTrace(1).time,1:32,reshape(PAMissz,[],size(PAMissz,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');caxis([-4 8]);

%% Percent Phase Locking
xgp = arrayfun(@(s) s.xgp, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
[PPLHit] = getPPL(xgp,parameters);
xgp = arrayfun(@(s) shuffle3DMatrix(s.xgp,3), IntanBehaviour.cueHitTrace, 'UniformOutput', false);
[PPLHitShuffle] = getPPL(xgp,parameters);
xgp = arrayfun(@(s) s.xgp, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
[PPLMiss] = getPPL(xgp,parameters);
xgp = arrayfun(@(s) s.xgp, IntanBehaviour.missTrace, 'UniformOutput', false);
[PPLFA] = getPPL(xgp,parameters);

figure();
subplot(4,1,1);
title("Percentage Phase across all electrodes - Hits")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:32,reshape(PPLHit,[],size(PPLHit,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');
subplot(4,1,2);
title("Percentage Phase across all electrodes - Misses")
imagesc(IntanBehaviour.cueMissTrace(1).time,1:32,reshape(PPLMiss,[],size(PPLMiss,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');
subplot(4,1,3);
title("Percentage Phase across all electrodes - False Alarms")
imagesc(IntanBehaviour.missTrace(1).time-parameters.windowBeforeCue,1:32,reshape(PPLFA,[],size(PPLFA,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');
subplot(4,1,4);
title("Percentage Phase across all electrodes - Hit Shuffled")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:32,reshape(PPLHitShuffle,[],size(PPLHitShuffle,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');


figure();
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PPLHit,[1 2])),'-r','LineWidth',1.2); hold on;
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PPLHitShuffle,[1 2])),'-k','LineWidth',1);
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PPLFA,[1 2])),'-b','LineWidth',1);
% plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PPLMiss,[1 2])),'-g','LineWidth',0.2);
ylabel("Percentage Phase Locking"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
title('Percentage Phase Locking for hits');box off;legend('Hits','Shuffled','False Alarms');%,'Misses');

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

%% Average PGD accross frequency bands
[PGDfreqHit,PGDfreqMiss] = getPGDFreqBand(LFP.LFPdatacube,IntanBehaviour,0,parameters);

figure();
PGDFreq = [5:5:100];
plot(PGDFreq,PGDfreqHit); hold on;
plot(PGDFreq,PGDfreqMiss);
xlabel('Frequency'); ylabel('Phase Gradient Directionality'); ylim([0.35 0.7]); box off;

%% Average LFP for Hits and Misses 
LFPHit = arrayfun(@(s) reshape(s.xf,[],size(s.xf,3)), IntanBehaviour.cueHitTrace,"UniformOutput",false);
LFPMiss = arrayfun(@(s) reshape(s.xf,[],size(s.xf,3)), IntanBehaviour.cueMissTrace,"UniformOutput",false);
LFPFA = arrayfun(@(s) reshape(s.xf,[],size(s.xf,3)), IntanBehaviour.missTrace,"UniformOutput",false);

avgLFPHit = zeros(parameters.rows*parameters.cols,size(LFPHit{1,1},2));
avgLFPMiss = zeros(parameters.rows*parameters.cols,size(LFPMiss{1,1},2));
avgLFPFA = zeros(parameters.rows*parameters.cols,size(LFPFA{1,1},2));

% Getting trial averaged LFP for each channel for Hits
a = cell2struct(LFPHit,'lfp',1);
for i=1:(parameters.rows*parameters.cols)
    avgLFPHit(i,:) = mean(cell2mat(arrayfun(@(s) s.lfp(i,:),a, 'UniformOutput',false)),1);
end

% Getting trial averaged LFP for each channel for Hits
a = cell2struct(LFPMiss,'lfp',1);
for i=1:(parameters.rows*parameters.cols)
    avgLFPMiss(i,:) = mean(cell2mat(arrayfun(@(s) s.lfp(i,:),a, 'UniformOutput',false)),1);
end

% Getting trial averaged LFP for each channel for False alarms
a = cell2struct(LFPFA,'lfp',1);
for i=1:(parameters.rows*parameters.cols)
    avgLFPFA(i,:) = mean(cell2mat(arrayfun(@(s) s.lfp(i,:),a, 'UniformOutput',false)),1);
end

figure(); % Top half is hits and bottom half is misses
title("Trial Average LFP for Hits and Misses")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:96,[avgLFPHit;avgLFPMiss;avgLFPFA]); colormap(jet);
ylabel("Electrodes");xlabel("Time (s)"); 
h = colorbar; h.Label.String = 'Amplitude (uV)';
xline(0,'-k','Cue','LabelVerticalAlignment','top');
yline(32.5,'-k');caxis([-20 20]);
yline(64.5,'-k');caxis([-20 20]);

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
