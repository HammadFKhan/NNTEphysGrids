function [wavesStat] = getWaveStats(Waves,parameters,plot)

rows = parameters.rows;
cols = parameters.cols;

% Number of waves
WavesPerTrial = mean( horzcat(Waves.wavesStart(1:end).nWaves),'all');
totalWaves =  sum( horzcat(Waves.wavesStart(1:end).nWaves));

WavesPerTrialStop = mean( horzcat(Waves.wavesStop(1:end).nWaves),'all');
totalWavesStop =  sum( horzcat(Waves.wavesStop(1:end).nWaves));

% Wave Speed stats
speedComb = horzcat(Waves.wavesStart(1:end).speed);
avgSpeed = mean(speedComb);

speedCombStop = horzcat(Waves.wavesStop(1:end).speed);
avgSpeedStop = mean(speedCombStop);

% Perform the t-test.
[p, t] = ranksum(speedComb, speedCombStop);
% Print the results.
disp('Wave Speed')
disp('h-statistic:');
disp(t);
disp('p-value:');
disp(p);


if plot == 1
    figure('Name','Histogram of wave speeds in motion Initiation and termination');
    subplot(2,1,1);
    histfit(speedComb,100,'kernel');
    xline(avgSpeed,'-r',{'Mean speed = ' num2str(avgSpeed) ' cm/s'});
    xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Motion Initiation');
    subplot(2,1,2);
    histfit(speedCombStop,100,'kernel');
    xline(avgSpeedStop,'-r',{'Mean speed = ' num2str(avgSpeedStop) ' cm/s'});
    xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Motion Termination');

    figure('Name','Wave speeds in motion Initiation and termination');
    group = [ones(size(speedComb')); 2.*ones(size(speedCombStop'))];
    boxplot([speedComb';speedCombStop'],group,'BoxStyle','filled','PlotStyle','compact');
    set(gca,'XTickLabel',{'Initiation','Termination'});
    ylabel('Wave speed in cm/s');
end


% Wavelength stats
lComb = horzcat(Waves.wavesStart(1:end).wavelength);
avgl = mean(lComb);

lCombStop = horzcat(Waves.wavesStop(1:end).wavelength);
avglStop = mean(lCombStop);

% Perform the t-test.
[p, t] = ranksum(lComb, lCombStop);
% Print the results.
disp('Wavelength')
disp('h-statistic:');
disp(t);
disp('p-value:');
disp(p);


if plot == 1
    figure('Name','Histogram of wavelength in motion Initiation and termination');
    subplot(2,1,1);
    histfit(lComb,100,'kernel');
    xline(avgl,'-r',{'Mean wavelength = ' num2str(avgl) ' cm'});
    xlabel('Wavelength in cm');ylabel('Frequency');title('Wavelength: Motion Initiation');
    subplot(2,1,2);
    histfit(lCombStop,100,'kernel');
    xline(avglStop,'-r',{'Mean wavelength = ' num2str(avglStop) ' cm'});
    xlabel('Wavelength in cm');ylabel('Frequency');title('Wavelength : Motion Termination');

    figure('Name','Wavelength in motion Initiation and termination');
    group = [ones(size(lComb')); 2.*ones(size(lCombStop'))];
    boxplot([lComb';lCombStop'],group,'BoxStyle','filled','PlotStyle','compact');
    set(gca,'XTickLabel',{'Initiation','Termination'});
    xlabel('Wavelength in cm');
end


% Wave direction stats
dirComb = horzcat(Waves.wavesStart(1:end).waveDir);
avgDir = mean(dirComb);
dirCombStop = horzcat(Waves.wavesStop(1:end).waveDir);
avgDirStop = mean(dirCombStop);

[p, t] = ranksum(dirComb, dirCombStop);
% Print the results.
disp('Wave Direction')
disp('h-statistic:');
disp(t);
disp('p-value:');
disp(p);


if plot == 1
    figure('Name','Polar Histogram for wave direction in Motion Initiation and Termination');
    subplot(2,1,1);
    polarhistogram(dirComb,30);
    title('Wave Direction : Motion Initiation');
    subplot(2,1,2);
    polarhistogram(dirCombStop,30);
    title('Wave Direction : Motion Termination');

    figure('Name','Wave Direction in motion Initiation and termination');
    group = [ones(size(dirComb'));2.*ones(size(dirCombStop'))];
    boxplot([rad2deg(dirComb)';rad2deg(dirCombStop)'],group,'BoxStyle','filled','PlotStyle','compact');
    set(gca,'XTickLabel',{'Initiation','Termination'});
end

% Wave source points stats
sourceComb = horzcat(Waves.wavesStart(1:end).source);
sourceDen = zeros(rows,cols);
sourceCombStop = horzcat(Waves.wavesStop(1:end).source);
sourceDenStop = zeros(rows,cols);

for j=1:size(sourceComb,2)
    sourceDen(sourceComb(2,j),sourceComb(1,j)) = sourceDen(sourceComb(2,j),sourceComb(1,j)) + 1;
end
maxSourcePoint = max(sourceComb);

for j=1:size(sourceCombStop,2)
    sourceDenStop(sourceCombStop(2,j),sourceCombStop(1,j)) = sourceDenStop(sourceCombStop(2,j),sourceCombStop(1,j)) + 1;
end
maxSourcePointStop = max(sourceCombStop);

if plot == 1
    figure('Name','Spatial map of source points in Motion Initiation and Termination'); 
    subplot(2,1,1);
    imagesc(sourceDen);set(gca,'YDir','normal');
    title('Spatial map of sources points : Motion Inititaion'); colorbar;
    subplot(2,1,2);
    imagesc(sourceDenStop);set(gca,'YDir','normal');
    title('Spatial map of sources points : Motion Termination'); colorbar;
end 

wavesStat.evaluationPoints =  horzcat(Waves.wavesStart(1:end).evaluationPoints);
wavesStat.evaluationPointsStop =  horzcat(Waves.wavesStop(1:end).evaluationPoints);
wavesStat.speed = speedComb;
wavesStat.avgSpeed = avgSpeed;
wavesStat.speedStop = speedCombStop;
wavesStat.avgSpeedStop = avgSpeedStop;

wavesStat.velDir = dirComb;
wavesStat.avgDir = avgDir;
wavesStat.velDirStop = dirCombStop;
wavesStat.avgDirStop = avgDirStop;

wavesStat.sourcePoints = sourceComb;
wavesStat.maxSourcePoint = maxSourcePoint;
wavesStat.sourcePointsStop = sourceCombStop;
wavesStat.maxSourcePointStop = maxSourcePointStop;

wavesStat.WavesPerTrial = WavesPerTrial;
wavesStat.totalWaves = totalWaves;
wavesStat.WavesPerTrialStop = WavesPerTrialStop;
wavesStat.totalWavesStop = totalWavesStop;

end

