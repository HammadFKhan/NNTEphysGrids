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
[t, p] = ttest2(speedComb, speedCombStop);
% Print the results.
disp('Wave Speed')
disp('t-statistic:');
disp(t);
disp('p-value:');
disp(p);


if plot == 1
    figure('Name','Histogram of wave speeds in motion Initiation and termination');
    subplot(2,1,1);
    histogram(speedComb,100);
    xline(avgSpeed,'-r',{'Mean speed = ' num2str(avgSpeed) ' cm/s'});
    xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Motion Initiation');
    subplot(2,1,2);
    histogram(speedCombStop,100);
    xline(avgSpeedStop,'-r',{'Mean speed = ' num2str(avgSpeedStop) ' cm/s'});
    xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Motion Termination');

    figure('Name','Wave speeds in motion Initiation and termination');
    group = [ones(size(speedComb')); 2.*ones(size(speedCombStop'))];
    boxplot([speedComb';speedCombStop'],group);
    set(gca,'XTickLabel',{'Initiation','Termination'});
end

% Wave direction stats
dirComb = horzcat(Waves.wavesStart(1:end).waveDir);
avgDir = mean(dirComb);
dirCombStop = horzcat(Waves.wavesStop(1:end).waveDir);
avgDirStop = mean(dirCombStop);

[t, p] = ttest2(dirComb, dirCombStop);
% Print the results.
disp('Wave Direction')
disp('t-statistic:');
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
    boxplot([rad2deg(dirComb)';rad2deg(dirCombStop)'],group);
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

