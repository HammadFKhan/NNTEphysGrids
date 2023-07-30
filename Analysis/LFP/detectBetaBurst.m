function [BetaBurst] = detectBetaBurst(LFP,behaviourTrace,parameters)

% Reference - https://academic.oup.com/brain/article/140/11/2968/4430808#112576655

betaBurstThres = 75; % 75 percenticle 
xamp = abs(LFP);

%% Getting LFP for each trial and channels 
rows = size(LFP,1);
cols = size(LFP,2);
nlength = size(behaviourTrace(2).trace,2);
nTrials = size(behaviourTrace,2);

%xTrials = zeros(rows,cols,nlength);

for i=1:nTrials
    xTrials(i).xamp = squeeze(mean(xamp(:,:,behaviourTrace(i).LFPIndex(1):behaviourTrace(i).LFPIndex(end)),[1,2]));
end

%% Getting thresholds
betaThres = prctile(mean(xamp,[1,2]),betaBurstThres);

%% Getting beta burst segments 
for i=1:nTrials
    BetaBurst(i).burstsegments = findThresSeg(xTrials(i).xamp,betaThres);
end

