function [MI] = getMI(IntanBehaviour,z_score,nPerm,plotFlag,parameters)

% Reference : Propagating waves mediate information transfer in the motor
% cortex, Nature Neuro, 2006 

nElectrodes = parameters.rows*parameters.cols;

% For phase - Cue encoding
xgpHit = arrayfun(@(s) angle(s.xgp), IntanBehaviour.cueHitTrace, 'UniformOutput', false);
xgpMiss = arrayfun(@(s) angle(s.xgp), IntanBehaviour.cueMissTrace, 'UniformOutput', false);
[MI.PhaseCue] = getMutualInformation(xgpHit,xgpMiss,parameters);

if plotFlag == 1
    figure();
    title("Mututal Information across all electrodes - Phase")
    a = peakSort2DArray(reshape(MI.PhaseCue,[],size(MI.PhaseCue,3)),'descend',2);
    a = removeNaNRows(a);
    imagesc(IntanBehaviour.cueHitTrace(1).time,1:size(a,1),a); colormap(hot);
    ylabel("Electrodes");xlabel("Time (s)"); 
    h = colorbar; h.Label.String = 'Information (bits)';
    xline(0,'-w','Cue','LabelVerticalAlignment','top');
    xline(mean(IntanBehaviour.reactionTime,'all'),'-w','Avg. Reaction Time','LabelVerticalAlignment','top');
end

% For phase - Reward encoding
xgpHitReward = arrayfun(@(s) angle(s.xgp), IntanBehaviour.hitTrace, 'UniformOutput', false);
xgpFA = arrayfun(@(s) angle(s.xgp), IntanBehaviour.missTrace, 'UniformOutput', false);
[MI.PhaseReward] = getMutualInformation(xgpHitReward,xgpFA,parameters);

if plotFlag == 1
    figure();
    title("Mututal Information (Reward) across all electrodes - Phase")
    a = peakSort2DArray(reshape(MI.PhaseReward,[],size(MI.PhaseReward,3)),'descend',2);
    a = removeNaNRows(a);
    imagesc(IntanBehaviour.hitTrace(1).time,1:size(a,1),a); colormap(hot);
    ylabel("Electrodes");xlabel("Time (s)"); 
    h = colorbar; h.Label.String = 'Information (bits)';
    xline(0,'-w','Reward','LabelVerticalAlignment','top');
    xline(-1*mean(IntanBehaviour.reactionTime,'all'),'-w','Avg. Cue Time','LabelVerticalAlignment','top');
end

% For phase - Motion encoding
xgpMIHit = arrayfun(@(s) angle(s.xgp), IntanBehaviour.MIHitTrace, 'UniformOutput', false);
xgpMIFA = arrayfun(@(s) angle(s.xgp), IntanBehaviour.MIFATrace, 'UniformOutput', false);
[MI.PhaseMI] = getMutualInformation(xgpMIHit,xgpMIFA,parameters);

if plotFlag == 1
    figure();
    title("Mututal Information (Motion Initiation) across all electrodes - Phase")
    a = peakSort2DArray(reshape(MI.PhaseMI,[],size(MI.PhaseMI,3)),'descend',2);
    a = removeNaNRows(a);
    imagesc(IntanBehaviour.MIHitTrace(1).time,1:size(a,1),a); colormap(sky);
    ylabel("Electrodes");xlabel("Time (s)"); 
    h = colorbar; h.Label.String = 'Information (bits)';
    xline(0,'-w','MI','LabelVerticalAlignment','top');
end

if z_score == 1
    xgpComb = [xgpHit xgpMiss];
    nHit = size(xgpHit,2);
    nMiss = size(xgpMiss,2);
    nTot = nHit + nMiss;
    
    nullDistCue = zeros(parameters.rows,parameters.cols,size(MI.PhaseCue,3),nPerm);
    for j=1:nPerm
        randIndex = randperm(nTot);
        xgpHitRand = xgpComb(randIndex(1:nHit));
        xgpMissRand = xgpComb(randIndex(nHit+1:end));
        nullDistCue(:,:,:,j) = getMutualInformation(xgpHitRand,xgpMissRand,parameters);
        j
    end
    MI.muCue = mean(nullDistCue,4); % Mean of the null distribution
    MI.sigmaCue = std(nullDistCue,0,4); % Standard deviation of null distribution
    
    MI.PhaseCuez = (MI.PhaseCue-MI.muCue)./MI.sigmaCue;

    if plotFlag == 1
        figure();
        title("Mututal Information in Phase for Cue - z-scored")
        imagesc(IntanBehaviour.cueHitTrace(1).time,1:32,peakSort2DArray(reshape(MI.PhaseCuez,[],size(MI.PhaseCuez,3)),'descend',2)); colormap(hot);
        ylabel("Electrodes");xlabel("Time (s)"); 
        h = colorbar; h.Label.String = 'Information (z-scored)';
        xline(0,'-w','Cue','LabelVerticalAlignment','top');
        xline(mean(IntanBehaviour.reactionTime,'all'),'-w','Avg. Reaction Time','LabelVerticalAlignment','top');
    end

    xgpComb = [xgpHitReward xgpFA];
    nHit = size(xgpHitReward,2);
    nMiss = size(xgpFA,2);
    nTot = nHit + nMiss;
    
    nullDistReward = zeros(parameters.rows,parameters.cols,size(MI.PhaseReward,3),nPerm);
    for j=1:nPerm
        randIndex = randperm(nTot);
        xgpHitRand = xgpComb(randIndex(1:nHit));
        xgpMissRand = xgpComb(randIndex(nHit+1:end));
        nullDistReward(:,:,:,j) = getMutualInformation(xgpHitRand,xgpMissRand,parameters);
        j
    end
    MI.muReward = mean(nullDistReward,4); % Mean of the null distribution
    MI.sigmaReward = std(nullDistReward,0,4); % Standard deviation of null distribution
    
    MI.PhaseRewardz = (MI.PhaseReward-MI.muReward)./MI.sigmaReward;

    if plotFlag == 1
        figure();
        title("Mututal Information in Phase for Reward - z-scored")
        imagesc(IntanBehaviour.cueHitTrace(1).time,1:32,peakSort2DArray(reshape(MI.PhaseRewardz,[],size(MI.PhaseRewardz,3)),'descend',2)); colormap(hot);
        ylabel("Electrodes");xlabel("Time (s)"); 
        h = colorbar; h.Label.String = 'Information (z-scored)';
        xline(0,'-w','Reward','LabelVerticalAlignment','top');
        xline(-mean(IntanBehaviour.reactionTime,'all'),'-w','Avg. Cue Time','LabelVerticalAlignment','top');
    end

end

