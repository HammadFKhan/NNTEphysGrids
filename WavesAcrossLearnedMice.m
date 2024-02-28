%% Waves across Naive mice
NaiveWaves = Waves.wavesMiss;
ExpertWaves = Waves1.wavesMiss;
%%
ExpertPGD = arrayfun(@(x) x.PGD,ExpertWaves,'UniformOutput',false); ExpertPGD = vertcat(ExpertPGD{:});
NaivePGD = arrayfun(@(x) x.PGD,NaiveWaves,'UniformOutput',false); NaivePGD = vertcat(NaivePGD{:});

ExpertWaveDen = arrayfun(@(x) x.nWaves,ExpertWaves)';
NaiveWaveDen = arrayfun(@(x) x.nWaves,NaiveWaves)';
%%
figure,histogram(ExpertWaveDen,1:1:15,'Normalization','probability')
hold on, histogram(NaiveWaveDen,1:1:15,'Normalization','probability')
%%
baseLinePGDExpert = smoothdata(std(ExpertPGD,[],1),'movmean',100);
baseLinePGDNaive = smoothdata(std(NaivePGD,[],1),'movmean',500);
%%
figure,plot(baseLinePGDExpert),hold on
plot(baseLinePGDNaive);box off,set(gca,'Tickdir','out','FontSize',16)
%% Compare by cue driven changes
ExpertPGDchange = mean(ExpertPGD(:,1501:3001),2)-mean(ExpertPGD(:,1:1500),2);
NaivePGDchange = mean(NaivePGD(:,1501:3001),2)-mean(NaivePGD(:,1:1500),2);