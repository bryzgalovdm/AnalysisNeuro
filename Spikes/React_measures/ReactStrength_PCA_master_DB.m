%%%%% ReactivationStrength_Peyrache2010_DB

%% Parameters
nmouse = [994];
Dir = PathForExperimentsERC_Dima('UMazePAG');
Dir = RestrictPathForExperiment(Dir,'nMice',nmouse);

% binsize for firing rate histogram (in timestamps!!!)
binsize = 2000; % =100ms (sampling rate = 20000Hz);

sav = 0;

%% Get Data

for j=1:length(Dir.path)
    cd(Dir.path{j}{1});
    load('SpikeData.mat','S','PlaceCells');
    load('behavResources.mat','SessionEpoch','CleanVtsd','CleanAlignedXtsd','CleanAlignedYtsd');
    load('Ripples.mat');
    if j == 7 || j==10 %%%% Mouse906 Mouse977
        load('SleepScoring_Accelero.mat','SWSEpoch','REMEpoch','Sleep');
    else
        load('SleepScoring_OBGamma.mat','SWSEpoch','REMEpoch','Sleep');
    end
end


%% Behavioral Epochs
% BaselineExplo Epoch
UMazeEpoch = or(SessionEpoch.Hab,SessionEpoch.TestPre1);
UMazeEpoch = or(UMazeEpoch,SessionEpoch.TestPre2);
UMazeEpoch = or(UMazeEpoch,SessionEpoch.TestPre3);
UMazeEpoch = or(UMazeEpoch,SessionEpoch.TestPre4);

CondEpoch = or(SessionEpoch.Cond1,SessionEpoch.Cond2);
CondEpoch = or(CondEpoch,SessionEpoch.Cond3);
CondEpoch = or(CondEpoch,SessionEpoch.Cond4);

% After Conditioning
AfterConditioningEpoch = or(SessionEpoch.TestPost1,SessionEpoch.TestPost2);
AfterConditioningEpoch = or(AfterConditioningEpoch,SessionEpoch.TestPost3);
AfterConditioningEpoch = or(AfterConditioningEpoch,SessionEpoch.TestPost4);

% Locomotion threshold
VtsdSmoothed  = tsd(Range(CleanVtsd),movmedian(Data(CleanVtsd),5));
LocomotionEpoch = thresholdIntervals(VtsdSmoothed,3,'Direction','Above');

% Get resulting epoch
UMazeMovingEpoch = and(LocomotionEpoch, UMazeEpoch);
AfterConditioningMovingEpoch = and(LocomotionEpoch, AfterConditioningEpoch);

%% Create template Epoch
% Ripples
ripSt = ts(ripples(:,1)*1e4);
PostRipplesEpoch = mergeCloseIntervals(intervalSet(Range(ripSt),Range(ripSt)+0.15*1e4),0.05*1e4);
PreRipplesEpoch = mergeCloseIntervals(intervalSet(Range(ripSt)-0.85*1e4,Range(ripSt)-0.7*1e4),0.05*1e4);

CondPostRipplesEpoch = and(CondEpoch,PostRipplesEpoch);
CondPreRipplesEpoch = and(CondEpoch,PreRipplesEpoch);

PSleepPostRipplesEpoch = and(SessionEpoch.PostSleep,PostRipplesEpoch);
PSleepPreRipplesEpoch = and(SessionEpoch.PostSleep,PreRipplesEpoch);



%% Create epochs to match
% PostSleep
PostSleepEpochSWS = and(SessionEpoch.PostSleep,SWSEpoch);
% CondEpoch
CondEpoch = CondEpoch;

%% Create an epoch to average

PostSleepRipplesEpochSWS = and(PostSleepEpochSWS,PostRipplesEpoch);

%% Make Q

Q = MakeQfromS(S,binsize);
% QTemplatePost = Restrict(Q,CondPostRipplesEpoch);
% QTemplatePre = Restrict(Q,CondPreRipplesEpoch);
QTemplatePost = Restrict(Q,PSleepPostRipplesEpoch);
QTemplatePre = Restrict(Q,PSleepPreRipplesEpoch);
DatTemplatePost = full(Data(QTemplatePost));
idx_zero_post = find(sum(DatTemplatePost)==0);
DatTemplatePre = full(Data(QTemplatePre));
idx_zero_pre = find(sum(DatTemplatePre)==0);
idx_zero = union(idx_zero_post,idx_zero_pre);
DatTemplatePost(:,idx_zero) = []; 
DatTemplatePre(:,idx_zero) = []; 



% QTemplate = Restrict(Q,UMazeMovingEpoch);
% DatTemplate = full(Data(QTemplate));
% idx_zero = find(sum(DatTemplate)==0);
% DatTemplate = DatTemplate(:,idx_zero);

QMatch = Restrict(Q,CondEpoch);
Dat_Match = full(Data(QMatch));
Dat_Match(:,idx_zero) = [];
% QMatch = Restrict(Q,SessionEpoch.PostSleep);

%% Do template

[templatesPost,correlations,eigenvalues,eigenvectors,lambdaMax] = ActivityTemplates_SB(DatTemplatePost,0);
[templatesPre,correlations,eigenvalues,eigenvectors,lambdaMax] = ActivityTemplates_SB(DatTemplatePre,0);

% [templates,correlations,eigenvalues,eigenvectors,lambdaMax] = ActivityTemplates_SB(Data(QTemplate),0);

%% Do matching
RStrengthPost = ReactivationStrength_SB(Dat_Match,templatesPost);
RStrengthPre = ReactivationStrength_SB(Dat_Match,templatesPre);

% RStrength = ReactivationStrength_SB(Dat_Match,templates);

%% See the result

for i=1:size(RStrengthPost,2)
    ToAvPost{i} = tsd(Range(QMatch),RStrengthPost(:,i));
end
for i=1:size(RStrengthPre,2)
    ToAvPre{i} = tsd(Range(QMatch),RStrengthPre(:,i));
end

%
% for i=1:size(RStrength,2)
%     ToAv{i} = tsd(Range(QMatch),RStrength(:,i));
% end

% % 
for i=1:size(RStrengthPost,2)
    [MPost{i},TPost{i}] = PlotRipRaw(ToAvPost{i},Start(CondPostRipplesEpoch,'s'),[-3000 3000]);
end
for i=1:size(RStrengthPre,2)
    [MPre{i},TPre{i}] = PlotRipRaw(ToAvPre{i},Start(CondPostRipplesEpoch,'s'),[-3000 3000]);
end
close all


% for i=1:size(RStrength,2)
%     [M{i},T{i}] = PlotRipRaw(ToAv{i},Start(CondPostRipplesEpoch,'s'),[-3000 3000]);
% end
% close all

%%%%%
for i=1:size(RStrengthPost,2)
    MPost_Av(1:length(MPost{1}),i) = MPost{i}(:,2);
end
for i=1:size(RStrengthPre,2)
    MPre_Av(1:length(MPre{1}),i) = MPre{i}(:,2);
end

% for i=1:size(RStrength,2)
%     M_Av(1:length(M{1}),i) = M{i}(:,2);
% end
% 
MPost_Av_Mean = mean(MPost_Av,2);
MPost_Av_std = std(MPost_Av,0,2);
MPre_Av_Mean = mean(MPre_Av,2);
MPre_Av_std = std(MPre_Av,0,2);

% M_Av_Mean = mean(M_Av,2);
% M_Av_std = std(M_Av,0,2);


%% Plot
% % Prepare x-axis
% for i=1:5:size(MPost{1},1)
%     XTicks{i} = num2str(round(MPost{1}(i,1),2));
% end
% XTicks = XTicks(~cellfun('isempty',XTicks));

figure

errorbar(MPost{1}(:,1),MPost_Av_Mean,MPost_Av_std);
hold on
errorbar(MPre{1}(:,1),MPre_Av_Mean,MPre_Av_std);
hold off
%     ylim([0 150])
set(gca,'XTick',[1:5:61],'XTickLabels',XTicks);

% errorbar(M{1}(:,1),M_Av_Mean,M_Av_std);
% %     ylim([0 150])
% set(gca,'XTick',[1:5:61],'XTickLabels',XTicks);
