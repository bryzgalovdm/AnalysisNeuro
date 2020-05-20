%%%% NumDurEpisodesStages
% 
% A script - go to parameters
% 
% Calulates some properties of sleep in PAG
% Particularly,
% - overall length
% - number of episodes
% - duration of episodes
% 
% For REM and NREM sleep, and separately for NREM in the second hour of sleep
% 
% 20/05/2020 Dima Bryzgalov, MOBs team, Paris France
% github.com/bryzgalovdm

%% Parameters
% Mice that go in the analysis
nmouse = [797 798 828 861 882 905 906 911 912 977 994];
% nmouse = [906 912]; % Had PreMazes
% nmouse = [905 911]; % Did not have PreMazes

% Get paths
Dir = PathForExperimentsERC_Dima('UMazePAG');
% Dir = PathForExperimentsERC_DimaMAC('UMazePAG');
Dir = RestrictPathForExperiment(Dir,'nMice',nmouse);

% Sleep time to restrict
SleepTimeToRestrict = 2*60*60*1e4; % 2 hours

% Do you want to save the figures?
savfig = true;

% Paths and names to save
pathfig = '/MOBS_workingON/Dima/Ongoing_results/Sleep/NumDur/'; % Without dropbox

%% Preallocation
% Data
behav = cell(length(Dir.path), 1);
sleep = cell(length(Dir.path), 1);
% Epochs
NREMPre = cell(length(Dir.path), 1);
REMPre = cell(length(Dir.path), 1);
NREMPost = cell(length(Dir.path), 1);
REMPost = cell(length(Dir.path), 1);
NREMPre_lastH = cell(length(Dir.path), 1);
NREMPost_lastH = cell(length(Dir.path), 1);
% Results - overall length
lenNREMPre = zeros(length(Dir.path), 1);
lenREMPre = zeros(length(Dir.path), 1);
lenNREMPost = zeros(length(Dir.path), 1);
lenREMPost = zeros(length(Dir.path), 1);
% Results - numEpisodes
NumNREMPre = zeros(length(Dir.path), 1);
NumREMPre = zeros(length(Dir.path), 1);
NumNREMPost = zeros(length(Dir.path), 1);
NumREMPost = zeros(length(Dir.path), 1);
% Results - duration Episodes
DurNREMPre = cell(length(Dir.path), 1);
DurREMPre = cell(length(Dir.path), 1);
DurNREMPost = cell(length(Dir.path), 1);
DurREMPost = cell(length(Dir.path), 1);
% Results - NREM #spikes last hour
lenNREMPre_lastH = nan(length(Dir.path), 1);
lenNREMPost_lastH = nan(length(Dir.path), 1);
numNREMPre_lastH = nan(length(Dir.path), 1);
numNREMPost_lastH = nan(length(Dir.path), 1);
DurNREMPre_lastH = cell(length(Dir.path), 1);
DurNREMPost_lastH = cell(length(Dir.path), 1);

%% Load data - here I handle the exceptions too

for i=1:length(Dir.path)
    behav{i} = load([Dir.path{i}{1} 'behavResources.mat'],'SessionEpoch');
    if strcmp(Dir.name{i}, 'Mouse906') || strcmp(Dir.name{i}, 'Mouse977') % Mice with bad OB-based sleep scoring
        sleep{i} = load([Dir.path{i}{1} 'SleepScoring_Accelero.mat'],'SWSEpoch','REMEpoch','Sleep'); % Sleep is not used
    else
        sleep{i} = load([Dir.path{i}{1} 'SleepScoring_OBGamma.mat'],'SWSEpoch','REMEpoch','Sleep');  % Sleep is not used
    end
end

%% Create epochs

for i=1:length(Dir.path)
    try % Restrict to the first two hours if applicable
        NREMPre{i} = and(RestrictToTime(behav{i}.SessionEpoch.PreSleep, SleepTimeToRestrict), sleep{i}.SWSEpoch);
        REMPre{i} = and(RestrictToTime(behav{i}.SessionEpoch.PreSleep, SleepTimeToRestrict), sleep{i}.REMEpoch);
    catch
        NREMPre{i} = and(behav{i}.SessionEpoch.PreSleep, sleep{i}.SWSEpoch);
        REMPre{i} = and(behav{i}.SessionEpoch.PreSleep, sleep{i}.REMEpoch);
    end
    try % Restrict to the first two hours if applicable
        NREMPost{i} = and(RestrictToTime(behav{i}.SessionEpoch.PostSleep, SleepTimeToRestrict), sleep{i}.SWSEpoch);
        REMPost{i} = and(RestrictToTime(behav{i}.SessionEpoch.PostSleep, SleepTimeToRestrict), sleep{i}.REMEpoch);
    catch
        NREMPost{i} = and(behav{i}.SessionEpoch.PostSleep, sleep{i}.SWSEpoch);
        REMPost{i} = and(behav{i}.SessionEpoch.PostSleep, sleep{i}.REMEpoch);
    end
end

% Create the last hour of sleep - NREM (exclude the mice with less than 2 hours of sleep session)
for i = 1:length(Dir.path)
    % Pre
    temp = SplitIntervals(behav{i}.SessionEpoch.PreSleep, SleepTimeToRestrict/2);
    if length(temp) > 1
        NREMPre_lastH{i} = and(temp{2}, sleep{i}.SWSEpoch); % Second hour
    end
    clear temp
    % Post
    temp = SplitIntervals(behav{i}.SessionEpoch.PostSleep, SleepTimeToRestrict/2);
    if length(temp) > 1
        NREMPost_lastH{i} = and(temp{2}, sleep{i}.SWSEpoch); % Second hour
    end
    clear temp
end

%% Calculate overall length of sleep stages
for i=1:length(Dir.path)
    lenNREMPre(i) = sum(End(NREMPre{i},'s') - Start(NREMPre{i}, 's'));
    lenREMPre(i) = sum(End(REMPre{i},'s') - Start(REMPre{i}, 's'));
    lenNREMPost(i) = sum(End(NREMPost{i},'s') - Start(NREMPost{i}, 's'));
    lenREMPost(i) = sum(End(REMPost{i},'s') - Start(REMPost{i}, 's'));
end

%% Calculate number of sleep stage episodes
for i = 1:length(Dir.path)
    NumNREMPre(i) = length(Start(NREMPre{i}));
    NumREMPre(i) = length(Start(REMPre{i}));
    NumNREMPost(i) = length(Start(NREMPost{i}));
    NumREMPost(i) = length(Start(REMPost{i}));
end

%% Calculate duration of episodes

for i = 1:length(Dir.path)
    % PreAllocate matrices
    DurNREMPre{i} = zeros(length(Start(NREMPre{i})), 1);
    DurREMPre{i} = zeros(length(Start(REMPre{i})), 1);
    DurNREMPost{i} = zeros(length(Start(NREMPost{i})), 1);
    DurREMPost{i} = zeros(length(Start(NREMPost{i})), 1);
    
    % PreNREM
    st = Start(NREMPre{i},'s');
    fin = End(NREMPre{i},'s');
    DurNREMPre{i} = arrayfun(@(x) fin(x)-st(x), [1:length(Start(NREMPre{i}))]);
    clear st fin
    
    % PreREM
    st = Start(REMPre{i},'s');
    fin = End(REMPre{i},'s');
    DurREMPre{i} = arrayfun(@(x) fin(x)-st(x), [1:length(Start(REMPre{i}))]);
    clear st fin
    
    % PostNREM
    st = Start(NREMPost{i},'s');
    fin = End(NREMPost{i},'s');
    DurNREMPost{i} = arrayfun(@(x) fin(x)-st(x), [1:length(Start(NREMPost{i}))]);
    clear st fin
    
    % PostNREM
    st = Start(REMPost{i},'s');
    fin = End(REMPost{i},'s');
    DurREMPost{i} = arrayfun(@(x) fin(x)-st(x), [1:length(Start(REMPost{i}))]);
    clear st fin
end

%% Calculate number of episodes in the last hour of NREM

for i = 1:length(Dir.path)
    
    if ~isempty(NREMPre_lastH{i})
        numNREMPre_lastH(i) = length(Start(NREMPre_lastH{i}));
    end
    if ~isempty(NREMPost_lastH{i})
        numNREMPost_lastH(i) = length(Start(NREMPost_lastH{i}));
    end
    
    if ~isempty(NREMPre_lastH{i})
        lenNREMPre_lastH(i) = sum(End(NREMPre_lastH{i},'s') - Start(NREMPre_lastH{i}, 's'));
    end
    if ~isempty(NREMPost_lastH{i})
        lenNREMPost_lastH(i) = sum(End(NREMPost_lastH{i},'s') - Start(NREMPost_lastH{i}, 's'));
    end
        
    
    if ~isempty(NREMPre_lastH{i})
        % Pre-Allocate
        DurNREMPre_lastH{i} = zeros(length(Start(NREMPre_lastH{i})), 1);
        % PreNREM
        st = Start(NREMPre_lastH{i},'s');
        fin = End(NREMPre_lastH{i},'s');
        DurNREMPre_lastH{i} = arrayfun(@(x) fin(x)-st(x), [1:length(Start(NREMPre_lastH{i}))]);
        clear st fin
    else
        DurNREMPre_lastH{i} = [];
    end
    if ~isempty(NREMPost_lastH{i})
        % PreAllocate
        DurNREMPost_lastH{i} = zeros(length(Start(NREMPost_lastH{i})), 1);
        % PostNREM
        st = Start(NREMPost_lastH{i},'s');
        fin = End(NREMPost_lastH{i},'s');
        DurNREMPost_lastH{i} = arrayfun(@(x) fin(x)-st(x), [1:length(Start(NREMPost_lastH{i}))]);
        clear st fin
    else
        DurNREMPost_lastH{i} = [];
    end
    
    
end

%% Pool durations

% FR of neurons
Dur.Pre.NREM = DurNREMPre{1}';
Dur.Pre.REM = DurREMPre{1}';
Dur.Post.NREM = DurNREMPost{1}';
Dur.Post.REM = DurREMPost{1}';

Dur_last.Pre.NREM = DurNREMPre_lastH{1}';
Dur_last.Post.NREM = DurNREMPost_lastH{1}';

if length(Dir.path) > 1
    for i = 2:length(Dir.path)
        Dur.Pre.NREM = [Dur.Pre.NREM; DurNREMPre{i}'];
        Dur.Pre.REM = [Dur.Pre.REM; DurREMPre{i}'];
        Dur.Post.NREM = [Dur.Post.NREM; DurNREMPost{i}'];
        Dur.Post.REM = [Dur.Post.REM; DurREMPost{i}'];
        
        Dur_last.Pre.NREM = [Dur_last.Pre.NREM; DurNREMPre_lastH{i}'];
        Dur_last.Post.NREM = [Dur_last.Post.NREM; DurNREMPost_lastH{i}'];
    end
end

%% Illustrate the data

%%%%%%%%%%% WHOLE SLEEP SESSION %%%%%%%%%%%%%%%%%%%%%%%%%
% NREM Sleep
Pl = {NumNREMPre, NumNREMPost; lenNREMPre, lenNREMPost; Dur.Pre.NREM, Dur.Post.NREM};
Cols = {[0.7 0.7 0.7], [0.2 0.2 0.2]};
xlabs = {'PreNREM','PostNREM'};
ylabs = {'# episodes', 'Duration (s)', 'Mean duration (s)'};
tits = {'Number of episodes', 'Overall length', 'Duration of episodes'};
f1 = figure('units', 'normalized', 'outerposition', [0 0 0.9 0.55]);
ax = arrayfun(@(i) subplot(1,3,i, 'NextPlot', 'add', 'Box', 'off'), [1:3]);
for i=1:2
    axes(ax(i));
    MakeBoxPlot_DB(Pl(i,:),Cols,1:2,[],1);
    set(gca,'XTick', [1:2], 'XTickLabel', xlabs,...
        'FontSize', 16, 'FontWeight', 'bold');
    ylabel(ylabs{i});
    title(tits{i});
end
axes(ax(3));
MakeBoxPlot_DB(Pl(3,:),Cols,1:2,[],0);
set(gca,'XTick', [1:2], 'XTickLabel', xlabs,...
        'FontSize', 16, 'FontWeight', 'bold');
ylabel(ylabs{3})
title(tits{3});
if savfig
    saveas(f1,[dropbox pathfig 'NREMSleepChar.fig']);
    saveFigure(f1,'NREMSleepChar',[dropbox pathfig]);
end

% REM Sleep
Pl = {NumREMPre, NumREMPost; lenREMPre, lenREMPost; Dur.Pre.REM, Dur.Post.REM};
Cols = {[0.7 0.7 0.7], [0.2 0.2 0.2]};
xlabs = {'PreREM','PostREM'};
ylabs = {'# episodes', 'Duration (s)', 'Mean duration (s)'};
tits = {'Number of episodes', 'Overall length', 'Duration of episodes'};
f2 = figure('units', 'normalized', 'outerposition', [0 0 0.9 0.55]);
ax = arrayfun(@(i) subplot(1,3,i, 'NextPlot', 'add', 'Box', 'off'), [1:3]);
for i=1:2
    axes(ax(i));
    MakeBoxPlot_DB(Pl(i,:),Cols,1:2,[],1);
    set(gca,'XTick', [1:2], 'XTickLabel', xlabs,...
        'FontSize', 16, 'FontWeight', 'bold');
    ylabel(ylabs{i});
    title(tits{i});
end
axes(ax(3));
MakeBoxPlot_DB(Pl(3,:),Cols,1:2,[],0);
set(gca,'XTick', [1:2], 'XTickLabel', xlabs,...
        'FontSize', 16, 'FontWeight', 'bold');
ylabel(ylabs{3})
title(tits{3});
if savfig
    saveas(f2,[dropbox pathfig 'REMSleepChar.fig']);
    saveFigure(f2,'REMSleepChar',[dropbox pathfig]);
end
%%%%%%%%%%% WHOLE SLEEP SESSION %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% SECOND HOUR HOUR %%%%%%%%%%%%%%%%%%%%%%%%%
Pl = {numNREMPre_lastH, numNREMPost_lastH; lenNREMPre_lastH, lenNREMPost_lastH; Dur_last.Pre.NREM, Dur_last.Post.NREM};
Cols = {[0.7 0.7 0.7], [0.2 0.2 0.2]};
xlabs = {'PreNREM','PostNREM'};
tits = {'Number of episodes', 'Overall length', 'Duration of episodes'};
f3 = figure('units', 'normalized', 'outerposition', [0 0 0.9 0.55]);
ax = arrayfun(@(i) subplot(1,3,i, 'NextPlot', 'add', 'Box', 'off'), [1:3]);
for i=1:2
    axes(ax(i));
    MakeBoxPlot_DB(Pl(i,:),Cols,1:2,[],1);
    set(gca,'XTick', [1:2], 'XTickLabel', xlabs,...
        'FontSize', 16, 'FontWeight', 'bold');
    ylabel(ylabs{i});
    title(tits{i});
end
axes(ax(3));
MakeBoxPlot_DB(Pl(3,:),Cols,1:2,[],0);
set(gca,'XTick', [1:2], 'XTickLabel', xlabs,...
        'FontSize', 16, 'FontWeight', 'bold');
ylabel(ylabs{3})
title(tits{3});
if savfig
    saveas(f3,[dropbox pathfig 'NREM_lasthour_SleepChar.fig']);
    saveFigure(f3,'NREM_lasthour_SleepChar',[dropbox pathfig]);
end
%%%%%%%%%%% SECOND HOUR HOUR %%%%%%%%%%%%%%%%%%%%%%%%%