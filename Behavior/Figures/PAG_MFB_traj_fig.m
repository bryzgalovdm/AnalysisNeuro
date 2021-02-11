function fh = PAG_MFB_traj_fig(numtest, time_restrict, varargin)
%
% It plots trajectories of mice from reward-learning exeriment and aversive
% learning experiment
%
%  INPUT
%
%       numtest         first <numtest> of tests will be taken to show
%                       trajectories
%       time_restrict   time to restrict each trajectory in time ([] - don't restrict)
%       IsSave          true if you want to save the resulting figure in
%                           Dropbox (boolean - default = false)
%       MatchTraj       true if you want to match number of trajectories (default - false)
%
%
%  OUTPUT
%
%       Figure
%
%       See
%
%       BehaviorERC, PathForExperimentERC_Dima
%
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 04/02/2021
% github.com/bryzgalovdm

%% Parameters management
if isempty(time_restrict)
    IsRestrict = false;
else
    IsRestrict = true;
end

IsSave = false;
IsMatch = false;
% Optional parameters
for i=1:2:length(varargin)
    switch(lower(varargin{i}))
        case 'issave'
            IsSave = varargin{i+1};
            if IsSave ~= 1 && IsSave ~= 0
                error('Incorrect value for property ''IsSave'' (type ''help ExampleTrajectory_DB'' for details).');
            end
        case 'matchedtraj'
            IsMatch = varargin{i+1};
            if IsMatch ~= 1 && IsMatch ~= 0
                error('Incorrect value for property ''MatchedTraj'' (type ''help ExampleTrajectory_DB'' for details).');
            end
    end
end

%% Get directories
DirMFB = PathForExperimentsERC_SL('StimMFBWake');
DirMFB = RestrictPathForExperiment(DirMFB,'nMice', [882 941 117 161 162 168]);
DirPAG = PathForExperimentsERC_Dima('UMazePAG');
DirPAG = RestrictPathForExperiment(DirPAG,'nMice', [797 798 828 861 882 905 906 911 912 977 994 1117 1161 1162 1168]);

%% Get data
dat_MFB = cell(length(DirMFB.path),1);
dat_PAG = cell(length(DirPAG.path),1);
for idir = 1:length(DirMFB.path)
    dat_MFB{idir} = load([DirMFB.path{idir}{1} 'behavResources.mat'], 'behavResources', 'SessionEpoch');
end
for idir = 1:length(DirPAG.path)
    dat_PAG{idir} = load([DirPAG.path{idir}{1} 'behavResources.mat'], 'behavResources', 'SessionEpoch');
end

%% Find indices of pre and post tests
id_Pre_MFB = cell(length(DirMFB.path),1);
id_Cond_MFB = cell(length(DirMFB.path),1);
id_Post_MFB = cell(length(DirMFB.path),1);
for idir = 1:length(DirMFB.path)
    id_Pre_MFB{idir} = FindSessionID_ERC(dat_MFB{idir}.behavResources,'TestPre');
    id_Cond_MFB{idir} = FindSessionID_ERC(dat_MFB{idir}.behavResources,'Cond');
    id_Post_MFB{idir} = FindSessionID_ERC(dat_MFB{idir}.behavResources,'TestPost');
end
id_Pre_PAG = cell(length(DirPAG.path),1);
id_Cond_PAG = cell(length(DirPAG.path),1);
id_Post_PAG = cell(length(DirPAG.path),1);
for idir = 1:length(DirPAG.path)
    id_Pre_PAG{idir} = FindSessionID_ERC(dat_PAG{idir}.behavResources,'TestPre');
    id_Cond_PAG{idir} = FindSessionID_ERC(dat_PAG{idir}.behavResources,'Cond');
    id_Post_PAG{idir} = FindSessionID_ERC(dat_PAG{idir}.behavResources,'TestPost');
end


%% Cut the data to 2 min and organize the data into cells
% MFB data
cnt_pre=1;
cnt_cond=1;
cnt_post = 1;
for idir = 1:length(DirMFB.path)
    for itest = 1:numtest
        if IsRestrict
            time = Range(dat_MFB{idir}.behavResources(id_Pre_MFB{idir}(itest)).AlignedXtsd);
            Epoch = intervalSet(time(1), time(1)+time_restrict*1e4);
            traj_pre_mfb{cnt_pre} = [Data(Restrict(dat_MFB{idir}.behavResources(id_Pre_MFB{idir}(itest)).AlignedXtsd, Epoch)) ...
                Data(Restrict(dat_MFB{idir}.behavResources(id_Pre_MFB{idir}(itest)).AlignedYtsd, Epoch))];
        else
            traj_pre_mfb{cnt_pre} = [Data(dat_MFB{idir}.behavResources(id_Pre_MFB{idir}(itest)).AlignedXtsd) ...
                Data(dat_MFB{idir}.behavResources(id_Pre_MFB{idir}(itest)).AlignedYtsd)];
        end
        cnt_pre=cnt_pre+1;
    end
    for itest = 1:numtest
        traj_cond_mfb{cnt_cond} = [Data(dat_MFB{idir}.behavResources(id_Cond_MFB{idir}(itest)).AlignedXtsd) ...
            Data(dat_MFB{idir}.behavResources(id_Cond_MFB{idir}(itest)).AlignedYtsd)];
        cnt_cond=cnt_cond+1;
    end
    for itest = 1:numtest
        if IsRestrict
            time = Range(dat_MFB{idir}.behavResources(id_Post_MFB{idir}(itest)).AlignedXtsd);
            Epoch = intervalSet(time(1), time(1)+time_restrict*1e4);
            traj_post_mfb{cnt_post} = [Data(Restrict(dat_MFB{idir}.behavResources(id_Post_MFB{idir}(itest)).AlignedXtsd, Epoch)) ...
                Data(Restrict(dat_MFB{idir}.behavResources(id_Post_MFB{idir}(itest)).AlignedYtsd, Epoch))];
        else
            traj_post_mfb{cnt_post} = [Data(dat_MFB{idir}.behavResources(id_Post_MFB{idir}(itest)).AlignedXtsd) ...
                Data(dat_MFB{idir}.behavResources(id_Post_MFB{idir}(itest)).AlignedYtsd)];
        end
        cnt_post=cnt_post+1;
    end
    
end

% PAG data
cnt_pre=1;
cnt_cond=1;
cnt_post=1;
for idir = 1:length(DirPAG.path)
    
    %% Store in the variable
    for itest = 1:numtest
        if IsRestrict
            time = Range(dat_PAG{idir}.behavResources(id_Pre_PAG{idir}(itest)).CleanAlignedXtsd);
            Epoch = intervalSet(time(1), time(1)+time_restrict*1e4);
            traj_pre_pag{cnt_pre} = [Data(Restrict(dat_PAG{idir}.behavResources(id_Pre_PAG{idir}(itest)).CleanAlignedXtsd, Epoch)) ...
                Data(Restrict(dat_PAG{idir}.behavResources(id_Pre_PAG{idir}(itest)).CleanAlignedYtsd, Epoch))];
        else
            traj_pre_pag{cnt_pre} = [Data(dat_PAG{idir}.behavResources(id_Pre_PAG{idir}(itest)).CleanAlignedXtsd) ...
                Data(dat_PAG{idir}.behavResources(id_Pre_PAG{idir}(itest)).CleanAlignedYtsd)];
        end
        cnt_pre=cnt_pre+1;
    end
    for itest = 1:numtest
        traj_cond_pag{cnt_cond} = [Data(dat_PAG{idir}.behavResources(id_Cond_PAG{idir}(itest)).CleanAlignedXtsd) ...
            Data(dat_PAG{idir}.behavResources(id_Cond_PAG{idir}(itest)).CleanAlignedYtsd)];
        cnt_cond=cnt_cond+1;
    end
    for itest = 1:numtest
        if IsRestrict
            time = Range(dat_PAG{idir}.behavResources(id_Post_PAG{idir}(itest)).CleanAlignedXtsd);
            Epoch = intervalSet(time(1), time(1)+time_restrict*1e4);
            traj_post_pag{cnt_post} = [Data(Restrict(dat_PAG{idir}.behavResources(id_Post_PAG{idir}(itest)).CleanAlignedXtsd, Epoch)) ...
                Data(Restrict(dat_PAG{idir}.behavResources(id_Post_PAG{idir}(itest)).CleanAlignedYtsd, Epoch))];
        else
            traj_post_pag{cnt_post} = [Data(dat_PAG{idir}.behavResources(id_Post_PAG{idir}(itest)).CleanAlignedXtsd) ...
                Data(dat_PAG{idir}.behavResources(id_Post_PAG{idir}(itest)).CleanAlignedYtsd)];
        end
        cnt_post=cnt_post+1;
    end
    
end

%% Set up the figure

% Axes
fh = figure('units', 'normalized', 'outerposition', [0 0 0.9 0.9]);
PreTestTrajPAG_Axes = axes('position', [0.07 0.07 0.28 0.38]);
CondTrajPAG_Axes = axes('position', [0.38 0.07 0.28 0.38]);
PostTestTrajPAG_Axes = axes('position', [0.69 0.07 0.28 0.38]);

PreTestTrajMFB_Axes = axes('position', [0.07 0.54 0.28 0.38]);
CondTrajMFB_Axes = axes('position', [0.38 0.54 0.28 0.38]);
PostTestTrajMFB_Axes = axes('position', [0.69 0.54 0.28 0.38]);

maze = [0 0; 0 1; 1 1; 1 0; 0.63 0; 0.63 0.75; 0.35 0.75; 0.35 0; 0 0];
shockZone = [0 0; 0 0.43; 0.35 0.43; 0.35 0; 0 0];

% Plot
% Trajectories
axes(PreTestTrajMFB_Axes);
hold on
for i=1:length(traj_pre_mfb)
    plot(traj_pre_mfb{i}(:,1), traj_pre_mfb{i}(:,2),'LineWidth',2, 'Color', [0.8 0.8 0.8]);
end
box off
set(gca,'XtickLabel',{},'YtickLabel',{}, 'XDir', 'reverse');
plot(maze(:,1),maze(:,2),'k','LineWidth',3)
plot(shockZone(:,1),shockZone(:,2),'g','LineWidth',3)
title('PreTests MFB','FontSize',18,'FontWeight','bold');
xlim([0 1])
ylim([0 1])
hold off

axes(CondTrajMFB_Axes);
hold on
for i=1:length(traj_cond_mfb)
    plot(traj_cond_mfb{i}(:,1), traj_cond_mfb{i}(:,2),'LineWidth',2, 'Color', [0.8 0.8 0.8]);
end
box off
set(gca,'XtickLabel',{},'YtickLabel',{}, 'XDir', 'reverse');
plot(maze(:,1),maze(:,2),'k','LineWidth',3)
plot(shockZone(:,1),shockZone(:,2),'g','LineWidth',3)
title('Conditioning MFB','FontSize',18,'FontWeight','bold');
xlim([0 1])
ylim([0 1])
hold off

% PostTests
axes(PostTestTrajMFB_Axes);
hold on
for i=1:length(traj_post_mfb)
    plot(traj_post_mfb{i}(:,1), traj_post_mfb{i}(:,2),'LineWidth',2, 'Color', [0.8 0.8 0.8]);
end
box off
set(gca,'XtickLabel',{},'YtickLabel',{}, 'XDir', 'reverse');
plot(maze(:,1),maze(:,2),'k','LineWidth',3)
plot(shockZone(:,1),shockZone(:,2),'g','LineWidth',3)
title('PostTests MFB','FontSize',18,'FontWeight','bold');
xlim([0 1])
ylim([0 1])
hold off

% get indices
if IsMatch
    idx = datasample(1:length(traj_pre_pag), length(traj_pre_mfb), 'Replace', false);
end

axes(PreTestTrajPAG_Axes);
hold on
if IsMatch
    for i=1:idx
        plot(traj_pre_pag{i}(:,1), traj_pre_pag{i}(:,2),'LineWidth',2, 'Color', [0.8 0.8 0.8]);
    end
else
    for i=1:length(traj_pre_pag)
        plot(traj_pre_pag{i}(:,1), traj_pre_pag{i}(:,2),'LineWidth',2, 'Color', [0.8 0.8 0.8]);
    end
end
box off
set(gca,'XtickLabel',{},'YtickLabel',{});
plot(maze(:,1),maze(:,2),'k','LineWidth',3)
plot(shockZone(:,1),shockZone(:,2),'r','LineWidth',3)
title('PreTests PAG','FontSize',18,'FontWeight','bold');
xlim([0 1])
ylim([0 1])
hold off

axes(CondTrajPAG_Axes);
hold on
if IsMatch
    for i=idx
        plot(traj_cond_pag{i}(:,1), traj_cond_pag{i}(:,2),'LineWidth',2, 'Color', [0.8 0.8 0.8]);
    end
else
    for i=length(traj_cond_pag)
        plot(traj_cond_pag{i}(:,1), traj_cond_pag{i}(:,2),'LineWidth',2, 'Color', [0.8 0.8 0.8]);
    end
end
box off
set(gca,'XtickLabel',{},'YtickLabel',{});
plot(maze(:,1),maze(:,2),'k','LineWidth',3)
plot(shockZone(:,1),shockZone(:,2),'r','LineWidth',3)
title('Conditioning PAG','FontSize',18,'FontWeight','bold');
xlim([0 1])
ylim([0 1])
hold off

% PostTests
axes(PostTestTrajPAG_Axes);
hold on
if IsMatch
    for i=idx
        plot(traj_post_pag{i}(:,1), traj_post_pag{i}(:,2),'LineWidth',2, 'Color', [0.8 0.8 0.8]);
    end
else
    for i=length(traj_post_pag)
        plot(traj_post_pag{i}(:,1), traj_post_pag{i}(:,2),'LineWidth',2, 'Color', [0.8 0.8 0.8]);
    end
end
box off
set(gca,'XtickLabel',{},'YtickLabel',{});
plot(maze(:,1),maze(:,2),'k','LineWidth',3)
plot(shockZone(:,1),shockZone(:,2),'r','LineWidth',3)
title('PostTests PAG','FontSize',18,'FontWeight','bold');
xlim([0 1])
ylim([0 1])
hold off

%% Save it
if IsSave
    dirout = ChooseFolderForFigures_DB('Behavior');
    if IsMatch
        nameout = 'Trajectories_matched_traj';
    else
        nameout = 'Trajectories_allmice';
    end
    saveas(fh, [dirout '/' nameout '.fig']);
    saveFigure(fh,nameout,dirout);
end

end