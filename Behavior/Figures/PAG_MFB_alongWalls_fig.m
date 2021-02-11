function fh = PAG_MFB_alongWalls_fig(numtest, varargin)
%
% It plots trajectories of mice from reward-learning exeriment and aversive
% learning experiment
%
%  INPUT
%
%       numtest         first <numtest> of tests will be taken to show
%                       trajectories
%       IsSave          true if you want to save the resulting figure in
%                           Dropbox (boolean - default = false)
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
% 05/02/2021
% github.com/bryzgalovdm


%%% TODO:
% * Do a statictical test

%% Parameters management
speed_thresh_mfb = 6;
speed_thresh_pag = 6;

IsSave = false;
IsMatchTraj = false;
% Optional parameters
for i=1:2:length(varargin)
    switch(lower(varargin{i}))
        case 'issave'
            IsSave = varargin{i+1};
            if IsSave ~= 1 && IsSave ~= 0
                error('Incorrect value for property ''IsSave'' (type ''help PAG_MFB_alongWalls_fig'' for details).');
            end
         case 'ismatchtraj'
            IsMatchTraj = varargin{i+1};
            if IsMatchTraj ~= 1 && IsMatchTraj ~= 0
                error('Incorrect value for property ''IsMatchTraj'' (type ''help PAG_MFB_alongWalls_fig'' for details).');
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

%% Calculate occupancy
% MFB data
cnt_pre=1;
cnt_cond=1;
cnt_post = 1;
for idir = 1:length(DirMFB.path)
    for itest = 1:numtest
        % Prepare variables
        VtsdSmoothed = tsd(Range(dat_MFB{idir}.behavResources(id_Pre_MFB{idir}(itest)).Vtsd),...
            movmedian(Data(dat_MFB{idir}.behavResources(id_Pre_MFB{idir}(itest)).Vtsd),5));
        LocomotionEpoch = thresholdIntervals(VtsdSmoothed,speed_thresh_mfb,'Direction','Above');
        MovingX = Restrict(dat_MFB{idir}.behavResources(id_Pre_MFB{idir}(itest)).Xtsd, LocomotionEpoch);
        MovingY = Restrict(dat_MFB{idir}.behavResources(id_Pre_MFB{idir}(itest)).Ytsd, LocomotionEpoch);
        % Calculate occupancies
        [~,~,ZoneIndices] = RecalcInOutZones(dat_MFB{idir}.behavResources(id_Pre_MFB{idir}(itest)).mask, MovingX, MovingY,...
            dat_MFB{idir}.behavResources(id_Pre_MFB{idir}(itest)).Ratio_IMAonREAL, 'Ratio', 0.33);
        for izone = 1:length(ZoneIndices)
            occup_pre_mfb(cnt_pre,izone) = size(ZoneIndices{izone},1)...
                ./size(Data(Restrict(dat_MFB{idir}.behavResources(id_Pre_MFB{idir}(itest)).AlignedXtsd, LocomotionEpoch)), 1);
        end
        cnt_pre=cnt_pre+1;
    end
    for itest = 1:numtest
        % Prepare variables
        VtsdSmoothed = tsd(Range(dat_MFB{idir}.behavResources(id_Cond_MFB{idir}(itest)).Vtsd),...
            movmedian(Data(dat_MFB{idir}.behavResources(id_Cond_MFB{idir}(itest)).Vtsd),5));
        LocomotionEpoch = thresholdIntervals(VtsdSmoothed,speed_thresh_mfb,'Direction','Above');
        MovingX = Restrict(dat_MFB{idir}.behavResources(id_Cond_MFB{idir}(itest)).Xtsd, LocomotionEpoch);
        MovingY = Restrict(dat_MFB{idir}.behavResources(id_Cond_MFB{idir}(itest)).Ytsd, LocomotionEpoch);
        % Calculate occupancies
        [~,~,ZoneIndices] = RecalcInOutZones(dat_MFB{idir}.behavResources(id_Cond_MFB{idir}(itest)).mask, MovingX, MovingY, ...
            dat_MFB{idir}.behavResources(id_Cond_MFB{idir}(itest)).Ratio_IMAonREAL, 'Ratio', 0.33);
        for izone = 1:length(ZoneIndices)
            occup_cond_mfb(cnt_cond, izone) = size(ZoneIndices{izone},1)...
                ./size(Data(Restrict(dat_MFB{idir}.behavResources(id_Cond_MFB{idir}(itest)).AlignedXtsd, LocomotionEpoch)), 1);
        end
        cnt_cond=cnt_cond+1;
    end
    for itest = 1:numtest
        % Prepare variables
        VtsdSmoothed = tsd(Range(dat_MFB{idir}.behavResources(id_Post_MFB{idir}(itest)).Vtsd),...
            movmedian(Data(dat_MFB{idir}.behavResources(id_Post_MFB{idir}(itest)).Vtsd),5));
        LocomotionEpoch = thresholdIntervals(VtsdSmoothed,speed_thresh_mfb,'Direction','Above');
        MovingX = Restrict(dat_MFB{idir}.behavResources(id_Post_MFB{idir}(itest)).Xtsd, LocomotionEpoch);
        MovingY = Restrict(dat_MFB{idir}.behavResources(id_Post_MFB{idir}(itest)).Ytsd, LocomotionEpoch);
        % Calculate occupancies
        [~,~,ZoneIndices] = RecalcInOutZones(dat_MFB{idir}.behavResources(id_Post_MFB{idir}(itest)).mask,  MovingX, MovingY, ...
            dat_MFB{idir}.behavResources(id_Post_MFB{idir}(itest)).Ratio_IMAonREAL, 'Ratio', 0.33);
        for izone = 1:length(ZoneIndices)
            occup_post_mfb(cnt_post,izone) = size(ZoneIndices{izone},1)...
                ./size(Data(Restrict(dat_MFB{idir}.behavResources(id_Post_MFB{idir}(itest)).AlignedXtsd, LocomotionEpoch)), 1);
        end
        cnt_post=cnt_post+1;
    end
    
end

% PAG data
cnt_pre=1;
cnt_cond=1;
cnt_post=1;
for idir = 1:length(DirPAG.path)
    
    % Store in the variable
    for itest = 1:numtest
        % Prepare variables
        VtsdSmoothed = tsd(Range(dat_PAG{idir}.behavResources(id_Pre_PAG{idir}(itest)).CleanVtsd),...
            movmedian(Data(dat_PAG{idir}.behavResources(id_Pre_PAG{idir}(itest)).CleanVtsd),5));
        LocomotionEpoch = thresholdIntervals(VtsdSmoothed,speed_thresh_pag,'Direction','Above');
        MovingX = Restrict(dat_PAG{idir}.behavResources(id_Pre_PAG{idir}(itest)).Xtsd, LocomotionEpoch);
        MovingY = Restrict(dat_PAG{idir}.behavResources(id_Pre_PAG{idir}(itest)).Ytsd, LocomotionEpoch);
        % Calculate occupancies
        [~,~,ZoneIndices] = RecalcInOutZones(dat_PAG{idir}.behavResources(id_Pre_PAG{idir}(itest)).mask, MovingX, MovingY, ...
            dat_PAG{idir}.behavResources(id_Pre_PAG{idir}(itest)).Ratio_IMAonREAL, 'Ratio', 0.33);
        for izone = 1:length(ZoneIndices)
            occup_pre_pag(cnt_pre, izone) = size(ZoneIndices{izone},1)...
                ./size(Data(Restrict(dat_PAG{idir}.behavResources(id_Pre_PAG{idir}(itest)).CleanAlignedXtsd, LocomotionEpoch)), 1);
        end
        cnt_pre=cnt_pre+1;
    end
    for itest = 1:numtest
        % Prepare variables
        VtsdSmoothed = tsd(Range(dat_PAG{idir}.behavResources(id_Cond_PAG{idir}(itest)).CleanVtsd),...
            movmedian(Data(dat_PAG{idir}.behavResources(id_Cond_PAG{idir}(itest)).CleanVtsd),5));
        LocomotionEpoch = thresholdIntervals(VtsdSmoothed,speed_thresh_pag,'Direction','Above');
        MovingX = Restrict(dat_PAG{idir}.behavResources(id_Cond_PAG{idir}(itest)).Xtsd, LocomotionEpoch);
        MovingY = Restrict(dat_PAG{idir}.behavResources(id_Cond_PAG{idir}(itest)).Ytsd, LocomotionEpoch);
        % Calculate occupancies
        [~,~,ZoneIndices] = RecalcInOutZones(dat_PAG{idir}.behavResources(id_Cond_PAG{idir}(itest)).mask, MovingX, MovingY, ...
            dat_PAG{idir}.behavResources(id_Cond_PAG{idir}(itest)).Ratio_IMAonREAL, 'Ratio', 0.33);
        for izone = 1:length(ZoneIndices)
            occup_cond_pag(cnt_cond, izone) = size(ZoneIndices{izone},1)...
                ./size(Data(Restrict(dat_PAG{idir}.behavResources(id_Cond_PAG{idir}(itest)).CleanAlignedXtsd, LocomotionEpoch)), 1);
        end
        cnt_cond=cnt_cond+1;
    end
    for itest = 1:numtest
        % Prepare variables
        VtsdSmoothed = tsd(Range(dat_PAG{idir}.behavResources(id_Post_PAG{idir}(itest)).CleanVtsd),...
            movmedian(Data(dat_PAG{idir}.behavResources(id_Post_PAG{idir}(itest)).CleanVtsd),5));
        LocomotionEpoch = thresholdIntervals(VtsdSmoothed,speed_thresh_pag,'Direction','Above');
        MovingX = Restrict(dat_PAG{idir}.behavResources(id_Post_PAG{idir}(itest)).Xtsd, LocomotionEpoch);
        MovingY = Restrict(dat_PAG{idir}.behavResources(id_Post_PAG{idir}(itest)).Ytsd, LocomotionEpoch);
        % Calculate occupancies
        [~,~,ZoneIndices] = RecalcInOutZones(dat_PAG{idir}.behavResources(id_Post_PAG{idir}(itest)).mask,  MovingX, MovingY, ...
            dat_PAG{idir}.behavResources(id_Post_PAG{idir}(itest)).Ratio_IMAonREAL, 'Ratio', 0.33);
        for izone = 1:length(ZoneIndices)
            occup_post_pag(cnt_post, izone) = size(ZoneIndices{izone},1)...
                ./size(Data(Restrict(dat_PAG{idir}.behavResources(id_Post_PAG{idir}(itest)).CleanAlignedXtsd, LocomotionEpoch)), 1);
        end
        cnt_post=cnt_post+1;
    end
    
end

%% Create arrays to plot - Ratio Walls/Center
comp_pre_mfb = occup_pre_mfb(:,2)./occup_pre_mfb(:,1);
comp_pre_pag = occup_pre_pag(:,2)./occup_pre_pag(:,1);

comp_cond_mfb = occup_cond_mfb(:,2)./occup_cond_mfb(:,1);
comp_cond_pag = occup_cond_pag(:,2)./occup_cond_pag(:,1);

comp_post_mfb = occup_post_mfb(:,2)./occup_post_mfb(:,1);
comp_post_pag = occup_post_pag(:,2)./occup_post_pag(:,1);

%% Set up the figure

% Axes
fh = figure('units', 'normalized', 'outerposition', [0 0 0.7 0.6]);
PreTest_Axes = axes('position', [0.07 0.07 0.25 0.8]);
CondAxes = axes('position', [0.38 0.07 0.25 0.8]);
PostTest_Axes = axes('position', [0.69 0.07 0.25 0.8]);
cols = {[0.6 1 0.6], [1 0.6 0.6]};

if IsMatchTraj
    idx = datasample(1:length(comp_pre_pag), length(comp_pre_mfb), 'Replace', false);
end

% Plot
% PreTests
axes(PreTest_Axes);
if ~IsMatchTraj
    MakeBoxPlot_DB({comp_pre_mfb, comp_pre_pag}, cols, 1:2, {'MFB', 'PAG'},0);
else
    MakeBoxPlot_DB({comp_pre_mfb, comp_pre_pag(idx)}, cols, 1:2, {'MFB', 'PAG'},0);
    % Rank sum test
    p = ranksum(comp_pre_mfb,comp_pre_pag(idx));
    if p<0.05
        sigstar_DB({[1 2]},p,0,'LineWigth',16,'StarSize',24);
    end
end
ylabel('Walls/Center ratio');
box off
title('PreTests','FontSize',18,'FontWeight','bold');
if IsMatchTraj
    ylim([0 8]) 
else
    ylim([0 6])
end
makepretty

axes(CondAxes);
if ~IsMatchTraj
    MakeBoxPlot_DB({comp_cond_mfb, comp_cond_pag}, cols, 1:2, {'MFB', 'PAG'},0);
else
    MakeBoxPlot_DB({comp_cond_mfb, comp_cond_pag(idx)}, cols, 1:2, {'MFB', 'PAG'},0);
    % Rank sum test
    p = ranksum(comp_cond_mfb,comp_cond_pag(idx));
    if p<0.05
        sigstar_DB({[1 2]},p,0,'LineWigth',16,'StarSize',24);
    end
end

ylabel('Walls/Center ratio');
box off
title('Conditioning','FontSize',18,'FontWeight','bold');
if IsMatchTraj
    ylim([0 8]) 
else
    ylim([0 6])
end
makepretty

% PostTests
axes(PostTest_Axes);
if ~IsMatchTraj
    MakeBoxPlot_DB({comp_post_mfb, comp_post_pag}, cols, 1:2, {'MFB', 'PAG'},0);
else
    MakeBoxPlot_DB({comp_post_mfb, comp_post_pag(idx)}, cols, 1:2, {'MFB', 'PAG'},0);
    % Rank sum test
    p = ranksum(comp_post_mfb,comp_post_pag(idx));
    if p<0.05
        sigstar_DB({[1 2]},p,0,'LineWigth',16,'StarSize',24);
    end
end
ylabel('Walls/Center ratio');
box off
title('PostTests','FontSize',18,'FontWeight','bold');
if IsMatchTraj
    ylim([0 8]) 
else
    ylim([0 6])
end
makepretty

%% Save it
if IsSave
    dirout = ChooseFolderForFigures_DB('Behavior');
    nameout = 'Along_theWalls_ratio';
    saveas(fh, [dirout '/' nameout '.fig']);
    saveFigure(fh,nameout,dirout);
end

end