function [EV, REV, states_sleep, states_wake] = ExplainedVariance_master_DB(nmice, experiment, varargin)
%
% The function calculates explained variance (EV) and reverse EV
% the same manner it was calculated in Kudrimoti et al., 1999
%
% EV = ((R_task,post - R_task,pre*R_post,pre)/sqrt((1-R_task,pre.^2)(1-R_post,pre.^2))).^2;
%
% REV = ((R_task,pre - R_task,post*R_post,pre)/sqrt((1-R_task,post.^2)(1-R_post,pre.^2))).^2;
%
%
% INPUT
%
%   nmice               array with mice numbers that will harvested from ERC
%                       PathForExperiments. Each mouse should contain
%                       PlaceCells structure in its SpikeData.mat
%   experiment          tyoe of experiment to analyse ('PAG' or 'MFB')
%   SpeedThresh         threshold (in cm/s) to calculate place fields
%                       (default = 4) (optional)
%   BinSize             binsize (in tsd units) used to creat spike-time
%                       histograms (default = 0.1*1e4). In RepeatOriginalPaper
%                       mode equals 10 (1 ms) (optional)
%   OverlapFactor       defines how many pixels in place fields of two neurons
%                       should overlap to be considered overlapping
%                       (default = 25) (optional)
%   SplitSleep          Splits 40 min of NREM sleep into n intervals, each
%                       SplitSleep min long. If empty, no split happens (default=20)
%   HighThresh          if not empty, removes all neurons that fire more
%                       than HighThresh Hz (default = 5)
%   IsSaveFig           if true, saves figures in dropbox (default = false)(optional)
%
% OUTPUT
%
%  None
%
% EXAMPLE
%
%   [EV, REV] = ExplainedVariance_master_DB([797 798 1117]);
%   PairwiseCorrelationsAcrossMice_WM1994_DB([797 798 1117], 'SaveFigure', true, 'BinSize', 0.05*1e4);
%
% By Dima Bryzgalov, MOBS team, Paris,
% 01/04/2021
% github.com/bryzgalovdm

%% Default optional parameters
% States
states_sleep = {'NREM', 'REM'};
if strcmp(experiment, 'PAG')
    states_wake = {'Explo', 'CondMov', 'CondFreeze', 'FullTask', 'PostTests'};
elseif strcmp(experiment, 'MFB')
    states_wake = {'Explo', 'CondMov', 'FullTask', 'PostTests'};
end
% Threshold on number of place cells
PCnum_thresh = 2;
% Max time to split sleep
maxtimeSplit = 40;
% Low threshold for taking out the neurons (in Hz)
lowThresh = .3;

% High threshold for taking out the neurons (in Hz)
highThresh = 5;
% Threshold on speed (im cm/s, epoch with lower are not considered)
speed_thresh = 4;
% Do you want to save the figures?
IsSaveFig = false; %%% Does not work now
% Binsize for firing rate histogram (in tsd units!!!)
binsize = 0.1*1e4;
% Do you want to split NREM sleep epochs into several consecutive intervals (in min)?
% If not - put []. Note that number of intervals will be floor(40/splitSleep)
% It cannot be more than 40
splitSleep = 20;

%% Parse parameters
for i=1:2:length(varargin)
    switch(lower(varargin{i}))
        case 'binsize'
            binsize = varargin{i+1};
            if ~isa(binsize, 'numeric') && speed_thresh <= 0
                error('Incorrect value for property ''BinSize'' (type ''help ExplainedVariance_master_DB'' for details).');
            end
        case 'issavefig'
            IsSaveFig = varargin{i+1};
            if IsSaveFig ~= 1 && IsSaveFig ~= 0
                error('Incorrect value for property ''IsSaveFig'' (type ''help ExplainedVariance_master_DB'' for details).');
            end
        case 'splitsleep'
            splitSleep = varargin{i+1};
            if ~isa(splitSleep, 'numeric') && speed_thresh <= 0
                error('Incorrect value for property ''BinSize'' (type ''help ExplainedVariance_master_DB'' for details).');
            elseif splitSleep > 20
                warning(['There is going to be only one interval: ' num2str(splitSleep) ' min']);
            end
        case 'speedthresh'
            speed_thresh = varargin{i+1};
            if ~isa(speed_thresh, 'numeric') && speed_thresh <= 0
                error('Incorrect value for property ''SpeedThresh'' (type ''help ExplainedVariance_master_DB'' for details).');
            end
        case 'hightresh'
            highThresh = varargin{i+1};
            if ~isa(speed_thresh, 'numeric')
                if ~isempty(highTresh) && highThresh <= 0
                    error('Incorrect value for property ''HighThresh'' (type ''help ExplainedVariance_master_DB'' for details).');
                end
            end
    end
end

%% Manage experiment
if strcmp(experiment, 'PAG')
    fetchpaths = 'UMazePAG';
elseif strcmp(experiment, 'MFB')
    fetchpaths = 'StimMFBWake';
end

%% Allocate memory
% Data
s = cell(length(nmice), 1);
b = cell(length(nmice), 1);
r = cell(length(nmice), 1);
sleep = cell(length(nmice), 1);
selected_cells = cell(length(nmice), 1);
% Epochs
SleepEpochs = cell(length(nmice), 1);
WakeEpochs = cell(length(nmice), 1);
% Firing histogram
QPreSleep = cell(length(nmice), 1);
QPostSleep = cell(length(nmice), 1);
QWake = cell(length(nmice), 1);
% Explained variance and reverse explained variance
EV = cell(length(states_sleep),1);
REV = cell(length(states_sleep),1);
for isleep = 1:length(EV)
    EV{isleep} = cell(length(states_wake), 1);
    REV{isleep} = cell(length(states_wake), 1);
end
for isleep = 1:length(EV)
    for iwake = 1:length(EV{isleep})
        EV{isleep}{iwake} = nan(length(nmice), 1);
        REV{isleep}{iwake} = nan(length(nmice), 1);
    end
end

%% Load the data
% Get paths of each individual mouse
Dir = PathForExperimentsERC(fetchpaths);
Dir = RestrictPathForExperiment(Dir,'nMice',nmice);

for imouse=1:length(Dir.path)
    s{imouse} = load([Dir.path{imouse}{1} 'SpikeData.mat'],'S','PlaceCells', 'RippleGroups', 'TT', 'BasicNeuronInfo');
    b{imouse} = load([Dir.path{imouse}{1} 'behavResources.mat'],'SessionEpoch', 'CleanVtsd', 'FreezeAccEpoch');
    try
        r{imouse} = load([Dir.path{imouse}{1} 'SWR.mat'],'ripples');
    catch
        r{imouse} = load([Dir.path{imouse}{1} 'Ripples.mat'],'ripples');
    end
    try
        sleep{imouse} = load([Dir.path{imouse}{1} 'SleepScoring_OBGamma.mat'],'SWSEpoch','REMEpoch','Sleep');
    catch
        sleep{imouse} = load([Dir.path{imouse}{1} 'SleepScoring_Accelero.mat'],'SWSEpoch','REMEpoch','Sleep'); % REM and Sleep are not used
    end
end

%% Manage splits
if ~isempty(splitSleep) && splitSleep <= 20
    split=true;
    numint = floor(maxtimeSplit/splitSleep);
    % Allocate epochs
    PreSWSSplitEpoch = cell(length(nmice), 1);
    PostSWSSplitEpoch = cell(length(nmice), 1);
    % Allocate FR histos
    QSWSPreSleep_split = cell(length(nmice), 1);
    QSWSPostSleep_split = cell(length(nmice), 1);
    % Allocate EV
    EV_SWSsplit = cell(numint, 1);
    REV_SWSsplit = cell(numint, 1);
    for isplit = 1:length(EV_SWSsplit)
        EV_SWSsplit{isplit} = cell(length(states_wake), 1);
        REV_SWSsplit{isplit} = cell(length(states_wake), 1);
    end
    for isplit = 1:length(EV_SWSsplit)
        for iwake = 1:length(EV{isplit})
            EV_SWSsplit{isplit}{iwake} = nan(length(nmice), 1);
            REV_SWSsplit{isplit}{iwake} = nan(length(nmice), 1);
        end
    end
end

%% Create epoch
for imouse = 1:length(Dir.path)
    
    % Get epochs
    [~, UMazeMovEpoch, CondMovEpoch, TaskMovEpoch, AfterConditioningMovEpoch] = ReturnMnemozyneEpochs(b{imouse}.SessionEpoch,...
        'Speed', b{imouse}.CleanVtsd, 'SpeedThresh', speed_thresh);
    if strcmp(experiment, 'PAG')
        [~, ~, CondEpoch, ~, ~] = ReturnMnemozyneEpochs(b{imouse}.SessionEpoch);
        CondFreezeEpoch = and(CondEpoch, b{imouse}.FreezeAccEpoch);
    end
    
    % Create iterable epochs
    SleepEpochs{imouse} = {sleep{imouse}.SWSEpoch, sleep{imouse}.REMEpoch};
    if strcmp(experiment, 'PAG')
        WakeEpochs{imouse} = {UMazeMovEpoch, CondMovEpoch, CondFreezeEpoch, TaskMovEpoch, AfterConditioningMovEpoch};
    elseif strcmp(experiment, 'MFB')
        WakeEpochs{imouse} = {UMazeMovEpoch, CondMovEpoch, TaskMovEpoch, AfterConditioningMovEpoch};
    end
    
    % Split epochs if necessary
    if split
        PreSWSSplitEpoch{imouse} = SplitIntervals(and(b{imouse}.SessionEpoch.PreSleep, sleep{imouse}.SWSEpoch),...
            splitSleep*60*1e4);
        PostSWSSplitEpoch{imouse} = SplitIntervals(and(b{imouse}.SessionEpoch.PostSleep, sleep{imouse}.SWSEpoch),...
            splitSleep*60*1e4);
    end
    
end

%% Create firing rate histograms
mask = zeros(length(nmice), 1); % true at the position where mouse has more than 2 place cells. Otherwise Qs are not created
for imouse = 1:length(Dir.path)
    
    if isfield(s{imouse}.PlaceCells,'idx') && length(s{imouse}.PlaceCells.idx) > PCnum_thresh
        mask(imouse) = 1;
        
        % Get pyramidal cells with sustainable firing
        id_ripples = PyramLayerSorting(s{imouse}.S, s{imouse}.RippleGroups, s{imouse}.TT);
        FR_threshold_lower = find(s{imouse}.BasicNeuronInfo.firingrate >= lowThresh);
        if ~isempty(highThresh)
            FR_threshold_upper = find(s{imouse}.BasicNeuronInfo.firingrate <= highThresh);
            id_selected_cells = intersect(intersect(s{imouse}.BasicNeuronInfo.idx_SUA, id_ripples), intersect(FR_threshold_lower, FR_threshold_upper));
        else
            id_selected_cells = intersect(intersect(s{imouse}.BasicNeuronInfo.idx_SUA, id_ripples), FR_threshold_lower);
        end
        
        selected_cells{imouse} = s{imouse}.S(id_selected_cells);
        fprintf(['\n There are ' num2str(length(id_selected_cells)) ' neurons in the analysis from the mouse ' num2str(Dir.name{imouse}) '\n']);
        
        % Bin the trains
        Q=MakeQfromS(selected_cells{imouse}, binsize);
        
        QPreSleep{imouse} = cell(length(states_sleep),1);
        QPostSleep{imouse} = cell(length(states_sleep),1);
        % Iterate over sleep states
        for istate = 1:length(states_sleep)
            QPreSleep{imouse}{istate} = zscore(full(Data(Restrict(Q,and(b{imouse}.SessionEpoch.PreSleep, SleepEpochs{imouse}{istate})))));
            QPostSleep{imouse}{istate} = zscore(full(Data(Restrict(Q,and(b{imouse}.SessionEpoch.PostSleep, SleepEpochs{imouse}{istate})))));
        end
        
        QWake{imouse} = cell(length(states_wake),1);
        % Iterate over wake states
        for istate = 1:length(states_wake)
            QWake{imouse}{istate} = zscore(full(Data(Restrict(Q, WakeEpochs{imouse}{istate}))));
        end
        
        
        if split
            for isplit = 1:numint
                QSWSPreSleep_split{imouse}{isplit} = zscore(full(Data(Restrict(Q,PreSWSSplitEpoch{imouse}{isplit}))));
                QSWSPostSleep_split{imouse}{isplit} = zscore(full(Data(Restrict(Q,PostSWSSplitEpoch{imouse}{isplit}))));
            end
        end
        
    end
    
end

%% Calculate EV and REV
for imouse = 1:length(nmice)
    if mask(imouse)
        % EV on full epoch
        for isleep = 1:length(states_sleep)
            for iwake = 1:length(states_wake)
                [EV{isleep}{iwake}(imouse), REV{isleep}{iwake}(imouse)] = ExplainedVariance(QPreSleep{imouse}{isleep},QWake{imouse}{iwake},QPostSleep{imouse}{isleep});
            end
        end
        
        % EV on split epochs (only SWS epochs are split)
        if split
            for isplit = 1:numint
                for iwake = 1:length(states_wake)
                    [EV_SWSsplit{isplit}{iwake}(imouse), REV_SWSsplit{isplit}{iwake}(imouse)] = ...
                        ExplainedVariance(QSWSPreSleep_split{imouse}{isplit},QWake{imouse}{iwake},QSWSPostSleep_split{imouse}{isplit});
                end
            end
        end
    end
end

%% Some message for the public
fprintf(['~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n',...
    '         ' num2str(sum(mask)) ' mice are in the analysis           \n' ,...
    '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n']);

%% Plot
if strcmp(experiment, 'PAG')
    states_wake_toplot = states_wake(1:4); % Everything except posttests
elseif strcmp(experiment, 'MFB')
    states_wake_toplot = states_wake(1:3); % Everything except posttests
end
cols = {[.8 .8 .8], [.2 .2 .2]};
xticks = {'EV', 'REV'};

for iwake_state = 1:length(states_wake_toplot)
    fh(iwake_state) = figure('units', 'normalized', 'outerposition', [0 0 0.4 0.5]);
    for isleep_state = 1:length(states_sleep)
        subplot(1, 2, isleep_state)
        id = find(contains(states_wake, states_wake_toplot{iwake_state}));
        data_toplot = {EV{isleep_state}{id}*100, REV{isleep_state}{id}*100};
        MakeBoxPlot_DB(data_toplot, cols, [1 2], xticks, 1, 'ConnectDots', 1);
        ylabel('% explained');
        title([states_wake{id} ' in ' states_sleep{isleep_state}]);
        
        p = DoWilcoxonOnArray(data_toplot, {[1 2]});
        if p < 0.05
            sigstar_DB({[1 2]},p,0,'LineWigth',16,'StarSize',24);
        end
        makepretty
    end
    
    if IsSaveFig
        foldertosave = ChooseFolderForFigures_DB('ReactReplay');
        if exist([foldertosave filesep 'EV'], 'dir') ~= 7
            mkdir([foldertosave filesep 'EV']);
        end
        saveas(fh(iwake_state), [foldertosave filesep 'EV' filesep states_wake_toplot{iwake_state} '_' experiment '.fig']);
        saveFigure(fh(iwake_state), [states_wake_toplot{iwake_state} '_' experiment] , [foldertosave filesep 'EV']);
    end
end

if split
    xt = union([1:3:3*numint], [2:3:3*numint+1]);
    for iwake = 1:length(states_wake_toplot)
        fsh(iwake) = figure('units', 'normalized', 'outerposition', [0 0 0.6 0.5]);
        
        id = find(contains(states_wake, states_wake_toplot{iwake}));
        data_toplot = {EV_SWSsplit{1}{id}*100, REV_SWSsplit{1}{id}*100};
        for isplit = 2:numint
            data_toplot{isplit*2-1} = EV_SWSsplit{isplit}{id}*100;
            data_toplot{isplit*2} = REV_SWSsplit{isplit}{id}*100;
        end
        MakeBoxPlot_DB(data_toplot, repmat(cols,1,numint), xt,  repmat({'EV', 'REV'}, 1, numint), 1, 'ConnectDots', 1);
        ylabel('% explained');
        title(['EV of ' states_wake{id} ' in SWS split in periods']);
        
        pairs_toCompare = {[1 2], [3 4]};
        p = DoWilcoxonOnArray(data_toplot, pairs_toCompare);
        pairs_stars = arrayfun(@(x) xt([x*2-1 x*2]), 1:numint, 'UniformOutput', false);
        for ip = 1:length(p)
            if p(ip) <= 0.05
                sigstar_DB(pairs_stars(ip),p(ip),0,'LineWigth',16,'StarSize',24);
            end
        end
        makepretty
        
        if IsSaveFig
            foldertosave = ChooseFolderForFigures_DB('ReactReplay');
            if exist([foldertosave filesep 'EV'], 'dir') ~= 7
                mkdir([foldertosave filesep 'EV']);
            end
            saveas(fsh(iwake), [foldertosave filesep 'EV' filesep states_wake_toplot{iwake} '_' experiment '_split.fig']);
            saveFigure(fsh(iwake), [states_wake_toplot{iwake} '_' experiment '_split'], [foldertosave filesep 'EV']);
        end
        
    end
end

end
