function [EV, REV, states_sleep, states_wake, num_mice] = ExplainedVariance_master_DB(nmice, experiment, varargin)
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
%   experiment          type of experiment to analyse ('PAG' or 'MFB')
%   SpeedThresh         threshold (in cm/s) to calculate place fields
%                       (default = 4) (optional)
%   PlotResults         whether to plot results or not (optional - default=true)
%   BinSize             binsize (in tsd units) used to creat spike-time
%                       histograms (default = 0.1*1e4). In RepeatOriginalPaper
%                       mode equals 10 (1 ms) (optional)
%   SplitSleep          Splits 40 min of NREM sleep into n intervals, each
%                       SplitSleep min long. If empty, no split happens (default=20)
%   IncludeInterneurons if false, removes all interneurons (default = true)
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
    states_wake = {'Explo', 'CondMov', 'CondFreeze', 'FullTask', 'RipplesEpoch', 'PostTests'};
elseif strcmp(experiment, 'MFB')
    states_wake = {'Explo', 'CondMov', 'FullTask', 'RipplesEpoch', 'PostTests'};
elseif strcmp(experiment, 'Novel') || strcmp(experiment, 'Known')
    states_wake = {'Explo', 'CondMov', 'FullTask', 'RipplesEpoch'};
end
% Threshold on number of place cells
PCnum_thresh = 0;
% Max time to split sleep
maxtimeSplit = 40;
% Low threshold for taking out the neurons (in Hz)
lowThresh = .3;
%
save_res = 1; % Put '1' if yes.
foldertosave = ChooseFolderForFigures_DB('Data');
IsPlot = true;

% Include interneurons or not
IsII = true;
% Threshold on speed (im cm/s, epoch with lower are not considered)
speed_thresh = 4;
% Do you want to save the figures?
IsSaveFig = false;
% Binsize for firing rate histogram (in tsd units!!!). Currently Binsize = 100 ms. 1 bin = 10 ms
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
        case 'plotresults'
            IsPlot = varargin{i+1};
            if IsPlot ~= 1 && IsPlot ~= 0
                error('Incorrect value for property ''PlotResults'' (type ''help ExplainedVariance_master_DB'' for details).');
            end
        case 'splitsleep'
            splitSleep = varargin{i+1};
            if ~isa(splitSleep, 'numeric') && speed_thresh <= 0
                error('Incorrect value for property ''SplitSleep'' (type ''help ExplainedVariance_master_DB'' for details).');
            elseif splitSleep > 20
                warning(['There is going to be only one interval: ' num2str(splitSleep) ' min']);
            end
        case 'speedthresh'
            speed_thresh = varargin{i+1};
            if ~isa(speed_thresh, 'numeric') && speed_thresh <= 0
                error('Incorrect value for property ''SpeedThresh'' (type ''help ExplainedVariance_master_DB'' for details).');
            end
        case 'includeinterneurons'
            IsII = varargin{i+1};
            if ~isa(IsII, 'numeric')
                if IsII ~= 1 && IsII ~= 0
                    error('Incorrect value for property ''IncludeInterneurons'' (type ''help ExplainedVariance_master_DB'' for details).');
                end
            end
    end
end

%% Manage experiment
if strcmp(experiment, 'PAG')
    fetchpaths = 'UMazePAG';
elseif strcmp(experiment, 'MFB')
    fetchpaths = 'StimMFBWake';
elseif strcmp(experiment, 'Novel')
    fetchpaths = 'Novel';
elseif strcmp(experiment, 'Known')
    fetchpaths = 'Known';   
end

%% Allocate memory
% Get paths of each individual mouse
Dir = PathForExperimentsERC(fetchpaths);
Dir = RestrictPathForExperiment(Dir,'nMice',nmice);
numsessions = CountNumSesionsERC(Dir);

% Data
s = cell(length(numsessions), 1);
b = cell(length(numsessions), 1);
r = cell(length(numsessions), 1);
sleep = cell(length(numsessions), 1);
selected_cells = cell(length(numsessions), 1);
% Epochs
SleepEpochs = cell(length(numsessions), 1);
WakeEpochs = cell(length(numsessions), 1);
% Firing histogram
QPreSleep = cell(length(numsessions), 1);
QPostSleep = cell(length(numsessions), 1);
QWake = cell(length(numsessions), 1);
% Explained variance and reverse explained variance
EV = cell(length(states_sleep),1);
REV = cell(length(states_sleep),1);
for isleep = 1:length(EV)
    EV{isleep} = cell(length(states_wake), 1);
    REV{isleep} = cell(length(states_wake), 1);
end
for isleep = 1:length(EV)
    for iwake = 1:length(EV{isleep})
        EV{isleep}{iwake} = nan(length(numsessions), 1);
        REV{isleep}{iwake} = nan(length(numsessions), 1);
    end
end

%% Load the data
cnt=1;
for imouse=1:length(Dir.path)
    for isession = 1:length(Dir.path{imouse})
        s{cnt} = load([Dir.path{imouse}{isession} 'SpikeData.mat'],'S','PlaceCells', 'RippleGroups', 'TT', 'BasicNeuronInfo');
        b{cnt} = load([Dir.path{imouse}{isession} 'behavResources.mat'],'SessionEpoch', 'Vtsd', 'FreezeAccEpoch');
        try
            r{cnt} = load([Dir.path{imouse}{isession} 'SWR.mat'],'ripples');
        catch
            r{cnt} = load([Dir.path{imouse}{isession} 'Ripples.mat'],'ripples');
        end
        try
            sleep{cnt} = load([Dir.path{imouse}{isession} 'SleepScoring_OBGamma.mat'],'SWSEpoch','REMEpoch','Sleep');
        catch
            sleep{cnt} = load([Dir.path{imouse}{isession} 'SleepScoring_Accelero.mat'],'SWSEpoch','REMEpoch','Sleep'); % REM and Sleep are not used
        end
        cnt=cnt+1;
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
else
    split = false;
end

%% Create epoch
for isession = 1:numsessions
    
    % Get epochs
    if strcmp(experiment, 'Novel')
        [~, ~, UMazeMovEpoch, CondMovEpoch, TaskMovEpoch] = ReturnMnemozyneEpochs(b{isession}.SessionEpoch,...
            'Speed', b{isession}.Vtsd, 'SpeedThresh', speed_thresh);
        [~, ~, ~, CondEpoch, TaskEpoch, ~] = ReturnMnemozyneEpochs(b{isession}.SessionEpoch);
    else
        [~, ~, UMazeMovEpoch, CondMovEpoch, TaskMovEpoch, AfterConditioningMovEpoch] = ReturnMnemozyneEpochs(b{isession}.SessionEpoch,...
            'Speed', b{isession}.Vtsd, 'SpeedThresh', speed_thresh);
        [~, ~, ~, CondEpoch, TaskEpoch, ~] = ReturnMnemozyneEpochs(b{isession}.SessionEpoch);
    end
    RipplesEpoch = intervalSet(r{isession}.ripples(:,2)*1e4 - 1e3, r{isession}.ripples(:,2)*1e4 + 1e3);
    RipplesEpoch = and(TaskEpoch, RipplesEpoch);
    if strcmp(experiment, 'PAG')
        CondFreezeEpoch = and(CondEpoch, b{isession}.FreezeAccEpoch);
    end
    
    % Create iterable epochs
    SleepEpochs{isession} = {sleep{isession}.SWSEpoch, sleep{isession}.REMEpoch};
    if strcmp(experiment, 'PAG')
        WakeEpochs{isession} = {UMazeMovEpoch, CondMovEpoch, CondFreezeEpoch, TaskMovEpoch, RipplesEpoch, AfterConditioningMovEpoch};
    elseif strcmp(experiment, 'MFB')
        WakeEpochs{isession} = {UMazeMovEpoch, CondMovEpoch, TaskMovEpoch, RipplesEpoch, AfterConditioningMovEpoch};
    elseif strcmp(experiment, 'Novel') || strcmp(experiment, 'Known')
        WakeEpochs{isession} = {UMazeMovEpoch, CondMovEpoch, TaskMovEpoch, RipplesEpoch};
    end
    
    
    % Split epochs if necessary
    if split
        PreSWSSplitEpoch{isession} = SplitIntervals(and(b{isession}.SessionEpoch.PreSleep, sleep{isession}.SWSEpoch),...
            splitSleep*60*1e4);
        PostSWSSplitEpoch{isession} = SplitIntervals(and(b{isession}.SessionEpoch.PostSleep, sleep{isession}.SWSEpoch),...
            splitSleep*60*1e4);
    end
    
end

%% Create firing rate histograms
mask = zeros(numsessions, 1); % true at the position where mouse has more than PCnum_thresh place cells. Otherwise Qs are not created
for isession = 1:numsessions
    
%     if isfield(s{isession}.PlaceCells,'idx') && length(s{isession}.PlaceCells.idx) > PCnum_thresh
        mask(isession) = 1;
        
        % Get pyramidal cells with sustainable firing
        id_ripples = PyramLayerSorting(s{isession}.RippleGroups, s{isession}.TT);
        FR_threshold_lower = find(s{isession}.BasicNeuronInfo.firingrate >= lowThresh);
        if ~IsII
            pyr_cells_id = find(s{isession}.BasicNeuronInfo.neuroclass > 0);
            id_selected_cells = intersect(intersect(s{isession}.BasicNeuronInfo.idx_SUA, id_ripples), intersect(FR_threshold_lower, pyr_cells_id));
        else
            id_selected_cells = intersect(intersect(s{isession}.BasicNeuronInfo.idx_SUA, id_ripples), FR_threshold_lower);
        end
        
        selected_cells{isession} = s{isession}.S(id_selected_cells);
        fprintf(['\n There are ' num2str(length(id_selected_cells)) ' neurons in the analysis from the mouse ' num2str(isession) '\n']);
        
        % Bin the trains
        Q=MakeQfromS(selected_cells{isession}, binsize);
        
        QPreSleep{isession} = cell(length(states_sleep),1);
        QPostSleep{isession} = cell(length(states_sleep),1);
        % Iterate over sleep states
        for istate = 1:length(states_sleep)
            if ~strcmp(experiment, 'Novel')
                QPreSleep{isession}{istate} = zscore(full(Data(Restrict(Q,and(b{isession}.SessionEpoch.PreSleep,...
                    SleepEpochs{isession}{istate})))));
                QPostSleep{isession}{istate} = zscore(full(Data(Restrict(Q,and(b{isession}.SessionEpoch.PostSleep,...
                    SleepEpochs{isession}{istate})))));
            else
                if istate == 2
                    QPreSleep{isession}{istate} = zscore(full(Data(Restrict(Q,and(b{isession}.SessionEpoch.PreSleep,...
                        SleepEpochs{isession}{istate})))));
                    QPostSleep{isession}{istate} = zscore(full(Data(Restrict(Q,and(b{isession}.SessionEpoch.PostSleep,...
                        SleepEpochs{isession}{istate})))));
                else
                    try
                        QPreSleep{isession}{istate} = zscore(full(Data(Restrict(Q,RestrictToTime(and(b{isession}.SessionEpoch.PreSleep,...
                            SleepEpochs{isession}{istate}), 30*60*1e4)))));
                        QPostSleep{isession}{istate} = zscore(full(Data(Restrict(Q,RestrictToTime(and(b{isession}.SessionEpoch.PostSleep,...
                            SleepEpochs{isession}{istate}), 30*60*1e4)))));
                    catch
                        QPreSleep{isession}{istate} = zscore(full(Data(Restrict(Q,and(b{isession}.SessionEpoch.PreSleep,...
                            SleepEpochs{isession}{istate})))));
                        QPostSleep{isession}{istate} = zscore(full(Data(Restrict(Q,and(b{isession}.SessionEpoch.PostSleep,...
                            SleepEpochs{isession}{istate})))));
                    end
                end
            end
        end
        
        QWake{isession} = cell(length(states_wake),1);
        % Iterate over wake states
        for istate = 1:length(states_wake)
            QWake{isession}{istate} = zscore(full(Data(Restrict(Q, WakeEpochs{isession}{istate}))));
        end
        
        
        if split
            for isplit = 1:numint
                QSWSPreSleep_split{isession}{isplit} = zscore(full(Data(Restrict(Q,PreSWSSplitEpoch{isession}{isplit}))));
                QSWSPostSleep_split{isession}{isplit} = zscore(full(Data(Restrict(Q,PostSWSSplitEpoch{isession}{isplit}))));
            end
        end
        
%     end
    
end

%% Calculate EV and REV
for isession = 1:numsessions
    if mask(isession)
        % EV on full epoch
        for isleep = 1:length(states_sleep)
            for iwake = 1:length(states_wake)
                [EV{isleep}{iwake}(isession), REV{isleep}{iwake}(isession)] = ExplainedVariance(QPreSleep{isession}{isleep},QWake{isession}{iwake},QPostSleep{isession}{isleep});
            end
        end
        
        % EV on split epochs (only SWS epochs are split)
        if split
            for isplit = 1:numint
                for iwake = 1:length(states_wake)
                    [EV_SWSsplit{isplit}{iwake}(isession), REV_SWSsplit{isplit}{iwake}(isession)] = ...
                        ExplainedVariance(QSWSPreSleep_split{isession}{isplit},QWake{isession}{iwake},QSWSPostSleep_split{isession}{isplit});
                end
            end
        end
    end
end

%% Some message for the public
fprintf(['~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n',...
    '         ' num2str(sum(mask)) ' mice are in the analysis           \n' ,...
    '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n']);

num_mice = find(mask);

%% Plot
if IsPlot
    if strcmp(experiment, 'PAG')
        states_wake_toplot = states_wake(1:5); % Everything except posttests
        cols = {[.9 0 0], [.2 .2 .2]};
    elseif strcmp(experiment, 'MFB')
        states_wake_toplot = states_wake(1:4); % Everything except posttests
        cols = {[0 .9 0], [.2 .2 .2]};
    elseif strcmp(experiment, 'Novel')
        states_wake_toplot = states_wake(1:4); % Everything except posttests
        cols = {[0.793 .902 0.184], [.2 .2 .2]};
%         cols = {[1 1 1], [.2 .2 .2]};
    elseif strcmp(experiment, 'Known')
        states_wake_toplot = states_wake(1:4); % Everything except posttests
        cols = {[0.793 .902 0.184], [.2 .2 .2]};
    end
    xticks = {'EV', 'REV'};
    
    for iwake_state = 1:length(states_wake_toplot)
        for isleep_state = 1:length(states_sleep)
            fh(iwake_state, isleep_state) = figure('units', 'normalized', 'outerposition', [0 0 0.45 0.8]);
            id = find(contains(states_wake, states_wake_toplot{iwake_state}));
            data_toplot = {EV{isleep_state}{id}*100, REV{isleep_state}{id}*100};
            [b, p] = MakeBoxPlot_DB(data_toplot, cols, [1 2], xticks, 1, 'ConnectDots', false);
            ylabel('% explained');
            title([states_wake{id} ' in ' states_sleep{isleep_state}]);
            for iplot = 1:length(b)
                b{iplot}.boxAlpha = .45;
            end
            for iplot = 1:length(p)
                set(p{iplot}{1}, 'MarkerSize', 40);
            end
            p = DoWilcoxonOnArray(data_toplot, {[1 2]});
            if p < 0.05
                sigstar_DB({[1 2]},p,0,'LineWigth',16,'StarSize',24);
            end
            if strcmp(experiment, 'Novel')
                ylim([-3 45]);
            else
                ylim([-3 36.5]);
            end
            makepretty_DB
            
            if IsSaveFig
                foldertosave = ChooseFolderForFigures_DB('ReactReplay');
                if exist([foldertosave filesep 'EV'], 'dir') ~= 7
                    mkdir([foldertosave filesep 'EV']);
                end
                saveas(fh(iwake_state, isleep_state), [foldertosave filesep 'EV' filesep states_wake_toplot{iwake_state} '_' ...
                    states_sleep{isleep_state} '_' experiment '.fig']);
                saveFigure(fh(iwake_state, isleep_state), [states_wake_toplot{iwake_state} '_' states_sleep{isleep_state}...
                    '_' experiment] , [foldertosave filesep 'EV']);
            end
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
%% save data
if save_res == 1
    save([foldertosave filesep 'EV_res_pyr.mat'],'nmice', 'EV', 'REV');
end


end
