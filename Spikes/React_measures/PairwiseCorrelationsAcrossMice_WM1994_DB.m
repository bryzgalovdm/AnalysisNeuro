%%%% PairwiseCorrelationsAcrossMice_WM1994_DB
%
% This script performs pairwise correlations on all overlapping
% and non-overlapping place cells of each particular mouse,
% averages them out and then does global mean across all mice
% 
% It repeats the analysis of Wilson&McNaughton, 1994, Science
% Clearly explained in Tingley&Peyrache,2020, Phil.Trans.B.
% 
% Please go inside the script and check the parameters
% 
% By Dima Bryzgalov, MOBS team, Paris, 
% 01/05/2020
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

% Threshold on speed (im cm/s, epoch with lower are not considered)
speed_thresh = 4;

% Mice with how many PCs are taken into analysis?
PCnum_thresh = 2;

% Parameters of cross-correlograms
binsize = 0.1*1e4; % (measured in tsd units!!!);

% How many pixels should overlap to define an overlap in the place field?
overlapFactor = 25;

% Do you want to correlate firing histograms (false) or short latency cross-correlations (true)?
do_as_in_paper = false; % true is very-very slow. Use with care

% Do you want to save the figures?
savfig = true;

% Do you want to plot mean over mice or over pooled neurons? (false - default - over pooled neurons)
plotmice = false; 

% Do you want to plot each single cross correlation of overlapping pairs?
plotexamples = true;

% Do you want to restrict your sleep time? If not put []
SleepTimeToRestrict = 20*60*1e4;

%% Allocate memory
spikes = cell(1,length(Dir.path));
behav = cell(1,length(Dir.path));
sleep = cell(1,length(Dir.path));
UMazeEpoch = cell(1,length(Dir.path));
ConditioningEpoch = cell(1, length(Dir.path));
AfterConditioningEpoch = cell(1, length(Dir.path));
LocomotionEpoch = cell(1, length(Dir.path));
UMazeMovingEpoch = cell(1, length(Dir.path));
AfterConditioningMovingEpoch = cell(1, length(Dir.path));
ConditioningMovingEpoch = cell(1, length(Dir.path));
ConditioningFreezingEpoch = cell(1, length(Dir.path));
PreSleepFinal = cell(1, length(Dir.path));
PostSleepFinal = cell(1, length(Dir.path));
overlappairs = cell(1, length(Dir.path));
distantpairs = cell(1, length(Dir.path));
QPRE = cell(1, length(Dir.path));
QTASK = cell(1, length(Dir.path));
QPOST = cell(1, length(Dir.path));
QCONDMOV = cell(1, length(Dir.path));
QCONDFREEZ = cell(1, length(Dir.path));
QPOSTTEST = cell(1, length(Dir.path));
rho_overlap = cell(1, length(Dir.path));
rho_distant = cell(1, length(Dir.path));

%% Load Data
for i=1:length(Dir.path)
    
    spikes{i} = load([Dir.path{i}{1} 'SpikeData.mat'],'S','PlaceCells');
    behav{i} = load([Dir.path{i}{1} 'behavResources.mat'],'SessionEpoch','CleanVtsd','CleanAlignedXtsd','CleanAlignedYtsd','FreezeAccEpoch');
    if strcmp(Dir.name{i}, 'Mouse906') || strcmp(Dir.name{i}, 'Mouse977') % Mice with bad OB-based sleep scoring
        sleep{i} = load([Dir.path{i}{1} 'SleepScoring_Accelero.mat'],'SWSEpoch','REMEpoch','Sleep'); % REM and Sleep are not used
    else
        sleep{i} = load([Dir.path{i}{1} 'SleepScoring_OBGamma.mat'],'SWSEpoch','REMEpoch','Sleep');  % REM and Sleep are not used
    end
    
    %% Epochs
    
    % BaselineExplo Epoch
    UMazeEpoch{i} = or(behav{i}.SessionEpoch.Hab,behav{i}.SessionEpoch.TestPre1);
    UMazeEpoch{i} = or(UMazeEpoch{i},behav{i}.SessionEpoch.TestPre2);
    UMazeEpoch{i} = or(UMazeEpoch{i},behav{i}.SessionEpoch.TestPre3);
    UMazeEpoch{i} = or(UMazeEpoch{i},behav{i}.SessionEpoch.TestPre4);
    
    % Conditioning
    ConditioningEpoch{i} = or(behav{i}.SessionEpoch.Cond1,behav{i}.SessionEpoch.Cond2);
    ConditioningEpoch{i} = or(ConditioningEpoch{i},behav{i}.SessionEpoch.Cond3);
    ConditioningEpoch{i} = or(ConditioningEpoch{i},behav{i}.SessionEpoch.Cond4);
    
    
    % After Conditioning
    AfterConditioningEpoch{i} = or(behav{i}.SessionEpoch.TestPost1,behav{i}.SessionEpoch.TestPost2);
    AfterConditioningEpoch{i} = or(AfterConditioningEpoch{i},behav{i}.SessionEpoch.TestPost3);
    AfterConditioningEpoch{i} = or(AfterConditioningEpoch{i},behav{i}.SessionEpoch.TestPost4);
    
    % Locomotion threshold
    VtsdSmoothed  = tsd(Range(behav{i}.CleanVtsd),movmedian(Data(behav{i}.CleanVtsd),5)); % SmoFac = 5
    LocomotionEpoch{i} = thresholdIntervals(VtsdSmoothed,speed_thresh,'Direction','Above');
    
    % Get resulting epoch
    UMazeMovingEpoch{i} = and(LocomotionEpoch{i}, UMazeEpoch{i});
    AfterConditioningMovingEpoch{i} = and(LocomotionEpoch{i}, AfterConditioningEpoch{i});
    ConditioningMovingEpoch{i} = and(LocomotionEpoch{i}, ConditioningEpoch{i});
    ConditioningFreezingEpoch{i} = and(behav{i}.FreezeAccEpoch, ConditioningEpoch{i});
    
    % SleepEpochs
    if isempty(SleepTimeToRestrict)
        PreSleepFinal{i} = and(behav{i}.SessionEpoch.PreSleep, sleep{i}.SWSEpoch);
        PostSleepFinal{i} = and(behav{i}.SessionEpoch.PostSleep, sleep{i}.SWSEpoch);
    else
        PreSleepFinal{i} = RestrictToTime(and(behav{i}.SessionEpoch.PreSleep, sleep{i}.SWSEpoch),SleepTimeToRestrict);
        PostSleepFinal{i} = RestrictToTime(and(behav{i}.SessionEpoch.PostSleep, sleep{i}.SWSEpoch),SleepTimeToRestrict);
    end
    
    
    %% Create binned FR vectors
    if do_as_in_paper
        
        Q = MakeQfromS(spikes{i}.S, 10); % 1ms nbin
        
        QPRE{i} = zscore(full(Data(Restrict(Q, PreSleepFinal{i}))));
        QTASK{i} = zscore(full(Data(Restrict(Q, UMazeMovingEpoch{i}))));
        QPOST{i} = zscore(full(Data(Restrict(Q, PostSleepFinal{i}))));
        QCONDMOV{i} = zscore(full(Data(Restrict(Q, ConditioningMovingEpoch{i}))));
        QCONDFREEZ{i} = zscore(full(Data(Restrict(Q, ConditioningFreezingEpoch{i}))));
        QPOSTTEST{i} = zscore(full(Data(Restrict(Q, AfterConditioningMovingEpoch{i}))));
        
        clear Q
        
    else
        Q = MakeQfromS(spikes{i}.S, binsize);
        
        QPRE{i} = zscore(full(Data(Restrict(Q, PreSleepFinal{i}))));
        QTASK{i} = zscore(full(Data(Restrict(Q, UMazeMovingEpoch{i}))));
        QPOST{i} = zscore(full(Data(Restrict(Q, PostSleepFinal{i}))));
        QCONDMOV{i} = zscore(full(Data(Restrict(Q, ConditioningMovingEpoch{i}))));
        QCONDFREEZ{i} = zscore(full(Data(Restrict(Q, ConditioningFreezingEpoch{i}))));
        QPOSTTEST{i} = zscore(full(Data(Restrict(Q, AfterConditioningMovingEpoch{i}))));
        
        clear Q
        
    end
    
end

%% Calculate overlap for each place field
for i=1:length(Dir.path)
    % Get the stats for each cell
    if isfield(spikes{i}.PlaceCells,'idx')
        if length(spikes{i}.PlaceCells.idx) > PCnum_thresh %%% Take only mice with number of place cells > PCnum_thresh
            
            S_PC = spikes{i}.S(spikes{i}.PlaceCells.idx);
            
            [overlappairs{i}, distantpairs{i}] = FindOverlappingPlaceFields(S_PC, behav{i}.CleanAlignedXtsd, behav{i}.CleanAlignedYtsd,...
                UMazeMovingEpoch{i}, overlapFactor);
            
        else
            overlappairs{i}={};
            distantpairs{i}={};
        end
    else
        overlappairs{i}={};
        distantpairs{i}={};
    end
end

%% Look at each pair overlapping
if plotexamples
    for i=1:length(Dir.path)
        
        if ~isempty(overlappairs{i})
            CPre = cell(length(overlappairs{i}),1);
            CTask = cell(length(overlappairs{i}),1);
            CPost = cell(length(overlappairs{i}),1);
            for j=1:length(overlappairs{i})
                [CPre{j},B] = CrossCorrDB(Data(Restrict(S_PC{overlappairs{i}{j}(1)},PreSleepFinal{i})),...
                    Data(Restrict(S_PC{overlappairs{i}{j}(2)},PreSleepFinal{i})),4,200);
                [CTask{j},B] = CrossCorrDB(Data(Restrict(S_PC{overlappairs{i}{j}(1)},UMazeMovingEpoch{i})),...
                    Data(Restrict(S_PC{overlappairs{i}{j}(2)},UMazeMovingEpoch{i})),4,200);
                [CPost{j},B] = CrossCorrDB(Data(Restrict(S_PC{overlappairs{i}{j}(1)},PostSleepFinal{i})),...
                    Data(Restrict(S_PC{overlappairs{i}{j}(2)},PostSleepFinal{i})),4,200);
            end
            
            % Figure
            if length(overlappairs{i}) < 15
                fff(i) = figure('units', 'normalized', 'outerposition', [0 0.55 0.7 0.4], 'Name', Dir.name{i});
                title(Dir.name{i});
                tabs = arrayfun(@(x) uitab('Title', ['Cl' num2str(spikes{i}.PlaceCells.idx(overlappairs{i}{x}(1)))...
                    ' vs Cl' num2str(spikes{i}.PlaceCells.idx(overlappairs{i}{x}(2)))]), 1:length(overlappairs{i}));
                ax = arrayfun(@(tab) axes(tab), tabs);
                arrayfun(@(k) title(ax(k), ['Cl' num2str(spikes{i}.PlaceCells.idx(overlappairs{i}{k}(1)))...
                    ' vs Cl' num2str(spikes{i}.PlaceCells.idx(overlappairs{i}{k}(2)))], 'FontSize', 16), 1:length(overlappairs{i}));
                arrayfun(@(k) set(ax(k), 'NextPlot','add', 'FontSize', 16, 'FontWeight', 'bold'), 1:length(overlappairs{i}));
                arrayfun(@(k) plot(ax(k), B, CPre{k}, 'k', 'LineWidth', 2),...
                    1:length(overlappairs{i}));
                arrayfun(@(k) plot(ax(k), B, CTask{k}, 'm', 'LineWidth', 2),...
                    1:length(overlappairs{i}));
                arrayfun(@(k) plot(ax(k), B, CPost{k}, 'b', 'LineWidth', 2),...
                    1:length(overlappairs{i}));
                arrayfun(@(k) line(ax(k),[0 0], ylim, 'Color', 'r', 'LineWidth',3), 1:length(overlappairs{i}));
                arrayfun(@(k) ylabel(ax(k), 'Count'), 1:length(overlappairs{i}));
                arrayfun(@(k) xlabel(ax(k), 'Time (ms)'), 1:length(overlappairs{i}));
                arrayfun(@(k) legend(ax(k), {'Pre', 'Task', 'Post'}), 1:length(overlappairs{i}));
            else
                for numfig = 1:ceil(length(overlappairs{i})/15)
                    try
                        rr = ((numfig-1)*15)+1:((numfig-1)*15)+15;
                        temp_o = overlappairs{i}(rr);
                    catch
                        rr = ((numfig-1)*15)+1:length(overlappairs{i});
                        temp_o = overlappairs{i}(rr);
                    end
                    fff1(numfig) = figure('units', 'normalized', 'outerposition', [0 0.55 0.7 0.4], 'Name', Dir.name{i});
                    title(Dir.name{i});
                    tabs = arrayfun(@(x) uitab('Title', ['Cl' num2str(spikes{i}.PlaceCells.idx(temp_o{x}(1)))...
                        ' vs Cl' num2str(spikes{i}.PlaceCells.idx(temp_o{x}(2)))]), 1:length(temp_o));
                    ax = arrayfun(@(tab) axes(tab), tabs);
                    arrayfun(@(k) title(ax(k), ['Cl' num2str(spikes{i}.PlaceCells.idx(temp_o{k}(1)))...
                        ' vs Cl' num2str(spikes{i}.PlaceCells.idx(temp_o{k}(2)))], 'FontSize', 16), 1:length(temp_o));
                    arrayfun(@(k) set(ax(k), 'NextPlot','add', 'FontSize', 16, 'FontWeight', 'bold'), 1:length(temp_o));
                    arrayfun(@(k) plot(ax(k+1-rr(1)), B, CPre{k}, 'k', 'LineWidth', 2),...
                        rr);
                    arrayfun(@(k) plot(ax(k+1-rr(1)), B, CTask{k}, 'm', 'LineWidth', 2),...
                        rr);
                    arrayfun(@(k) plot(ax(k+1-rr(1)), B, CPost{k}, 'b', 'LineWidth', 2),...
                        rr);
                    arrayfun(@(k) line(ax(k),[0 0], ylim, 'Color', 'r', 'LineWidth',3), 1:length(temp_o));
                    arrayfun(@(k) ylabel(ax(k), 'Count'), 1:length(temp_o));
                    arrayfun(@(k) xlabel(ax(k), 'Time (ms)'), 1:length(temp_o));
                    arrayfun(@(k) legend(ax(k), {'Pre', 'Task', 'Post'}), 1:length(temp_o));
                end
            end
        end
    end
end


%% Calculate cross-correlation between spike trains
for i=1:length(Dir.path)
    
    if do_as_in_paper
        
         % Perform calculations - overlapping cells
        if ~isempty(overlappairs{i})
            
            rho_overlap{i}.PRE = nan(1, length(overlappairs{i}));
            rho_overlap{i}.TASK = nan(1, length(overlappairs{i}));
            rho_overlap{i}.POST = nan(1, length(overlappairs{i}));
            rho_overlap{i}.CONDMOV = nan(1, length(overlappairs{i}));
            rho_overlap{i}.CONDFREEZ = nan(1, length(overlappairs{i}));
            rho_overlap{i}.POSTTEST = nan(1, length(overlappairs{i}));
        
            for j=1:length(overlappairs{i})
                pair = overlappairs{i}{j};
                % PreSleep
                rho_overlap{i}.PRE(j) = mean(xcorr(QPRE{i}(:,spikes{i}.PlaceCells.idx(pair(1))),...
                    QPRE{i}(:,spikes{i}.PlaceCells.idx(pair(2))), binsize/10, 'coeff'));
                % Task
                rho_overlap{i}.TASK(j) = mean(xcorr(QTASK{i}(:,spikes{i}.PlaceCells.idx(pair(1))),...
                    QTASK{i}(:,spikes{i}.PlaceCells.idx(pair(2))), binsize/10, 'coeff'));
                % CondMoving
                rho_overlap{i}.CONDMOV(j) = mean(xcorr(QCONDMOV{i}(:,spikes{i}.PlaceCells.idx(pair(1))),...
                    QCONDMOV{i}(:,spikes{i}.PlaceCells.idx(pair(2))), binsize/10, 'coeff'));
                % CondFreezing
                rho_overlap{i}.CONDFREEZ(j) = mean(xcorr(QCONDFREEZ{i}(:,spikes{i}.PlaceCells.idx(pair(1))),...
                    QCONDFREEZ{i}(:,spikes{i}.PlaceCells.idx(pair(2))), binsize/10, 'coeff'));
                % PostSleep
                rho_overlap{i}.POST(j) = mean(xcorr(QPOST{i}(:,spikes{i}.PlaceCells.idx(pair(1))),...
                    QPOST{i}(:,spikes{i}.PlaceCells.idx(pair(2))), binsize/10, 'coeff'));
                % PostTests
                rho_overlap{i}.POSTTEST(j) = mean(xcorr(QPOSTTEST{i}(:,spikes{i}.PlaceCells.idx(pair(1))),...
                    QPOSTTEST{i}(:,spikes{i}.PlaceCells.idx(pair(2))), binsize/10, 'coeff'));
            end
        else
            rho_overlap{i}.PRE = [];
            rho_overlap{i}.TASK = [];
            rho_overlap{i}.POST = [];
            rho_overlap{i}.CONDMOV = [];
            rho_overlap{i}.CONDFREEZ = [];
            rho_overlap{i}.POSTTEST = [];
        end
        clear pair
        
         % Perform calculations - non-overlapping cells
        if ~isempty(distantpairs{i})
            
            rho_distant{i}.PRE = nan(1, length(distantpairs{i}));
            rho_distant{i}.TASK = nan(1, length(distantpairs{i}));
            rho_distant{i}.POST = nan(1, length(distantpairs{i}));
            rho_distant{i}.CONDMOV = nan(1, length(distantpairs{i}));
            rho_distant{i}.CONDFREEZ = nan(1, length(distantpairs{i}));
            rho_distant{i}.POSTTEST = nan(1, length(distantpairs{i}));
            
            % Non-Overlapping
            for j=1:length(distantpairs{i})
                pair = distantpairs{i}{j};
                % PreSleep
                rho_distant{i}.PRE(j) = mean(xcorr(QPRE{i}(:,spikes{i}.PlaceCells.idx(pair(1))),...
                    QPRE{i}(:,spikes{i}.PlaceCells.idx(pair(2))), binsize/10, 'coeff'));
                % Task
                rho_distant{i}.TASK(j) = mean(xcorr(QTASK{i}(:,spikes{i}.PlaceCells.idx(pair(1))),...
                    QTASK{i}(:,spikes{i}.PlaceCells.idx(pair(2))), binsize/10, 'coeff'));
                % CondMoving
                rho_distant{i}.CONDMOV(j) = mean(xcorr(QCONDMOV{i}(:,spikes{i}.PlaceCells.idx(pair(1))),...
                    QCONDMOV{i}(:,spikes{i}.PlaceCells.idx(pair(2))), binsize/10, 'coeff'));
                % CondFreezing
                rho_distant{i}.CONDFREEZ(j) = mean(xcorr(QCONDFREEZ{i}(:,spikes{i}.PlaceCells.idx(pair(1))),...
                    QCONDFREEZ{i}(:,spikes{i}.PlaceCells.idx(pair(2))), binsize/10, 'coeff'));
                % PostSleep
                rho_distant{i}.POST(j) = mean(xcorr(QPOST{i}(:,spikes{i}.PlaceCells.idx(pair(1))),...
                    QPOST{i}(:,spikes{i}.PlaceCells.idx(pair(2))), binsize/10, 'coeff'));
                % PostTests
                rho_distant{i}.POSTTEST(j) = mean(xcorr(QPOSTTEST{i}(:,spikes{i}.PlaceCells.idx(pair(1))),...
                    QPOSTTEST{i}(:,spikes{i}.PlaceCells.idx(pair(2))), binsize/10, 'coeff'));
            end
        else
            rho_distant{i}.PRE = [];
            rho_distant{i}.TASK = [];
            rho_distant{i}.POST = [];
            rho_distant{i}.CONDMOV = [];
            rho_distant{i}.CONDFREEZ = [];
            rho_distant{i}.POSTTEST = [];
        end
         
    else % do_as_in_paper 
        % Perform calculations - overlapping cells
        if ~isempty(overlappairs{i})
            
            rho_overlap{i}.PRE = nan(1, length(overlappairs{i}));
            rho_overlap{i}.TASK = nan(1, length(overlappairs{i}));
            rho_overlap{i}.POST = nan(1, length(overlappairs{i}));
            rho_overlap{i}.CONDMOV = nan(1, length(overlappairs{i}));
            rho_overlap{i}.CONDFREEZ = nan(1, length(overlappairs{i}));
            rho_overlap{i}.POSTTEST = nan(1, length(overlappairs{i}));
            
            
            for j=1:length(overlappairs{i})
                pair = overlappairs{i}{j};
                % PreSleep
                temp = corrcoef(QPRE{i}(:,spikes{i}.PlaceCells.idx(pair(1))),QPRE{i}(:,spikes{i}.PlaceCells.idx(pair(2))));
                rho_overlap{i}.PRE(j) = temp(1,2); clear temp
                % Task
                temp = corrcoef(QTASK{i}(:,spikes{i}.PlaceCells.idx(pair(1))),QTASK{i}(:,spikes{i}.PlaceCells.idx(pair(2))));
                rho_overlap{i}.TASK(j) = temp(1,2); clear temp
                % CondMoving
                temp = corrcoef(QCONDMOV{i}(:,spikes{i}.PlaceCells.idx(pair(1))),QCONDMOV{i}(:,spikes{i}.PlaceCells.idx(pair(2))));
                rho_overlap{i}.CONDMOV(j) = temp(1,2); clear temp
                % CondFreezing
                temp = corrcoef(QCONDFREEZ{i}(:,spikes{i}.PlaceCells.idx(pair(1))),QCONDFREEZ{i}(:,spikes{i}.PlaceCells.idx(pair(2))));
                rho_overlap{i}.CONDFREEZ(j) = temp(1,2); clear temp
                % PostSleep
                temp = corrcoef(QPOST{i}(:,spikes{i}.PlaceCells.idx(pair(1))),QPOST{i}(:,spikes{i}.PlaceCells.idx(pair(2))));
                rho_overlap{i}.POST(j) = temp(1,2); clear temp
                % PostTests
                temp = corrcoef(QPOSTTEST{i}(:,spikes{i}.PlaceCells.idx(pair(1))),QPOSTTEST{i}(:,spikes{i}.PlaceCells.idx(pair(2))));
                rho_overlap{i}.POSTTEST(j) = temp(1,2); clear temp
            end
        else
            rho_overlap{i}.PRE = [];
            rho_overlap{i}.TASK = [];
            rho_overlap{i}.POST = [];
            rho_overlap{i}.CONDMOV = [];
            rho_overlap{i}.CONDFREEZ = [];
            rho_overlap{i}.POSTTEST = [];
        end
        clear pair
        
        % Perform calculations - non-overlapping cells
        if ~isempty(distantpairs{i})
            
            rho_distant{i}.PRE = nan(1, length(distantpairs{i}));
            rho_distant{i}.TASK = nan(1, length(distantpairs{i}));
            rho_distant{i}.POST = nan(1, length(distantpairs{i}));
            rho_distant{i}.CONDMOV = nan(1, length(distantpairs{i}));
            rho_distant{i}.CONDFREEZ = nan(1, length(distantpairs{i}));
            rho_distant{i}.POSTTEST = nan(1, length(distantpairs{i}));
            
            % Non-Overlapping
            for j=1:length(distantpairs{i})
                pair = distantpairs{i}{j};
                % PreSleep
                temp = corrcoef(QPRE{i}(:,spikes{i}.PlaceCells.idx(pair(1))),QPRE{i}(:,spikes{i}.PlaceCells.idx(pair(2))));
                rho_distant{i}.PRE(j) = temp(1,2); clear temp
                % Task
                temp = corrcoef(QTASK{i}(:,spikes{i}.PlaceCells.idx(pair(1))),QTASK{i}(:,spikes{i}.PlaceCells.idx(pair(2))));
                rho_distant{i}.TASK(j) = temp(1,2); clear temp
                % CondMoving
                temp = corrcoef(QCONDMOV{i}(:,spikes{i}.PlaceCells.idx(pair(1))),QCONDMOV{i}(:,spikes{i}.PlaceCells.idx(pair(2))));
                rho_distant{i}.CONDMOV(j) = temp(1,2); clear temp
                % CondFreezing
                temp = corrcoef(QCONDFREEZ{i}(:,spikes{i}.PlaceCells.idx(pair(1))),QCONDFREEZ{i}(:,spikes{i}.PlaceCells.idx(pair(2))));
                rho_distant{i}.CONDFREEZ(j) = temp(1,2); clear temp
                % PostSleep
                temp = corrcoef(QPOST{i}(:,spikes{i}.PlaceCells.idx(pair(1))),QPOST{i}(:,spikes{i}.PlaceCells.idx(pair(2))));
                rho_distant{i}.POST(j) = temp(1,2); clear temp
                % PostTests
                temp = corrcoef(QPOSTTEST{i}(:,spikes{i}.PlaceCells.idx(pair(1))),QPOSTTEST{i}(:,spikes{i}.PlaceCells.idx(pair(2))));
                rho_distant{i}.POSTTEST(j) = temp(1,2); clear temp
            end
        else
            rho_distant{i}.PRE = [];
            rho_distant{i}.TASK = [];
            rho_distant{i}.POST = [];
            rho_distant{i}.CONDMOV = [];
            rho_distant{i}.CONDFREEZ = [];
            rho_distant{i}.POSTTEST = [];
        end
        
    end
    
end

%% Average

% Allocate
CorrPre_mean = nan(length(Dir.path),2);
CorrPre_std = nan(length(Dir.path),2);
CorrTask_mean = nan(length(Dir.path),2);
CorrTask_std = nan(length(Dir.path),2);
CorrCondMov_mean = nan(length(Dir.path),2);
CorrCondMov_std = nan(length(Dir.path),2);
CorrCondFreez_mean = nan(length(Dir.path),2);
CorrCondFreez_std = nan(length(Dir.path),2);
CorrPost_mean = nan(length(Dir.path),2);
CorrPost_std = nan(length(Dir.path),2);
CorrPostTest_mean = nan(length(Dir.path),2);
CorrPostTest_std = nan(length(Dir.path),2);

for i = 1:length(Dir.path)
    
    CorrPre_mean(i,1) = nanmean(rho_overlap{i}.PRE); % overlapping
    CorrPre_mean(i,2) = nanmean(rho_distant{i}.PRE); % non-overlapping
    CorrPre_std(i,1) = nanstd(rho_overlap{i}.PRE); % overlapping
    CorrPre_std(i,2) = nanstd(rho_distant{i}.PRE); % non-overlapping
    
    CorrTask_mean(i,1) = nanmean(rho_overlap{i}.TASK); % overlapping
    CorrTask_mean(i,2) = nanmean(rho_distant{i}.TASK); % non-overlapping
    CorrTask_std(i,1) = nanstd(rho_overlap{i}.TASK); % overlapping
    CorrTask_std(i,2) = nanstd(rho_distant{i}.TASK); % non-overlapping
    
    CorrCondMov_mean(i,1) = nanmean(rho_overlap{i}.CONDMOV); % overlapping
    CorrCondMov_mean(i,2) = nanmean(rho_distant{i}.CONDMOV); % non-overlapping
    CorrCondMov_std(i,1) = nanstd(rho_overlap{i}.CONDMOV); % overlapping
    CorrCondMov_std(i,2) = nanstd(rho_distant{i}.CONDMOV); % non-overlapping
    
    CorrCondFreez_mean(i,1) = nanmean(rho_overlap{i}.CONDFREEZ); % overlapping
    CorrCondFreez_mean(i,2) = nanmean(rho_distant{i}.CONDFREEZ); % non-overlapping
    CorrCondFreez_std(i,1) = nanstd(rho_overlap{i}.CONDFREEZ); % overlapping
    CorrCondFreez_std(i,2) = nanstd(rho_distant{i}.CONDFREEZ); % non-overlapping
    
    CorrPost_mean(i,1) = nanmean(rho_overlap{i}.POST); % overlapping
    CorrPost_mean(i,2) = nanmean(rho_distant{i}.POST); % non-overlapping
    CorrPost_std(i,1) = nanstd(rho_overlap{i}.POST); % overlapping
    CorrPost_std(i,2) = nanstd(rho_distant{i}.POST); % non-overlapping
    
    CorrPostTest_mean(i,1) = nanmean(rho_overlap{i}.POSTTEST); % overlapping
    CorrPostTest_mean(i,2) = nanmean(rho_distant{i}.POSTTEST); % non-overlapping
    CorrPostTest_std(i,1) = nanstd(rho_overlap{i}.POSTTEST); % overlapping
    CorrPostTest_std(i,2) = nanstd(rho_distant{i}.POSTTEST); % non-overlapping
    
end

%% Pool the pairs together

CorrOPre_pooled = rho_overlap{1}.PRE';
CorrDPre_pooled = rho_distant{1}.PRE';

CorrOTask_pooled = rho_overlap{1}.TASK';
CorrDTask_pooled = rho_distant{1}.TASK';

CorrOCondMov_pooled = rho_overlap{1}.CONDMOV';
CorrDCondMov_pooled = rho_distant{1}.CONDMOV';

CorrOCondFreez_pooled = rho_overlap{1}.CONDFREEZ';
CorrDCondFreez_pooled = rho_distant{1}.CONDFREEZ';

CorrOPost_pooled = rho_overlap{1}.POST';
CorrDPost_pooled = rho_distant{1}.POST';

CorrOPostTest_pooled = rho_overlap{1}.POSTTEST';
CorrDPostTest_pooled = rho_distant{1}.POSTTEST';

if length(Dir.path) > 1
    for i = 2:length(Dir.path)
        CorrOPre_pooled = [CorrOPre_pooled; rho_overlap{i}.PRE'];
        CorrDPre_pooled = [CorrDPre_pooled; rho_distant{i}.PRE'];
        
        CorrOTask_pooled = [CorrOTask_pooled; rho_overlap{i}.TASK'];
        CorrDTask_pooled = [CorrDTask_pooled; rho_distant{i}.TASK'];
        
        CorrOCondMov_pooled = [CorrOCondMov_pooled; rho_overlap{i}.CONDMOV'];
        CorrDCondMov_pooled = [CorrDCondMov_pooled; rho_distant{i}.CONDMOV'];
        
        CorrOCondFreez_pooled = [CorrOCondFreez_pooled; rho_overlap{i}.CONDFREEZ'];
        CorrDCondFreez_pooled = [CorrDCondFreez_pooled; rho_distant{i}.CONDFREEZ'];
        
        CorrOPost_pooled = [CorrOPost_pooled; rho_overlap{i}.POST'];
        CorrDPost_pooled = [CorrDPost_pooled; rho_distant{i}.POST'];
        
        CorrOPostTest_pooled = [CorrOPostTest_pooled; rho_overlap{i}.POSTTEST'];
        CorrDPostTest_pooled = [CorrDPostTest_pooled; rho_distant{i}.POSTTEST'];
        
    end
end

%% Plot averaged over mice

% Full figure
fh = figure('units', 'normalized', 'outerposition', [0 0 0.7 0.7]);
if plotmice
    % Average over mice
    pl_pre_std = nanstd(CorrPre_mean,0,1);
    pl_pre_mean = nanmean(CorrPre_mean,1);
    
    pl_task_std = nanstd(CorrTask_mean,0,1);
    pl_task_mean = nanmean(CorrTask_mean,1);
    
    pl_condmov_std = nanstd(CorrCondMov_mean,0,1);
    pl_condmov_mean = nanmean(CorrCondMov_mean,1);
    
    pl_condfreez_std = nanstd(CorrCondFreez_mean,0,1);
    pl_condfreez_mean = nanmean(CorrCondFreez_mean,1);
    
    pl_post_std = nanstd(CorrPost_mean,0,1);
    pl_post_mean = nanmean(CorrPost_mean,1);
    
    pl_posttest_std = nanstd(CorrPostTest_mean,0,1);
    pl_posttest_mean = nanmean(CorrPostTest_mean,1);
    
    b = barwitherr([pl_pre_std; pl_task_std; pl_condmov_std; pl_condfreez_std; pl_post_std; pl_posttest_std],...
        [pl_pre_mean; pl_task_mean; pl_condmov_mean; pl_condfreez_mean; pl_post_mean; pl_posttest_mean]);
else
    % Average over all pooled neurons
%     pl_pre_std = [nanstd(CorrOPre_pooled) nanstd(CorrDPre_pooled)];
    pl_pre_sem = [nansem_db(CorrOPre_pooled) nansem_db(CorrDPre_pooled)];
    pl_pre_mean = [nanmean(CorrOPre_pooled) nanmean(CorrDPre_pooled)];
    
%     pl_task_std = [nanstd(CorrOTask_pooled) nanstd(CorrDTask_pooled)];
    pl_task_sem = [nansem_db(CorrOTask_pooled) nansem_db(CorrDTask_pooled)];
    pl_task_mean = [nanmean(CorrOTask_pooled) nanmean(CorrDTask_pooled)];
    
%     pl_condmov_std = [nanstd(CorrOCondMov_pooled) nanstd(CorrDCondMov_pooled)];
    pl_condmov_sem = [nansem_db(CorrOCondMov_pooled) nansem_db(CorrDCondMov_pooled)];
    pl_condmov_mean = [nanmean(CorrOCondMov_pooled) nanmean(CorrDCondMov_pooled)];
    
%     pl_condfreez_std = [nanstd(CorrOCondFreez_pooled) nanstd(CorrDCondFreez_pooled)];
    pl_condfreez_sem = [nansem_db(CorrOCondFreez_pooled) nansem_db(CorrDCondFreez_pooled)];
    pl_condfreez_mean = [nanmean(CorrOCondFreez_pooled) nanmean(CorrDCondFreez_pooled)];
    
%     pl_post_std = [nanstd(CorrOPost_pooled) nanstd(CorrDPost_pooled)];
    pl_post_sem = [nansem_db(CorrOPost_pooled) nansem_db(CorrDPost_pooled)];
    pl_post_mean = [nanmean(CorrOPost_pooled) nanmean(CorrDPost_pooled)];
    
%     pl_posttest_std = [nanstd(CorrOPostTest_pooled) nanstd(CorrDPostTest_pooled)];
    pl_posttest_sem = [nansem_db(CorrOPostTest_pooled) nansem_db(CorrDPostTest_pooled)];
    pl_posttest_mean = [nanmean(CorrOPostTest_pooled) nanmean(CorrDPostTest_pooled)];
    
%     b = barwitherr([pl_pre_std; pl_task_std; pl_condmov_std; pl_condfreez_std; pl_post_std; pl_posttest_std],...
%         [pl_pre_mean; pl_task_mean; pl_condmov_mean; pl_condfreez_mean; pl_post_mean; pl_posttest_mean]);
    b = barwitherr([pl_pre_sem; pl_task_sem; pl_condmov_sem; pl_condfreez_sem; pl_post_sem; pl_posttest_sem],...
        [pl_pre_mean; pl_task_mean; pl_condmov_mean; pl_condfreez_mean; pl_post_mean; pl_posttest_mean]);
end
set(gca, 'FontSize', 14, 'FontWeight',  'bold');
set(gca, 'LineWidth', 3);
b(1).BarWidth = 0.8;
b(1).FaceColor = 'k';
b(2).FaceColor = 'w';
b(1).LineWidth = 3;
b(2).LineWidth = 3;
x = [b(1).XData + [b(1).XOffset]; b(1).XData - [b(1).XOffset]];
hold on
set(gca,'Xtick',[1:6],'XtickLabel',{'PreSleep', 'PreExplorations',...
    'Cond running', 'Cond Freezing', 'PostSleep', 'PostTests'})
ylabel('Corellation')
if plotmice
    title('Averaged over mice');
else
    title('All neurons pooled');
end
hold off
box off
set(gca, 'FontSize', 16, 'FontWeight',  'bold');
set(gca, 'LineWidth', 3);
% title
lg = legend('Overlapping PCs', 'Non-overlapping PCs');
lg.FontSize = 14;

% Short figure
f1 = figure('units', 'normalized', 'outerposition', [0 0 0.4 0.7]);
if plotmice
    % Average over mice
    pl_pre_std = nanstd(CorrPre_mean,0,1);
    pl_pre_mean = nanmean(CorrPre_mean,1);
    
    pl_task_std = nanstd(CorrTask_mean,0,1);
    pl_task_mean = nanmean(CorrTask_mean,1); 
    
    pl_post_std = nanstd(CorrPost_mean,0,1);
    pl_post_mean = nanmean(CorrPost_mean,1);
    
    b = barwitherr([pl_pre_std; pl_task_std; pl_post_std],...
        [pl_pre_mean; pl_task_mean; pl_post_mean]);
else
    % Average over all pooled neurons
    pl_pre_std = [nanstd(CorrOPre_pooled) nanstd(CorrDPre_pooled)];
%     pl_pre_sem = [nansem_db(CorrOPre_pooled) nansem_db(CorrDPre_pooled)];
    pl_pre_mean = [nanmean(CorrOPre_pooled) nanmean(CorrDPre_pooled)];
    
    pl_task_std = [nanstd(CorrOTask_pooled) nanstd(CorrDTask_pooled)];
%     pl_task_sem = [nansem_db(CorrOTask_pooled) nansem_db(CorrDTask_pooled)];
    pl_task_mean = [nanmean(CorrOTask_pooled) nanmean(CorrDTask_pooled)];
    
    pl_post_std = [nanstd(CorrOPost_pooled) nanstd(CorrDPost_pooled)];
%     pl_post_sem = [nansem_db(CorrOPost_pooled) nansem_db(CorrDPost_pooled)];
    pl_post_mean = [nanmean(CorrOPost_pooled) nanmean(CorrDPost_pooled)]; 
    
    b = barwitherr([pl_pre_std; pl_task_std; pl_post_std],...
        [pl_pre_mean; pl_task_mean; pl_post_mean]);
%     b = barwitherr([pl_pre_sem; pl_task_sem; pl_post_sem],...
%         [pl_pre_mean; pl_task_mean; pl_post_mean]);
end
set(gca, 'FontSize', 18, 'FontWeight',  'bold');
set(gca, 'LineWidth', 3);
b(1).BarWidth = 0.8;
b(1).FaceColor = 'k';
b(2).FaceColor = 'w';
b(1).LineWidth = 3;
b(2).LineWidth = 3;
x = [b(1).XData + [b(1).XOffset]; b(1).XData - [b(1).XOffset]];
hold on
set(gca,'Xtick',[1:6],'XtickLabel',{'PreSleep', 'NoStimTask', 'PostSleep'})
ylabel('Corellation')
if plotmice
    title('Averaged over mice');
else
    title('All neurons pooled');
end
hold off
box off
set(gca, 'FontSize', 18, 'FontWeight',  'bold');
set(gca, 'LineWidth', 3);
% title
lg = legend('Overlapping PCs', 'Non-overlapping PCs');
lg.FontSize = 14;
% Stats - TODO
% [p,h5,stats] = signrank([CorrOPre_pooled; CorrDPre_pooled]);
% if p < 0.05
%     sigstar_DB({{1,2}},p,0, 'StarSize',14);
% end
% [p,h5,stats] = signrank([CorrOTask_pooled; CorrDTask_pooled]);
% if p < 0.05
%     sigstar_DB({{3,4}},p,0, 'StarSize',14);
% end
% [p,h5,stats] = signrank([CorrOPost_pooled; CorrDPost_pooled]);
% if p < 0.05
%     sigstar_DB({{5,6}},p,0, 'StarSize',14);
% end

% Save the figure
if savfig
    saveas(f1,[dropbox '/MOBS_workingON/Dima/Ongoing_results/React-Replay/Pairwise_Small.fig']);
    saveFigure(f1,'Pairwise_Small',...
        [dropbox '/MOBS_workingON/Dima/Ongoing_results/React-Replay/']);
end
