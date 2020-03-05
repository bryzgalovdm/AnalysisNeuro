%%%% PairwiseCorrelationsAcrossMice_WM1994_DB
% 
% 
% This script performs pairwise correlations on all overlapping
% and non-overlapping place cells of each particular mouse,
% averages them out and then does global mean across all mice
% 
% Please go inside the script and check the parameters


%% Parameters
% Mice that go in the analysis
nmouse = [797 798 828 861 882 905 906 911 912 977 994];
% nmouse = [906 912]; % Had PreMazes
% nmouse = [905 911]; % Did not have PreMazes

% Get paths of each individual mouse
Dir = PathForExperimentsERC_Dima('UMazePAG');
Dir = RestrictPathForExperiment(Dir,'nMice',nmouse);

% Parameters of cross-correlograms
binsize = 5; % How much ms do you wnat a bin to contain?
CCtime = 100; % How much time around a 0 you want to have (+-, in ms)?

% How many pixels should overlap to define an overlap in the place field?
overlapFactor = 50;

% Do you want to save the figures?
sav = 0;

% Do you want to restrict your sleep time? If not put []
SleepTimeToRestrict = 20*60*1e4;

%% Allocate memory
CPRE_O_Mean = zeros(1,length(Dir.path));
CPRE_O_std = zeros(1,length(Dir.path));
CPRE_D_Mean = zeros(1,length(Dir.path));
CPRE_D_std = zeros(1,length(Dir.path));

CTASK_O_Mean = zeros(1,length(Dir.path));
CTASK_O_std = zeros(1,length(Dir.path));
CTASK_D_Mean = zeros(1,length(Dir.path));
CTASK_D_std = zeros(1,length(Dir.path));

CPOST_O_Mean = zeros(1,length(Dir.path));
CPOST_O_std = zeros(1,length(Dir.path));
CPOST_D_Mean = zeros(1,length(Dir.path));
CPOST_D_std = zeros(1,length(Dir.path));

CCONDMOV_O_Mean = zeros(1,length(Dir.path));
CCONDMOV_O_std = zeros(1,length(Dir.path));
CCONDMOV_D_Mean = zeros(1,length(Dir.path));
CCONDMOV_D_std = zeros(1,length(Dir.path));

CCONDFREEZ_O_Mean = zeros(1,length(Dir.path));
CCONDFREEZ_O_std = zeros(1,length(Dir.path));
CCONDFREEZ_D_Mean = zeros(1,length(Dir.path));
CCONDFREEZ_D_std = zeros(1,length(Dir.path));

CPOSTTEST_O_Mean = zeros(1,length(Dir.path));
CPOSTTEST_O_std = zeros(1,length(Dir.path));
CPOSTTEST_D_Mean = zeros(1,length(Dir.path));
CPOSTTEST_D_std = zeros(1,length(Dir.path));

%% Load Data
for j=1:length(Dir.path)
    
    cd(Dir.path{j}{1});
    load('SpikeData.mat','S','PlaceCells');
    load('behavResources.mat','SessionEpoch','CleanVtsd','CleanAlignedXtsd','CleanAlignedYtsd','FreezeAccEpoch');
    if strcmp(Dir.name{j}, 'Mouse906') || strcmp(Dir.name{j}, 'Mouse977') % Mice with bad OB-based sleep scoring
        load('SleepScoring_Accelero.mat','SWSEpoch','REMEpoch','Sleep'); % REM and Sleep are not used
    else
        load('SleepScoring_OBGamma.mat','SWSEpoch','REMEpoch','Sleep');  % REM and Sleep are not used
    end
    
    
    
    %% Epochs
    
    % BaselineExplo Epoch
    UMazeEpoch{j} = or(SessionEpoch.Hab,SessionEpoch.TestPre1);
    UMazeEpoch{j} = or(UMazeEpoch{j},SessionEpoch.TestPre2);
    UMazeEpoch{j} = or(UMazeEpoch{j},SessionEpoch.TestPre3);
    UMazeEpoch{j} = or(UMazeEpoch{j},SessionEpoch.TestPre4);
    
    % Conditioning
    ConditioningEpoch{j} = or(SessionEpoch.Cond1,SessionEpoch.Cond2);
    ConditioningEpoch{j} = or(ConditioningEpoch{j},SessionEpoch.Cond3);
    ConditioningEpoch{j} = or(ConditioningEpoch{j},SessionEpoch.Cond4);
    
    
    % After Conditioning
    AfterConditioningEpoch{j} = or(SessionEpoch.TestPost1,SessionEpoch.TestPost2);
    AfterConditioningEpoch{j} = or(AfterConditioningEpoch{j},SessionEpoch.TestPost3);
    AfterConditioningEpoch{j} = or(AfterConditioningEpoch{j},SessionEpoch.TestPost4);
    
    % Locomotion threshold
    VtsdSmoothed  = tsd(Range(CleanVtsd),movmedian(Data(CleanVtsd),5));
    LocomotionEpoch{j} = thresholdIntervals(VtsdSmoothed,3,'Direction','Above');
    
    % Get resulting epoch
    UMazeMovingEpoch{j} = and(LocomotionEpoch{j}, UMazeEpoch{j});
    AfterConditioningMovingEpoch{j} = and(LocomotionEpoch{j}, AfterConditioningEpoch{j});
    ConditioningMovingEpoch{j} = and(LocomotionEpoch{j}, ConditioningEpoch{j});
    ConditioningFreezingEpoch{j} = and(FreezeAccEpoch, ConditioningEpoch{j});
    
    % SleepEpochs
    if isempty(SleepTimeToRestrict)
        PreSleepFinal = and(SessionEpoch.PreSleep, SWSEpoch);
        PostSleepFinal = and(SessionEpoch.PostSleep, SWSEpoch);
    else
        PreSleepFinal = RestrictToTime(and(SessionEpoch.PreSleep, SWSEpoch),SleepTimeToRestrict);
        PostSleepFinal = RestrictToTime(and(SessionEpoch.PostSleep, SWSEpoch),SleepTimeToRestrict);
    end
        
    
    %% Calculate overlap for each place field
    
    % Get the stats for each cell
    if isfield(PlaceCells,'idx')
        if length(PlaceCells.idx)>2 %%% Take only mice with number of place cells > 2
            for i=1:length(PlaceCells.idx)
                try % Get the place field
                    [map{i},mapS,stats{i},px,py,FR(i),sizeFinal,PrField{i},C,ScField,pfH,pf] =...
                        PlaceField_DB(Restrict(S{PlaceCells.idx(i)},UMazeMovingEpoch{j}),...
                        Restrict(CleanAlignedXtsd, UMazeMovingEpoch{j}),...
                        Restrict(CleanAlignedYtsd, UMazeMovingEpoch{j}),'threshold',0.5, 'plotresults',0);
                    close all;
                catch
                    map{i} = [];
                    stats{i}=[];
                    PrField{i}=[];
                end
            end
            
            % FindOverlappingPlaceCells
            o = 1;
            n = 1;
            for i=1:length(PlaceCells.idx)
                for k = i+1:length(PlaceCells.idx)
                    cell1 = PlaceCells.idx(i);
                    cell2 = PlaceCells.idx(k);
                    if iscell(stats{i}.field) && ~iscell(stats{k}.field) % If cell1 has two fields, and cell2 has one
                        OverlappedFields{1} = stats{i}.field{1} & stats{k}.field;
                        OverlappedFields{2} = stats{i}.field{2} & stats{k}.field;
                        numOverlap{1}(i,k) = nnz(OverlappedFields{1});
                        numOverlap{2}(i,k) = nnz(OverlappedFields{2});
                        % If either one ot the other place field overlaps, consider overlap
                        if numOverlap{1}(i,k) > overlapFactor || numOverlap{2}(i,k) > overlapFactor
                            overlappairs{j}{o} = [cell1, cell2];
                            o=o+1;
                        else
                            distantpairs{j}{n} = [cell1, cell2];
                            n=n+1;
                        end
                    elseif ~iscell(stats{i}.field) && iscell(stats{k}.field) % If cell1 has one field, and cell2 has two
                        OverlappedFields{1} = stats{i}.field & stats{k}.field{1};
                        OverlappedFields{2} = stats{i}.field & stats{k}.field{2};
                        numOverlap{1}(i,k) = nnz(OverlappedFields{1});
                        numOverlap{2}(i,k) = nnz(OverlappedFields{2});
                        if numOverlap{1}(i,k) > overlapFactor || numOverlap{2}(i,k) > overlapFactor
                            overlappairs{j}{o} = [cell1, cell2];
                            o=o+1;
                        else
                            distantpairs{j}{n} = [cell1, cell2];
                            n=n+1;
                        end
                    elseif iscell(stats{i}.field) && iscell(stats{i}.field) % If both cell1 and cell2 have two fields
                        OverlappedFields{1} = stats{i}.field{1} & stats{k}.field{1};
                        OverlappedFields{2} = stats{i}.field{1} & stats{k}.field{2};
                        OverlappedFields{3} = stats{i}.field{2} & stats{k}.field{1};
                        OverlappedFields{4} = stats{i}.field{2} & stats{k}.field{2};
                        numOverlap{1}(i,k) = nnz(OverlappedFields{1});
                        numOverlap{2}(i,k) = nnz(OverlappedFields{2});
                        numOverlap{3}(i,k) = nnz(OverlappedFields{3});
                        numOverlap{4}(i,k) = nnz(OverlappedFields{4});
                        if numOverlap{1}(i,k) > overlapFactor || numOverlap{2}(i,k) > overlapFactor ||...
                                numOverlap{3}(i,k) > overlapFactor || numOverlap{4}(i,k) > overlapFactor
                            overlappairs{j}{o} = [cell1, cell2];
                            o=o+1;
                        else
                            distantpairs{j}{n} = [cell1, cell2];
                            n=n+1;
                        end
                    else
                        OverlappedFields = stats{i}.field & stats{k}.field; % If both cell1 and cell2 have one field
                        numOverlap(i,j) = nnz(OverlappedFields);
                        if numOverlap(i,j) > overlapFactor
                            overlappairs{j}{o} = [cell1, cell2];
                            o=o+1;
                        else
                            distantpairs{j}{n} = [cell1, cell2];
                            n=n+1;
                        end
                    end
                    clear OverlappedFields numOverlap
                end
            end
            if ~exist('overlappairs','var') || numel(overlappairs) < j
                overlappairs{j} = [];
            end
            if ~exist('distantpairs','var') || numel(distantpairs) < j
                distantpairs{j} = [];
            end
        else
            overlappairs{j}=[];
            distantpairs{j}=[];
        end
    else
        overlappairs{j}=[];
        distantpairs{j}=[];
    end

    %% Calculate cross-correlation
    
    % Allocate memory - overlapping cells
    C_PreSleep_O{j} = zeros(1,length(overlappairs{j}));
    C_Task_O{j} = zeros(1,length(overlappairs{j}));
    C_CondMov_O{j} = zeros(1,length(overlappairs{j}));
    C_CondFreez_O{j} = zeros(1,length(overlappairs{j}));
    C_PostSleep_O{j} = zeros(1,length(overlappairs{j}));
    C_PostTest_O{j} = zeros(1,length(overlappairs{j}));
    
    % Perform calculations - overlapping cells
    if ~isempty(overlappairs{j})
        % Overlapped
        for i=1:length(overlappairs{j})
            pair = overlappairs{j}{i};
            % PreSleep
            clear C_PreSleep C_Task C_PostSleep C_CondMov C_CondFreez C_PostTest
            [C_PreSleep,B]=CrossCorrDB(Data(Restrict(S{pair(1)},PreSleepFinal)),...
                Data(Restrict(S{pair(2)},PreSleepFinal)),binsize,CCtime/binsize*2);
            C_PreSleep_O{j}(i) = mean(C_PreSleep);
            if C_PreSleep_O{j}(i) == 0
                C_PreSleep_O{j}(i) = NaN;
            end
            % Task
            [C_Task,B]=CrossCorrDB(Data(Restrict(S{pair(1)},UMazeMovingEpoch{j})),...
                Data(Restrict(S{pair(2)},UMazeMovingEpoch{j})),binsize,CCtime/binsize*2);
            C_Task_O{j}(i) = mean(C_Task);
            if C_Task_O{j}(i) == 0
                C_Task_O{j}(i) = NaN;
            end
            % CondMoving
            [C_CondMov,B]=CrossCorrDB(Data(Restrict(S{pair(1)},ConditioningMovingEpoch{j})),...
                Data(Restrict(S{pair(2)},ConditioningMovingEpoch{j})),binsize,CCtime/binsize*2);
            C_CondMov_O{j}(i) = mean(C_CondMov);
            if C_CondMov_O{j}(i) == 0
                C_CondMov_O{j}(i) = NaN;
            end
            % CondFreezing
            [C_CondFreez,B]=CrossCorrDB(Data(Restrict(S{pair(1)},ConditioningFreezingEpoch{j})),...
                Data(Restrict(S{pair(2)},ConditioningFreezingEpoch{j})),binsize,CCtime/binsize*2);
            C_CondFreez_O{j}(i) = mean(C_CondFreez);
            if C_CondFreez_O{j}(i) == 0
                C_CondFreez_O{j}(i) = NaN;
            end
            % PostSleep
            [C_PostSleep,B]=CrossCorrDB(Data(Restrict(S{pair(1)},PostSleepFinal)),...
                Data(Restrict(S{pair(2)},PostSleepFinal)),binsize,CCtime/binsize*2);
            C_PostSleep_O{j}(i) = mean(C_PostSleep);
            if C_PostSleep_O{j}(i) == 0
                C_PostSleep_O{j}(i) = NaN;
            end
            % PostTests
            [C_PostTest,B]=CrossCorrDB(Data(Restrict(S{pair(1)},AfterConditioningMovingEpoch{j})),...
                Data(Restrict(S{pair(2)},AfterConditioningMovingEpoch{j})),binsize,CCtime/binsize*2);
            C_PostTest_O{j}(i) = mean(C_PostTest);
            if C_PostTest_O{j}(i) == 0
                C_PostTest_O{j}(i) = NaN;
            end
        end
    else
        C_PreSleep_O{j} = [];
        C_Task_O{j} = [];
        C_CondMov_O{j} = [];
        C_CondFreez_O{j} = [];
        C_PostSleep_O{j} = [];
        C_PostTest_O{j} = [];
    end
    
    % Allocate memory - non-overlapping cells
    C_PreSleep_D{j} = zeros(1,length(distantpairs{j}));
    C_Task_D{j} = zeros(1,length(distantpairs{j}));
    C_CondMov_D{j} = zeros(1,length(distantpairs{j}));
    C_CondFreez_D{j} = zeros(1,length(distantpairs{j}));
    C_PostSleep_D{j} = zeros(1,length(distantpairs{j}));
    C_PostTest_D{j} = zeros(1,length(distantpairs{j}));
    
    % Perform calculations - non-overlapping cells
    if ~isempty(distantpairs{j})
        % Non-Overlapped
        for i=1:length(distantpairs{j})
            pair = distantpairs{j}{i};
            % PreSleep
            clear C_PreSleep C_Task C_PostSleep C_CondMov C_CondFreez C_PostTest
            [C_PreSleep,B]=CrossCorrDB(Data(Restrict(S{pair(1)},PreSleepFinal)),...
                Data(Restrict(S{pair(2)},PreSleepFinal)),binsize,CCtime/binsize*2);
            C_PreSleep_D{j}(i) = mean(C_PreSleep);
            if C_PreSleep_D{j}(i) == 0
                C_PreSleep_D{j}(i) = NaN;
            end
            % Task
            [C_Task,B]=CrossCorrDB(Data(Restrict(S{pair(1)},UMazeMovingEpoch{j})),...
                Data(Restrict(S{pair(2)},UMazeMovingEpoch{j})),binsize,CCtime/binsize*2);
            C_Task_D{j}(i) = mean(C_Task);
            if C_Task_D{j}(i) == 0
                C_Task_D{j}(i) = NaN;
            end
            % CondMoving
            [C_CondMov,B]=CrossCorrDB(Data(Restrict(S{pair(1)},ConditioningMovingEpoch{j})),...
                Data(Restrict(S{pair(2)},ConditioningMovingEpoch{j})),binsize,CCtime/binsize*2);
            C_CondMov_D{j}(i) = mean(C_CondMov);
            if C_CondMov_D{j}(i) == 0
                C_CondMov_D{j}(i) = NaN;
            end
            % CondFreezing
            [C_CondFreez,B]=CrossCorrDB(Data(Restrict(S{pair(1)},ConditioningFreezingEpoch{j})),...
                Data(Restrict(S{pair(2)},ConditioningFreezingEpoch{j})),binsize,CCtime/binsize*2);
            C_CondFreez_D{j}(i) = mean(C_CondFreez);
            if C_CondFreez_D{j}(i) == 0
                C_CondFreez_D{j}(i) = NaN;
            end
            % PostSleep
            [C_PostSleep,B]=CrossCorrDB(Data(Restrict(S{pair(1)},PostSleepFinal)),...
                Data(Restrict(S{pair(2)},PostSleepFinal)),binsize,CCtime/binsize*2);
            C_PostSleep_D{j}(i) = mean(C_PostSleep);
            if C_PostSleep_D{j}(i) == 0
                C_PostSleep_D{j}(i) = NaN;
            end
            % PostTests
            [C_PostTest,B]=CrossCorrDB(Data(Restrict(S{pair(1)},AfterConditioningMovingEpoch{j})),...
                Data(Restrict(S{pair(2)},AfterConditioningMovingEpoch{j})),binsize,CCtime/binsize*2);
            C_PostTest_D{j}(i) = mean(C_PostTest);
            if C_PostTest_D{j}(i) == 0
                C_PostTest_D{j}(i) = NaN;
            end
        end
    else
        C_PreSleep_D{j} = [];
        C_Task_D{j} = [];
        C_CondMov_D{j} = [];
        C_CondFreez_D{j} = [];
        C_PostSleep_D{j} = [];
        C_PostTest_D{j} = [];
    end

    
    
    %% Average
   
        CPRE_O_Mean(j) = nanmean(C_PreSleep_O{j});
        CPRE_O_std(j) = nanstd(C_PreSleep_O{j});
        CPRE_D_Mean(j) = nanmean(C_PreSleep_D{j});
        CPRE_D_std(j) = nanstd(C_PreSleep_D{j});
    
        CTASK_O_Mean(j) = nanmean(C_Task_O{j});
        CTASK_O_std(j) = nanstd(C_Task_O{j});
        CTASK_D_Mean(j) = nanmean(C_Task_D{j});
        CTASK_D_std(j) = nanstd(C_Task_D{j});
        
        CCONDMOV_O_Mean(j) = nanmean(C_CondMov_O{j});
        CCONDMOV_O_std(j) = nanstd(C_CondMov_O{j});
        CCONDMOV_D_Mean(j) = nanmean(C_CondMov_D{j});
        CCONDMOV_D_std(j) = nanstd(C_CondMov_D{j});
        
        CCONDFREEZ_O_Mean(j) = nanmean(C_CondFreez_O{j});
        CCONDFREEZ_O_std(j) = nanstd(C_CondFreez_O{j});
        CCONDFREEZ_D_Mean(j) = nanmean(C_CondFreez_D{j});
        CCONDFREEZ_D_std(j) = nanstd(C_CondFreez_D{j});

        CPOST_O_Mean(j) = nanmean(C_PostSleep_O{j});
        CPOST_O_std(j) = nanstd(C_PostSleep_O{j});
        CPOST_D_Mean(j) = nanmean(C_PostSleep_D{j});
        CPOST_D_std(j) = nanstd(C_PostSleep_D{j});

        CPOSTTEST_O_Mean(j) = nanmean(C_PostTest_O{j});
        CPOSTTEST_O_std(j) = nanstd(C_PostTest_O{j});
        CPOSTTEST_D_Mean(j) = nanmean(C_PostTest_D{j});
        CPOSTTEST_D_std(j) = nanstd(C_PostTest_D{j});
    
    
   clear CleanVtsd S SessionEpoch PlaceCells Sleep SWSEpoch REMEpoch map stats PrField CleanAlignedXtsd
   clear CleanAlignedYtsd PreSleepFinal PostSleepFinal
end

%% Pool the Neurons

C_PreSleep_O_Pooled = cell2mat(C_PreSleep_O);
C_Task_O_Pooled = cell2mat(C_Task_O);
C_CondMov_O_Pooled = cell2mat(C_CondMov_O);
C_CondFreez_O_Pooled = cell2mat(C_CondFreez_O);
C_PostSleep_O_Pooled = cell2mat(C_PostSleep_O);
C_PostTest_O_Pooled = cell2mat(C_PostTest_O);

C_PreSleep_D_Pooled = cell2mat(C_PreSleep_D);
C_Task_D_Pooled = cell2mat(C_Task_D);
C_CondMov_D_Pooled = cell2mat(C_CondMov_D);
C_CondFreez_D_Pooled = cell2mat(C_CondFreez_D);
C_PostSleep_D_Pooled = cell2mat(C_PostSleep_D);
C_PostTest_D_Pooled = cell2mat(C_PostTest_D);

%% Inter-subject average
CPRE_O_Mean_Av = nanmean(CPRE_O_Mean);
CPRE_O_std_Av = nanstd(CPRE_O_Mean);
CPRE_D_Mean_Av = nanmean(CPRE_D_Mean);
CPRE_D_std_Av = nanstd(CPRE_D_Mean);
CPRE_OD_Mean_Av = [CPRE_O_Mean_Av CPRE_D_Mean_Av];
CPRE_OD_std_Av = [CPRE_O_std_Av CPRE_D_std_Av];

CTASK_O_Mean_Av = nanmean(CTASK_O_Mean);
CTASK_O_std_Av = nanstd(CTASK_O_Mean);
CTASK_D_Mean_Av = nanmean(CTASK_D_Mean);
CTASK_D_std_Av = nanstd(CTASK_D_Mean);
CTASK_OD_Mean_Av = [CTASK_O_Mean_Av CTASK_D_Mean_Av];
CTASK_OD_std_Av = [CTASK_O_std_Av CTASK_D_std_Av];

CCONDMOV_O_Mean_Av = nanmean(CCONDMOV_O_Mean);
CCONDMOV_O_std_Av = nanstd(CCONDMOV_O_Mean);
CCONDMOV_D_Mean_Av = nanmean(CCONDMOV_D_Mean);
CCONDMOV_D_std_Av = nanstd(CCONDMOV_D_Mean);
CCONDMOV_OD_Mean_Av = [CCONDMOV_O_Mean_Av CCONDMOV_D_Mean_Av];
CCONDMOV_OD_std_Av = [CCONDMOV_O_std_Av CCONDMOV_D_std_Av];

CCONDFREEZ_O_Mean_Av = nanmean(CCONDFREEZ_O_Mean);
CCONDFREEZ_O_std_Av = nanstd(CCONDFREEZ_O_Mean);
CCONDFREEZ_D_Mean_Av = nanmean(CCONDFREEZ_D_Mean);
CCONDFREEZ_D_std_Av = nanstd(CCONDFREEZ_D_Mean);
CCONDFREEZ_OD_Mean_Av = [CCONDFREEZ_O_Mean_Av CCONDFREEZ_D_Mean_Av];
CCONDFREEZ_OD_std_Av = [CCONDFREEZ_O_std_Av CCONDFREEZ_D_std_Av];

CPOST_O_Mean_Av = nanmean(CPOST_O_Mean);
CPOST_O_std_Av = nanstd(CPOST_O_Mean);
CPOST_D_Mean_Av = nanmean(CPOST_D_Mean);
CPOST_D_std_Av = nanstd(CPOST_D_Mean);
CPOST_OD_Mean_Av = [CPOST_O_Mean_Av CPOST_D_Mean_Av];
CPOST_OD_std_Av = [CPOST_O_std_Av CPOST_D_std_Av];

CPOSTTEST_O_Mean_Av = nanmean(CPOSTTEST_O_Mean);
CPOSTTEST_O_std_Av = nanstd(CPOSTTEST_O_Mean);
CPOSTTEST_D_Mean_Av = nanmean(CPOSTTEST_D_Mean);
CPOSTTEST_D_std_Av = nanstd(CPOSTTEST_D_Mean);
CPOSTTEST_OD_Mean_Av = [CPOSTTEST_O_Mean_Av CPOSTTEST_D_Mean_Av];
CPOSTTEST_OD_std_Av = [CPOSTTEST_O_std_Av CPOSTTEST_D_std_Av];
    
%%%% Stats
PreStats = signrank(CPRE_O_Mean, CPRE_D_Mean);
TaskStats = signrank(CTASK_O_Mean, CTASK_D_Mean);
CondMovStats = signrank(CCONDMOV_O_Mean, CCONDMOV_D_Mean);
CondFreezStats = signrank(CCONDFREEZ_O_Mean, CCONDFREEZ_D_Mean);
PostStats = signrank(CPOST_O_Mean, CPOST_D_Mean);
PosTestStats = signrank(CPOSTTEST_O_Mean, CPOSTTEST_D_Mean);



%% Inter-neuron average
CPRE_O_Mean_Pooled = nanmean(C_PreSleep_O_Pooled);
CPRE_O_std_Pooled = nanstd(C_PreSleep_O_Pooled);
CPRE_D_Mean_Pooled = nanmean(C_PreSleep_D_Pooled);
CPRE_D_std_Pooled = nanstd(C_PreSleep_D_Pooled);
CPRE_OD_Mean_Pooled = [CPRE_O_Mean_Pooled CPRE_D_Mean_Pooled];
CPRE_OD_std_Pooled = [CPRE_O_std_Pooled CPRE_D_std_Pooled];

CTASK_O_Mean_Pooled = nanmean(C_Task_O_Pooled);
CTASK_O_std_Pooled = nanstd(C_Task_O_Pooled);
CTASK_D_Mean_Pooled = nanmean(C_Task_D_Pooled);
CTASK_D_std_Pooled = nanstd(C_Task_D_Pooled);
CTASK_OD_Mean_Pooled = [CTASK_O_Mean_Pooled CTASK_D_Mean_Pooled];
CTASK_OD_std_Pooled = [CTASK_O_std_Pooled CTASK_D_std_Pooled];

CCONDMOV_O_Mean_Pooled = nanmean(C_CondMov_O_Pooled);
CCONDMOV_O_std_Pooled = nanstd(C_CondMov_O_Pooled);
CCONDMOV_D_Mean_Pooled = nanmean(C_CondMov_D_Pooled);
CCONDMOV_D_std_Pooled = nanstd(C_CondMov_D_Pooled);
CCONDMOV_OD_Mean_Pooled = [CCONDMOV_O_Mean_Pooled CCONDMOV_D_Mean_Pooled];
CCONDMOV_OD_std_Pooled = [CCONDMOV_O_std_Pooled CCONDMOV_D_std_Pooled];

CCONDFREEZ_O_Mean_Pooled = nanmean(C_CondFreez_O_Pooled);
CCONDFREEZ_O_std_Pooled = nanstd(C_CondFreez_O_Pooled);
CCONDFREEZ_D_Mean_Pooled = nanmean(C_CondFreez_D_Pooled);
CCONDFREEZ_D_std_Pooled = nanstd(C_CondFreez_D_Pooled);
CCONDFREEZ_OD_Mean_Pooled = [CCONDFREEZ_O_Mean_Pooled CCONDFREEZ_D_Mean_Pooled];
CCONDFREEZ_OD_std_Pooled = [CCONDFREEZ_O_std_Pooled CCONDFREEZ_D_std_Pooled];


CPOST_O_Mean_Pooled = nanmean(C_PostSleep_O_Pooled);
CPOST_O_std_Pooled = nanstd(C_PostSleep_O_Pooled);
CPOST_D_Mean_Pooled = nanmean(C_PostSleep_D_Pooled);
CPOST_D_std_Pooled = nanstd(C_PostSleep_D_Pooled);
CPOST_OD_Mean_Pooled = [CPOST_O_Mean_Pooled CPOST_D_Mean_Pooled];
CPOST_OD_std_Pooled = [CPOST_O_std_Pooled CPOST_D_std_Pooled];

CPOSTTEST_O_Mean_Pooled = nanmean(C_PostTest_O_Pooled);
CPOSTTEST_O_std_Pooled = nanstd(C_PostTest_O_Pooled);
CPOSTTEST_D_Mean_Pooled = nanmean(C_PostTest_D_Pooled);
CPOSTTEST_D_std_Pooled = nanstd(C_PostTest_D_Pooled);
CPOSTTEST_OD_Mean_Pooled = [CPOSTTEST_O_Mean_Pooled CPOSTTEST_D_Mean_Pooled];
CPOSTTEST_OD_std_Pooled = [CPOSTTEST_O_std_Pooled CPOSTTEST_D_std_Pooled];

%%%% Stats
[PreStats_Pooled,ppre] = ttest2(C_PreSleep_O_Pooled, C_PreSleep_D_Pooled);
[TaskStats_Pooled,ptask] = ttest2(C_Task_O_Pooled, C_Task_D_Pooled);
[CondMovStats_Pooled,pcondmov] = ttest2(C_CondMov_O_Pooled,C_CondMov_D_Pooled);
[CondFreezStats_Pooled,pcondfreez] = ttest2(C_CondFreez_O_Pooled, C_CondFreez_D_Pooled);
[PostStats_Pooled,ppost] = ttest2(C_PostSleep_O_Pooled,C_PostSleep_D_Pooled);
[PosTestStats_Pooled,pposttest] = ttest2(C_PostTest_O_Pooled, C_PostTest_D_Pooled);


%% Plot averaged over mice

fh = figure('units', 'normalized', 'outerposition', [0 0 0.7 0.7]);
b = barwitherr([CPRE_OD_std_Av;CTASK_OD_std_Av; CPOST_OD_std_Av; CPOSTTEST_OD_std_Av],...
    [CPRE_OD_Mean_Av;CTASK_OD_Mean_Av; CPOST_OD_Mean_Av; CPOSTTEST_OD_Mean_Av]);
set(gca, 'FontSize', 14, 'FontWeight',  'bold');
set(gca, 'LineWidth', 3);
b(1).BarWidth = 0.8;
b(1).FaceColor = 'w';
b(2).FaceColor = 'k';
b(1).LineWidth = 3;
b(2).LineWidth = 3;
x = [b(1).XData + [b(1).XOffset]; b(1).XData - [b(1).XOffset]];
hold on
% set(gca,'Xtick',[1:6],'XtickLabel',{'PreSleep', 'PreExplorations', 'Cond running', 'Cond Freezing', 'PostSleep', 'PostTests'})
set(gca,'Xtick',[1:6],'XtickLabel',{'PreSleep', 'PreExplorations', 'PostSleep', 'PostTests'})
ylabel('Cross-Corellation')
title('Averaged over mice')
hold off
box off
set(gca, 'FontSize', 14, 'FontWeight',  'bold');
set(gca, 'LineWidth', 3);
% title 
lg = legend('Overlapping PCs', 'Non-overlapping PCs');
lg.FontSize = 14;
[p,h5,stats] = signrank([CPRE_O_Mean CPRE_D_Mean]);
if p < 0.05
    sigstar_DB({{1,2}},p,0, 'StarSize',14);
end
[p,h5,stats] = signrank([CTASK_O_Mean CTASK_D_Mean]);
if p < 0.05
    sigstar_DB({{3,4}},p,0, 'StarSize',14);
end

if sav
    saveas(gcf,[dropbox '/MOBS_workingON/Dima/Ongoing results/PlaceField_Final/Pairwise_Small.fig']);
    saveFigure(gcf,'Pairwise_Small',...
        [dropbox '/MOBS_workingON/Dima/Ongoing results/PlaceField_Final/']);
end

%% Plot2 for all cells pooled together

fh = figure('units', 'normalized', 'outerposition', [0 0 0.7 0.7]);
b = barwitherr([CPRE_OD_std_Pooled;CTASK_OD_std_Pooled; CPOST_OD_Mean_Pooled; CPOSTTEST_OD_Mean_Pooled],...
    [CPRE_OD_Mean_Pooled;CTASK_OD_Mean_Pooled; CPOST_OD_Mean_Pooled; CPOSTTEST_OD_Mean_Pooled]);
set(gca, 'FontSize', 14, 'FontWeight',  'bold');
set(gca, 'LineWidth', 3);
b(1).BarWidth = 0.8;
b(1).FaceColor = 'w';
b(2).FaceColor = 'k';
b(1).LineWidth = 3;
b(2).LineWidth = 3;
x = [b(1).XData + [b(1).XOffset]; b(1).XData - [b(1).XOffset]];
hold on
% set(gca,'Xtick',[1:6],'XtickLabel',{'PreSleep', 'PreExplorations', 'Cond running', 'Cond Freezing', 'PostSleep', 'PostTests'})
set(gca,'Xtick',[1:6],'XtickLabel',{'PreSleep', 'PreExplorations', 'PostSleep', 'PostTests'})
ylabel('Cross-Corellation')
title('All cells pooled together')
hold off
box off
set(gca, 'FontSize', 14, 'FontWeight',  'bold');
set(gca, 'LineWidth', 3);
% title 
lg = legend('Overlapping PCs', 'Non-overlapping PCs');
lg.FontSize = 14;
[p,h5,stats] = signrank(CPRE_OD_Mean_Pooled);
if p < 0.05
    sigstar_DB({{1,2}},p,0, 'StarSize',14);
end
[p,h5,stats] = signrank(CTASK_OD_Mean_Pooled);
if p < 0.05
    sigstar_DB({{1,2}},p,0, 'StarSize',14);
end
% sigstar_DB(groups,stats,0,'LineWigth',16,'StarSize',24);


if sav
    saveas(fh,[dropbox '/MOBS_workingON/Dima/Ongoing results/PlaceField_Final/Pairwise_Small_Pooled.fig']);
    saveFigure(fh,'Pairwise_Small_Pooled',...
        [dropbox '/MOBS_workingON/Dima/Ongoing results/PlaceField_Final/']);
end
