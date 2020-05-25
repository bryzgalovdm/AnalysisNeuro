%%%% FRSleepcharac
% 
% A script - go to parameters
% 
% Calculates
% - Firing rate
% - Number of spikes
% - Number of spikes per episode
% - Number of potential stimulation (triggered by a spike with some ISI)
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

% Interstimulus intervals to test for potential stimulations
ISI = [0.1*1e4; 0.5*1e4; 1e4; 2e4; 3e4; 5e4];

% Do you want to save the figures?
savfig = true;

% Paths and names to save
pathfig = '/MOBS_workingON/Dima/Ongoing_results/Spikes/BasicCharacteristics/Number_spikes/'; % Without dropbox

%% Preallocation

% Data
spikes = cell(length(Dir.path), 1);
behav = cell(length(Dir.path), 1);
sleep = cell(length(Dir.path), 1);
% Neurons sets
PCs = cell(length(Dir.path), 1);
nonPCs = cell(length(Dir.path), 1);
% Epochs
NREMPre = cell(length(Dir.path), 1);
REMPre = cell(length(Dir.path), 1);
NREMPost = cell(length(Dir.path), 1);
REMPost = cell(length(Dir.path), 1);
NREMPre_lastH = cell(length(Dir.path), 1);
NREMPost_lastH = cell(length(Dir.path), 1);
% Results - FR
FRNREMPre = cell(length(Dir.path), 1);
FRREMPre = cell(length(Dir.path), 1);
FRNREMPost = cell(length(Dir.path), 1);
FRREMPost = cell(length(Dir.path), 1);
% Results - REM #spikes
numSpREMPre = cell(length(Dir.path), 1);
numSpREMPost = cell(length(Dir.path), 1);
% Results - REM - # spikes per episode
numSpREMPre_ep = cell(length(Dir.path), 1);
numSpREMPost_ep = cell(length(Dir.path), 1);
% Results - NREM #spikes last hour
numSpNREMPre_lastH = cell(length(Dir.path), 1);
numSpNREMPost_lastH = cell(length(Dir.path), 1);
% Results - NREM #FR last hour
FRNREMPre_lastH = cell(length(Dir.path), 1);
FRNREMPost_lastH = cell(length(Dir.path), 1);
% Results - number of stimulations
numStimREM_Pre = cell(length(Dir.path),1);
numStimREM_Post = cell(length(Dir.path),1);
numStimNREM_Pre_lastH = cell(length(Dir.path),1);
numStimNREM_Post_lastH = cell(length(Dir.path),1);

%% Load data - here I handle the exceptions too

for i=1:length(Dir.path)
    spikes{i} = load([Dir.path{i}{1} 'SpikeData.mat'],'S','PlaceCells', 'BasicNeuronInfo');
    behav{i} = load([Dir.path{i}{1} 'behavResources.mat'],'SessionEpoch','CleanVtsd','CleanAlignedXtsd','CleanAlignedYtsd','FreezeAccEpoch');
    if strcmp(Dir.name{i}, 'Mouse906') || strcmp(Dir.name{i}, 'Mouse977') % Mice with bad OB-based sleep scoring
        sleep{i} = load([Dir.path{i}{1} 'SleepScoring_Accelero.mat'],'SWSEpoch','REMEpoch','Sleep'); % Sleep is not used
    else
        sleep{i} = load([Dir.path{i}{1} 'SleepScoring_OBGamma.mat'],'SWSEpoch','REMEpoch','Sleep');  % Sleep is not used
    end
end

%% Split spike arrays into PCs and nonPCs

for i = 1:length(Dir.path)
    if ~isempty(spikes{i}.PlaceCells.idx)
        PCs{i} = spikes{i}.PlaceCells.idx;
        nonPCs{i} = setdiff(find(spikes{i}.BasicNeuronInfo.neuroclass > 0),spikes{i}.BasicNeuronInfo.idx_MUA);
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

%% Calculate firing rate in PCs

for i = 1:length(Dir.path)
    if ~isempty(PCs{i})
        % PreAllocate matrices
        FRNREMPre{i}.PCs = zeros(length(PCs{i}), 1);
        FRREMPre{i}.PCs = zeros(length(PCs{i}), 1);
        FRNREMPost{i}.PCs = zeros(length(PCs{i}), 1);
        FRREMPost{i}.PCs = zeros(length(PCs{i}), 1);
        
        FRNREMPre{i}.nonPCs = zeros(length(nonPCs{i}), 1);
        FRREMPre{i}.nonPCs = zeros(length(nonPCs{i}), 1);
        FRNREMPost{i}.nonPCs = zeros(length(nonPCs{i}), 1);
        FRREMPost{i}.nonPCs = zeros(length(nonPCs{i}), 1);
        
        FRNREMPre_lastH{i}.PCs = nan(length(PCs{i}), 1);
        FRNREMPre_lastH{i}.nonPCs = nan(length(nonPCs{i}), 1);
        FRNREMPost_lastH{i}.PCs = nan(length(PCs{i}), 1);
        FRNREMPost_lastH{i}.nonPCs = nan(length(nonPCs{i}), 1);
        
        % PlaceCells
        for j = 1:length(PCs{i})
            FRNREMPre{i}.PCs(j) = length(Restrict(spikes{i}.S{PCs{i}(j)}, NREMPre{i}))/sum(End(NREMPre{i}, 's') - Start(NREMPre{i}, 's'));
            FRREMPre{i}.PCs(j) = length(Restrict(spikes{i}.S{PCs{i}(j)}, REMPre{i}))/sum(End(REMPre{i}, 's') - Start(REMPre{i}, 's'));
            FRNREMPost{i}.PCs(j) = length(Restrict(spikes{i}.S{PCs{i}(j)}, NREMPost{i}))/sum(End(NREMPost{i}, 's') - Start(NREMPost{i}, 's'));
            FRREMPost{i}.PCs(j) = length(Restrict(spikes{i}.S{PCs{i}(j)}, REMPost{i}))/sum(End(REMPost{i}, 's') - Start(REMPost{i}, 's'));
            if ~isempty(NREMPre_lastH{i})
                FRNREMPre_lastH{i}.PCs(j) = length(Restrict(spikes{i}.S{PCs{i}(j)}, NREMPre_lastH{i}))/sum(End(NREMPre_lastH{i}, 's') - Start(NREMPre_lastH{i}, 's'));
            end
            if ~isempty(NREMPost_lastH{i})
                FRNREMPost_lastH{i}.PCs(j) = length(Restrict(spikes{i}.S{PCs{i}(j)}, NREMPost_lastH{i}))/sum(End(NREMPost_lastH{i}, 's') - Start(NREMPost_lastH{i}, 's'));
            end
        end
        
        % non-PlaceCells
        for j = 1:length(nonPCs{i})
            FRNREMPre{i}.nonPCs(j) = length(Restrict(spikes{i}.S{nonPCs{i}(j)}, NREMPre{i}))/sum(End(NREMPre{i}, 's') - Start(NREMPre{i}, 's'));
            FRREMPre{i}.nonPCs(j) = length(Restrict(spikes{i}.S{nonPCs{i}(j)}, REMPre{i}))/sum(End(REMPre{i}, 's') - Start(REMPre{i}, 's'));
            FRNREMPost{i}.nonPCs(j) = length(Restrict(spikes{i}.S{nonPCs{i}(j)}, NREMPost{i}))/sum(End(NREMPost{i}, 's') - Start(NREMPost{i}, 's'));
            FRREMPost{i}.nonPCs(j) = length(Restrict(spikes{i}.S{nonPCs{i}(j)}, REMPost{i}))/sum(End(REMPost{i}, 's') - Start(REMPost{i}, 's'));
            if ~isempty(NREMPre_lastH{i})
                FRNREMPre_lastH{i}.nonPCs(j) = length(Restrict(spikes{i}.S{nonPCs{i}(j)}, NREMPre_lastH{i}))/sum(End(NREMPre_lastH{i}, 's') - Start(NREMPre_lastH{i}, 's'));
            end
            if ~isempty(NREMPost_lastH{i})
                FRNREMPost_lastH{i}.nonPCs(j) = length(Restrict(spikes{i}.S{nonPCs{i}(j)}, NREMPost_lastH{i}))/sum(End(NREMPost_lastH{i}, 's') - Start(NREMPost_lastH{i}, 's'));
            end
        end
    end
end

%% Calculate number of spikes in REM - first two hours of session

for i = 1:length(Dir.path)
    if ~isempty(PCs{i})
        % PreAllocate matrices
        numSpREMPre{i}.PCs = zeros(length(PCs{i}), 1);
        numSpREMPost{i}.PCs = zeros(length(PCs{i}), 1);
        numSpREMPre_ep{i}.PCs = zeros(length(PCs{i}), 1);
        numSpREMPost_ep{i}.PCs = zeros(length(PCs{i}), 1);
        
        numSpREMPre{i}.nonPCs = zeros(length(nonPCs{i}), 1);
        numSpREMPost{i}.nonPCs = zeros(length(nonPCs{i}), 1);
        numSpREMPre_ep{i}.nonPCs = zeros(length(nonPCs{i}), 1);
        numSpREMPost_ep{i}.nonPCs = zeros(length(nonPCs{i}), 1);
        
        % PlaceCells
        for j = 1:length(PCs{i})
            numSpREMPre{i}.PCs(j) = length(Restrict(spikes{i}.S{PCs{i}(j)}, REMPre{i}));
            numSpREMPost{i}.PCs(j) = length(Restrict(spikes{i}.S{PCs{i}(j)}, REMPost{i}));
            
            numSpREMPre_ep{i}.PCs(j) = length(Restrict(spikes{i}.S{PCs{i}(j)}, REMPre{i}))/length(Start(REMPre{i}));
            numSpREMPost_ep{i}.PCs(j) = length(Restrict(spikes{i}.S{PCs{i}(j)}, REMPost{i}))/length(Start(REMPost{i}));
        end
        
        % non-PlaceCells
        for j = 1:length(nonPCs{i})
            numSpREMPre{i}.nonPCs(j) = length(Restrict(spikes{i}.S{nonPCs{i}(j)}, REMPre{i}));
            numSpREMPost{i}.nonPCs(j) = length(Restrict(spikes{i}.S{nonPCs{i}(j)}, REMPost{i}));
            
            numSpREMPre_ep{i}.nonPCs(j) = length(Restrict(spikes{i}.S{nonPCs{i}(j)}, REMPre{i}))/length(Start(REMPre{i}));
            numSpREMPost_ep{i}.nonPCs(j) = length(Restrict(spikes{i}.S{nonPCs{i}(j)}, REMPost{i}))/length(Start(REMPost{i}));
        end
    end
end

%% Calculate number of spikes in the last hour of NREM

for i = 1:length(Dir.path)
    if ~isempty(PCs{i})
        % PreAllocate matrices
        numSpNREMPre_lastH{i}.PCs = zeros(length(PCs{i}), 1);
        numSpNREMPost_lastH{i}.PCs = zeros(length(PCs{i}), 1);
        
        numSpNREMPre_lastH{i}.nonPCs = zeros(length(nonPCs{i}), 1);
        numSpNREMPost_lastH{i}.nonPCs = zeros(length(nonPCs{i}), 1);
        
        % PlaceCells
        for j = 1:length(PCs{i})
            numSpNREMPre_lastH{i}.PCs(j) = length(Restrict(spikes{i}.S{PCs{i}(j)}, NREMPre_lastH{i}));
            numSpNREMPost_lastH{i}.PCs(j) = length(Restrict(spikes{i}.S{PCs{i}(j)}, NREMPost_lastH{i}));
        end
        
        % non-PlaceCells
        for j=1:length(nonPCs{i})
            numSpNREMPre_lastH{i}.nonPCs(j) = length(Restrict(spikes{i}.S{nonPCs{i}(j)}, NREMPre_lastH{i}));
            numSpNREMPost_lastH{i}.nonPCs(j) = length(Restrict(spikes{i}.S{nonPCs{i}(j)}, NREMPost_lastH{i}));
        end
    end
end

%% Calculate number of stimulations given the activity of each cell
% REM
for i=1:length(Dir.path)
    % PreAllocate
    numStimREM_Pre{i}.PCs = zeros(length(PCs{i}), length(ISI));
    numStimREM_Pre{i}.nonPCs = zeros(length(nonPCs{i}), length(ISI));
    numStimREM_Post{i}.PCs = zeros(length(PCs{i}), length(ISI));
    numStimREM_Post{i}.nonPCs = zeros(length(nonPCs{i}), length(ISI));
    
    numStimNREM_Pre_lastH{i}.PCs = nan(length(PCs{i}), length(ISI));
    numStimNREM_Pre_lastH{i}.nonPCs = nan(length(nonPCs{i}), length(ISI));
    numStimNREM_Post_lastH{i}.PCs = nan(length(PCs{i}), length(ISI));
    numStimNREM_Post_lastH{i}.nonPCs = nan(length(nonPCs{i}), length(ISI));
    
    
    % PCs
    for j = 1:length(PCs{i})
        spiketrain_preREM = Restrict(spikes{i}.S{PCs{i}(j)}, REMPre{i});
        spiketrain_postREM = Restrict(spikes{i}.S{PCs{i}(j)}, REMPost{i});
        spiketrain_preNREM_lastH = Restrict(spikes{i}.S{PCs{i}(j)}, NREMPre_lastH{i});
        spiketrain_postNREM_lastH = Restrict(spikes{i}.S{PCs{i}(j)}, NREMPost_lastH{i});
        for k = 1:length(ISI)
            % PreREM
            tempdiff = diff(Data(spiketrain_preREM));
            numStimREM_Pre{i}.PCs(j,k) = 1; % Stimulation to the first stimulus
            count = 0;
            for l = 1:length(tempdiff)
                count = count + tempdiff(l);
                if count > ISI(k)
                    count = 0;
                    numStimREM_Pre{i}.PCs(j,k) = numStimREM_Pre{i}.PCs(j,k) + 1;
                end
            end
            clear tempdiff count
            % PostREM
            tempdiff = diff(Data(spiketrain_postREM));
            numStimREM_Post{i}.PCs(j,k) = 1;
            count = 0;
            for l = 1:length(tempdiff)
                count = count + tempdiff(l);
                if count > ISI(k)
                    count = 0;
                    numStimREM_Post{i}.PCs(j,k) = numStimREM_Post{i}.PCs(j,k) + 1;
                end
            end
            clear tempdiff count
            
            % PreNREM
            if ~isempty(NREMPre_lastH{i})
                tempdiff = diff(Data(spiketrain_preNREM_lastH));
                numStimNREM_Pre_lastH{i}.PCs(j,k) = 1;
                count = 0;
                for l = 1:length(tempdiff)
                    count = count + tempdiff(l);
                    if count > ISI(k)
                        count = 0;
                        numStimNREM_Pre_lastH{i}.PCs(j,k) = numStimNREM_Pre_lastH{i}.PCs(j,k) + 1;
                    end
                end
                clear tempdiff count
            end
            % PostNREM
            if ~isempty(NREMPost_lastH{i})
                tempdiff = diff(Data(spiketrain_postNREM_lastH));
                numStimNREM_Post_lastH{i}.PCs(j,k) = 1;
                count = 0;
                for l = 1:length(tempdiff)
                    count = count + tempdiff(l);
                    if count > ISI(k)
                        count = 0;
                        numStimNREM_Post_lastH{i}.PCs(j,k) = numStimNREM_Post_lastH{i}.PCs(j,k) + 1;
                    end
                end
                clear tempdiff count
            end
        end
        clear spiketrain_preREM spiketrain_postREM spiketrain_preNREM_lastH spiketrain_postNREM_lastH
    end
    
    % non-PCs
    for j = 1:length(nonPCs{i})
        spiketrain_preREM = Restrict(spikes{i}.S{nonPCs{i}(j)}, REMPre{i});
        spiketrain_postREM = Restrict(spikes{i}.S{nonPCs{i}(j)}, REMPost{i});
        spiketrain_preNREM_lastH = Restrict(spikes{i}.S{nonPCs{i}(j)}, NREMPre_lastH{i});
        spiketrain_postNREM_lastH = Restrict(spikes{i}.S{nonPCs{i}(j)}, NREMPost_lastH{i});
        for k = 1:length(ISI)
            % PreREM
            tempdiff = diff(Data(spiketrain_preREM));
            numStimREM_Pre{i}.nonPCs(j,k) = 1; % Stimulation to the first stimulus
            count = 0;
            for l = 1:length(tempdiff)
                count = count + tempdiff(l);
                if count > ISI(k)
                    count = 0;
                    numStimREM_Pre{i}.nonPCs(j,k) = numStimREM_Pre{i}.nonPCs(j,k) + 1;
                end
            end
            clear tempdiff count
            % PostREM
            tempdiff = diff(Data(spiketrain_postREM));
            numStimREM_Post{i}.nonPCs(j,k) = 1;
            count = 0;
            for l = 1:length(tempdiff)
                count = count + tempdiff(l);
                if count > ISI(k)
                    count = 0;
                    numStimREM_Post{i}.nonPCs(j,k) = numStimREM_Post{i}.nonPCs(j,k) + 1;
                end
            end
            clear tempdiff count
            
            % PreNREM
            if ~isempty(NREMPre_lastH{i})
                tempdiff = diff(Data(spiketrain_preNREM_lastH));
                numStimNREM_Pre_lastH{i}.nonPCs(j,k) = 1;
                count = 0;
                for l = 1:length(tempdiff)
                    count = count + tempdiff(l);
                    if count > ISI(k)
                        count = 0;
                        numStimNREM_Pre_lastH{i}.nonPCs(j,k) = numStimNREM_Pre_lastH{i}.nonPCs(j,k) + 1;
                    end
                end
                clear tempdiff count
            end
            % PostNREM
            if ~isempty(NREMPost_lastH{i})
                tempdiff = diff(Data(spiketrain_postNREM_lastH));
                numStimNREM_Post_lastH{i}.nonPCs(j,k) = 1;
                count = 0;
                for l = 1:length(tempdiff)
                    count = count + tempdiff(l);
                    if count > ISI(k)
                        count = 0;
                        numStimNREM_Post_lastH{i}.nonPCs(j,k) = numStimNREM_Post_lastH{i}.nonPCs(j,k) + 1;
                    end
                end
                clear tempdiff count
            end
        end
        clear spiketrain_preREM spiketrain_postREM spiketrain_preNREM_lastH spiketrain_postNREM_lastH
    end
    
end

%% Pool neurons

% FR of neurons
FR.Pre.NREM.PCs = FRNREMPre{1}.PCs;
FR.Pre.NREM.nonPCs = FRNREMPre{1}.nonPCs;
FR.Pre.REM.PCs = FRREMPre{1}.PCs;
FR.Pre.REM.nonPCs = FRREMPre{1}.nonPCs;

FR.Post.NREM.PCs = FRNREMPost{1}.PCs;
FR.Post.NREM.nonPCs = FRNREMPost{1}.nonPCs;
FR.Post.REM.PCs = FRREMPost{1}.PCs;
FR.Post.REM.nonPCs = FRREMPost{1}.nonPCs;

% # spikes in REM
numSp.Pre.REM.PCs = numSpREMPre{1}.PCs;
numSp.Pre.REM.nonPCs = numSpREMPre{1}.nonPCs;
numSp_ep.Pre.REM.PCs = numSpREMPre_ep{1}.PCs;
numSp_ep.Pre.REM.nonPCs = numSpREMPre_ep{1}.nonPCs;

numSp.Post.REM.PCs = numSpREMPost{1}.PCs;
numSp.Post.REM.nonPCs = numSpREMPost{1}.nonPCs;
numSp_ep.Post.REM.PCs = numSpREMPost_ep{1}.PCs;
numSp_ep.Post.REM.nonPCs = numSpREMPost_ep{1}.nonPCs;

% Spikes in NREM
numSp.Pre.NREM.PCs = numSpNREMPre_lastH{1}.PCs;
numSp.Pre.NREM.nonPCs = numSpNREMPre_lastH{1}.nonPCs;
numSp.Post.NREM.PCs = numSpNREMPost_lastH{1}.PCs;
numSp.Post.NREM.nonPCs = numSpNREMPost_lastH{1}.nonPCs;

% Spikes in NREM (last hour)
FR.Pre.NREM_lastH.PCs = FRNREMPre_lastH{1}.PCs;
FR.Pre.NREM_lastH.nonPCs = FRNREMPre_lastH{1}.nonPCs;
FR.Post.NREM_lastH.PCs = FRNREMPost_lastH{1}.PCs;
FR.Post.NREM_lastH.nonPCs = FRNREMPost_lastH{1}.nonPCs;

% Stimulation numbers
numStim.Pre.REM.PCs = numStimREM_Pre{1}.PCs;
numStim.Post.REM.PCs = numStimREM_Post{1}.PCs;
numStim.Pre.REM.nonPCs = numStimREM_Pre{1}.nonPCs;
numStim.Post.REM.nonPCs = numStimREM_Post{1}.nonPCs;

numStim.Pre.NREM_last.PCs = numStimNREM_Pre_lastH{1}.PCs;
numStim.Post.NREM_last.PCs = numStimNREM_Post_lastH{1}.PCs;
numStim.Pre.NREM_last.nonPCs = numStimNREM_Pre_lastH{1}.nonPCs;
numStim.Post.NREM_last.nonPCs = numStimNREM_Post_lastH{1}.nonPCs;

if length(Dir.path) > 1
    for i = 2:length(Dir.path)
        if ~isempty(PCs{i})
            FR.Pre.NREM.PCs = [FR.Pre.NREM.PCs; FRNREMPre{i}.PCs];
            FR.Pre.NREM.nonPCs = [FR.Pre.NREM.nonPCs; FRNREMPre{i}.nonPCs];
            FR.Pre.REM.PCs = [FR.Pre.REM.PCs; FRREMPre{i}.PCs];
            FR.Pre.REM.nonPCs = [FR.Pre.REM.nonPCs; FRREMPre{i}.nonPCs];
            
            FR.Post.NREM.PCs = [FR.Post.NREM.PCs; FRNREMPost{i}.PCs];
            FR.Post.NREM.nonPCs = [FR.Post.NREM.nonPCs; FRNREMPost{i}.nonPCs];
            FR.Post.REM.PCs = [FR.Post.REM.PCs; FRREMPost{i}.PCs];
            FR.Post.REM.nonPCs = [FR.Post.REM.nonPCs; FRREMPost{i}.nonPCs];
            
            numSp.Pre.REM.PCs = [numSp.Pre.REM.PCs; numSpREMPre{i}.PCs];
            numSp.Pre.REM.nonPCs = [numSp.Pre.REM.nonPCs; numSpREMPre{i}.nonPCs];
            numSp_ep.Pre.REM.PCs = [numSp_ep.Pre.REM.PCs; numSpREMPre_ep{i}.PCs];
            numSp_ep.Pre.REM.nonPCs = [numSp_ep.Pre.REM.nonPCs; numSpREMPre_ep{i}.nonPCs];
            
            numSp.Post.REM.PCs = [numSp.Post.REM.PCs; numSpREMPost{i}.PCs];
            numSp.Post.REM.nonPCs = [numSp.Post.REM.nonPCs; numSpREMPost{i}.nonPCs];
            numSp_ep.Post.REM.PCs = [numSp_ep.Post.REM.PCs; numSpREMPost_ep{i}.PCs];
            numSp_ep.Post.REM.nonPCs = [numSp_ep.Post.REM.nonPCs; numSpREMPost_ep{i}.nonPCs];
            
            numSp.Pre.NREM.PCs = [numSp.Pre.NREM.PCs; numSpNREMPre_lastH{i}.PCs];
            numSp.Pre.NREM.nonPCs = [numSp.Pre.NREM.nonPCs; numSpNREMPre_lastH{i}.nonPCs];
            numSp.Post.NREM.PCs = [numSp.Post.NREM.PCs; numSpNREMPost_lastH{i}.PCs];
            numSp.Post.NREM.nonPCs = [numSp.Post.NREM.nonPCs; numSpNREMPost_lastH{i}.nonPCs];
            
            FR.Pre.NREM_lastH.PCs = [FR.Pre.NREM_lastH.PCs; FRNREMPre_lastH{i}.PCs];
            FR.Pre.NREM_lastH.nonPCs = [FR.Pre.NREM_lastH.nonPCs; FRNREMPre_lastH{i}.nonPCs];
            FR.Post.NREM_lastH.PCs = [FR.Post.NREM_lastH.PCs; FRNREMPost_lastH{i}.PCs];
            FR.Post.NREM_lastH.nonPCs = [FR.Post.NREM_lastH.PCs; FRNREMPost_lastH{i}.nonPCs];
            
            numStim.Pre.REM.PCs = [numStim.Pre.REM.PCs; numStimREM_Pre{i}.PCs];
            numStim.Post.REM.PCs = [numStim.Post.REM.PCs; numStimREM_Post{i}.PCs];
            numStim.Pre.REM.nonPCs = [numStim.Pre.REM.nonPCs; numStimREM_Pre{i}.nonPCs];
            numStim.Post.REM.nonPCs = [numStim.Post.REM.nonPCs; numStimREM_Post{i}.nonPCs];
            
            numStim.Pre.NREM_last.PCs = [numStim.Pre.NREM_last.PCs; numStimNREM_Pre_lastH{i}.PCs];
            numStim.Post.NREM_last.PCs = [numStim.Post.NREM_last.PCs; numStimNREM_Post_lastH{i}.PCs];
            numStim.Pre.NREM_last.nonPCs = [numStim.Pre.NREM_last.nonPCs; numStimNREM_Pre_lastH{i}.nonPCs];
            numStim.Post.NREM_last.nonPCs = [numStim.Post.NREM_last.nonPCs; numStimNREM_Post_lastH{i}.nonPCs];
        end
        
    end
end

%% Get shock zone PCs idxs for plotting
[numPCs,numPCs_M] = CountPlaceCells(Dir, 'Verbose', false);
idx_SZ = (arrayfun(@(x)find(spikes{1}.PlaceCells.idx==x,1),spikes{1}.PlaceCells.SZ))';
for i=2:length(Dir.path)
    idx_SZ = [idx_SZ; (arrayfun(@(x)find(spikes{i}.PlaceCells.idx==x,1),spikes{i}.PlaceCells.SZ)+sum(numPCs_M(1:i-1)))'];
end

%% Illustrate the data

%%%%%%%%%%% REM SLEEP %%%%%%%%%%%%%%%%%%%%%%%%%
% FR during REM
Pl = {FR.Pre.REM.PCs; FR.Post.REM.PCs; []; FR.Pre.REM.nonPCs; FR.Post.REM.nonPCs};
Cols = {[0.7 0.7 0.9], [0.2 0.2 0.9],[], [0.45 0.6 0.3], [0.3 0.4 0.2]};
ffr = figure('units', 'normalized', 'outerposition', [0 0 0.5 0.95]);
MakeViolinPlot_DB(Pl,Cols,1:5,[],0);
hold on
for i=1:2
    if ~isempty(Pl{i})
        handlesplot=plotSpread(Pl{i}(idx_SZ),'distributionColors',[1 0 0],'xValues',i,'spreadWidth',0.8);
        set(handlesplot{1},'MarkerSize',15);
    end
end
hold off
set(gca,'XTick', [1:5], 'XTickLabel', {'PreREM','PostREM', '', 'PreREM', 'PostREM'}, 'FontSize', 19, 'FontWeight', 'bold');
yl = ylim;
%//
% Add groups - code stolen from <stackoverflow.com/questions/33165830/double-ticklabel-in-matlab>
groupX = [1.5 4.5]; %// central value of each group
groupY = yl(1) - yl(2)*0.08; %// vertical position of texts. Adjust as needed
deltaY = .03; %// controls vertical compression of axis. Adjust as needed
groupNames = {'PlaceCells', 'non Place Cells'};
for g = 1:numel(groupX)
    h = text(groupX(g), groupY, groupNames{g}, 'Fontsize', 23, 'Fontweight', 'bold');
    %// create text for group with appropriate font size and weight
    pos = get(h, 'Position');
    ext = get(h, 'Extent');
    pos(1) = pos(1) - ext(3)/2; %// horizontally correct position to make it centered
    set(h, 'Position', pos); %// set corrected position for text
end
pos = get(gca, 'position');
pos(2) = pos(2) + deltaY; %// vertically compress axis to make room for texts
set(gca, 'Position', pos); %/ set corrected position for axis
%//
ylabel('Firing rate (Hz)')
title(['Firing rate during REM sleep, N(PCs)=' num2str(numPCs) ', N(nonPCs)=' num2str(length(FR.Pre.REM.nonPCs))],...
    'FontSize',16);
if savfig
    saveas(ffr,[dropbox pathfig 'FR_REM.fig']);
    saveFigure(ffr,'FR_REM',[dropbox pathfig]);
end

% Number of spikes during REM
Pl = {numSp.Pre.REM.PCs; numSp.Post.REM.PCs; []; numSp.Pre.REM.nonPCs; numSp.Post.REM.nonPCs};
Cols = {[0.7 0.7 0.9], [0.2 0.2 0.9],[], [0.45 0.6 0.3], [0.3 0.4 0.2]};
fh = figure('units', 'normalized', 'outerposition', [0 0 0.5 0.95]);
MakeViolinPlot_DB(Pl,Cols,1:5,[],0);
hold on
for i=1:2
    if ~isempty(Pl{i})
        handlesplot=plotSpread(Pl{i}(idx_SZ),'distributionColors',[1 0 0],'xValues',i,'spreadWidth',0.8);
        set(handlesplot{1},'MarkerSize',15);
    end
end
hold off
set(gca,'XTick', [1:5], 'XTickLabel', {'PreREM','PostREM', '', 'PreREM', 'PostREM'}, 'FontSize', 19, 'FontWeight', 'bold');
yl = ylim;
%//
% Add groups - code stolen from <stackoverflow.com/questions/33165830/double-ticklabel-in-matlab>
groupX = [1.5 4.5]; %// central value of each group
groupY = yl(1) - yl(2)*0.08; %// vertical position of texts. Adjust as needed
deltaY = .03; %// controls vertical compression of axis. Adjust as needed
groupNames = {'PlaceCells', 'non Place Cells'};
for g = 1:numel(groupX)
    h = text(groupX(g), groupY, groupNames{g}, 'Fontsize', 23, 'Fontweight', 'bold');
    %// create text for group with appropriate font size and weight
    pos = get(h, 'Position');
    ext = get(h, 'Extent');
    pos(1) = pos(1) - ext(3)/2; %// horizontally correct position to make it centered
    set(h, 'Position', pos); %// set corrected position for text
end
pos = get(gca, 'position');
pos(2) = pos(2) + deltaY; %// vertically compress axis to make room for texts
set(gca, 'Position', pos); %/ set corrected position for axis
%//
ylabel('Number of spikes')
title(['Number of spikes during REM sleep, N(PCs)=' num2str(numPCs) ', N(nonPCs)=' num2str(length(FR.Pre.REM.nonPCs))],...
    'FontSize',16);
if savfig
    saveas(fh,[dropbox pathfig 'numSp_REM.fig']);
    saveFigure(fh,'numSp_REM',[dropbox pathfig]);
end

% Number of spikes per episode
Pl = {numSp_ep.Pre.REM.PCs; numSp_ep.Post.REM.PCs; []; numSp_ep.Pre.REM.nonPCs; numSp_ep.Post.REM.nonPCs};
Cols = {[0.7 0.7 0.9], [0.2 0.2 0.9],[], [0.45 0.6 0.3], [0.3 0.4 0.2]};
fh_ep = figure('units', 'normalized', 'outerposition', [0 0 0.5 0.95]);
MakeViolinPlot_DB(Pl,Cols,1:5,[],0);
hold on
for i=1:2
    if ~isempty(Pl{i})
        handlesplot=plotSpread(Pl{i}(idx_SZ),'distributionColors',[1 0 0],'xValues',i,'spreadWidth',0.8);
        set(handlesplot{1},'MarkerSize',15);
    end
end
hold off
set(gca,'XTick', [1:5], 'XTickLabel', {'PreREM','PostREM', '', 'PreREM', 'PostREM'}, 'FontSize', 19, 'FontWeight', 'bold');
yl = ylim;
%//
% Add groups - code stolen from <stackoverflow.com/questions/33165830/double-ticklabel-in-matlab>
groupX = [1.5 4.5]; %// central value of each group
groupY = yl(1) - yl(2)*0.08; %// vertical position of texts. Adjust as needed
deltaY = .03; %// controls vertical compression of axis. Adjust as needed
groupNames = {'PlaceCells', 'non Place Cells'};
for g = 1:numel(groupX)
    h = text(groupX(g), groupY, groupNames{g}, 'Fontsize', 23, 'Fontweight', 'bold');
    %// create text for group with appropriate font size and weight
    pos = get(h, 'Position');
    ext = get(h, 'Extent');
    pos(1) = pos(1) - ext(3)/2; %// horizontally correct position to make it centered
    set(h, 'Position', pos); %// set corrected position for text
end
pos = get(gca, 'position');
pos(2) = pos(2) + deltaY; %// vertically compress axis to make room for texts
set(gca, 'Position', pos); %/ set corrected position for axis
%//
ylabel('Number of spikes')
title(['Number of spikes per REM episode, N(PCs)=' num2str(numPCs) ', N(nonPCs)=' num2str(length(FR.Pre.REM.nonPCs))],...
    'FontSize',16);
if savfig
    saveas(fh_ep,[dropbox pathfig 'numSpperEp_REM.fig']);
    saveFigure(fh_ep,'numSpperEp_REM',[dropbox pathfig]);
end


% Number of stimulations in REM
Pl = {numStim.Pre.REM.PCs, numStim.Pre.REM.PCs};
Cols = {[0.7 0.7 0.9], [0.2 0.2 0.9]};
tits = {['PreREM: all PCs (N=' num2str(numPCs) ')'],...
    ['PostREM: all PCs (N=' num2str(numPCs) ')']};
f1 = figure('units', 'normalized', 'outerposition', [0 0 0.25 0.75]);
ax = arrayfun(@(i) subplot(2,1,i, 'NextPlot', 'add', 'Box', 'off'), [1:2]);
for i = 1:length(ax)
    axes(ax(i));
    shadedErrorBar(ISI,nanmean(Pl{i},1),nanstd(Pl{i},1),...
        {'-ok','markerfacecolor',[0.7,0.4,0.3], 'Color', Cols{i}});
    xlim([ISI(1) ISI(end)]);
    set(gca, 'FontSize', 16, 'FontWeight', 'bold','XTick', ISI, 'XTickLabel',...
        {'100 ms', '500 ms', '1 s', '2s', '3s', '5s'});
    set(gca,'xscale','log');
    xlabel('log ISI');
    ylabel('# Stim');
    title(tits{i}, 'FontSize', 14);
end
yl = ylim;
if savfig
    saveas(f1,[dropbox pathfig 'numStim_REM.fig']);
    saveFigure(f1,'numStim_REM',[dropbox pathfig]);
end
% Number of stimulations in REM - separate figure for only SZ cells
Pl = {numStim.Pre.REM.PCs(idx_SZ,:), numStim.Pre.REM.PCs(idx_SZ,:)};
Cols = {[0.9 0.7 0.7], [0.9 0.2 0.2]};
tits = {['PreREM: shock zone PCs (N=' num2str(length(idx_SZ)) ')'],...
    ['PostREM: shock zone PCs (N=' num2str(length(idx_SZ)) ')']};
f2 = figure('units', 'normalized', 'outerposition', [0 0 0.25 0.75]);
ax = arrayfun(@(i) subplot(2,1,i, 'NextPlot', 'add', 'Box', 'off'), [1:2]);
for i = 1:length(ax)
    axes(ax(i));
    shadedErrorBar(ISI,nanmean(Pl{i},1),nanstd(Pl{i},1),...
        {'-ok','markerfacecolor',[0.7,0.4,0.3], 'Color', Cols{i}});
    xlim([ISI(1) ISI(end)]);
    set(gca, 'FontSize', 16, 'FontWeight', 'bold','XTick', ISI, 'XTickLabel',...
        {'100 ms', '500 ms', '1 s', '2s', '3s', '5s'});
    set(gca,'xscale','log');
    xlabel('log ISI');
    ylabel('# Stim');
    title(tits{i}, 'FontSize', 14);
    ylim(yl);
end
if savfig
    saveas(f2,[dropbox pathfig 'numStim_REM_SZ.fig']);
    saveFigure(f2,'numStim_REM_SZ',[dropbox pathfig]);
end

%%%%%%%%%%% REM SLEEP %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% NREM SLEEP %%%%%%%%%%%%%%%%%%%%%%%%%
% FR during NREM
Pl = {FR.Pre.NREM.PCs; FR.Post.NREM.PCs; []; FR.Pre.NREM.nonPCs; FR.Post.NREM.nonPCs};
Cols = {[0.7 0.7 0.9], [0.2 0.2 0.9],[], [0.45 0.6 0.3], [0.3 0.4 0.2]};
ffr_nrem = figure('units', 'normalized', 'outerposition', [0 0 0.5 0.95]);
MakeViolinPlot_DB(Pl,Cols,1:5,[],0);
hold on
for i=1:2
    if ~isempty(Pl{i})
        handlesplot=plotSpread(Pl{i}(idx_SZ),'distributionColors',[1 0 0],'xValues',i,'spreadWidth',0.8);
        set(handlesplot{1},'MarkerSize',15);
    end
end
hold off
set(gca,'XTick', [1:5], 'XTickLabel', {'PreREM','PostREM', '', 'PreREM', 'PostREM'}, 'FontSize', 19, 'FontWeight', 'bold');
yl = ylim;
%//
% Add groups - code stolen from <stackoverflow.com/questions/33165830/double-ticklabel-in-matlab>
groupX = [1.5 4.5]; %// central value of each group
groupY = yl(1) - yl(2)*0.08; %// vertical position of texts. Adjust as needed
deltaY = .03; %// controls vertical compression of axis. Adjust as needed
groupNames = {'PlaceCells', 'non Place Cells'};
for g = 1:numel(groupX)
    h = text(groupX(g), groupY, groupNames{g}, 'Fontsize', 23, 'Fontweight', 'bold');
    %// create text for group with appropriate font size and weight
    pos = get(h, 'Position');
    ext = get(h, 'Extent');
    pos(1) = pos(1) - ext(3)/2; %// horizontally correct position to make it centered
    set(h, 'Position', pos); %// set corrected position for text
end
pos = get(gca, 'position');
pos(2) = pos(2) + deltaY; %// vertically compress axis to make room for texts
set(gca, 'Position', pos); %/ set corrected position for axis
%//
ylabel('Firing rate (Hz)')
title(['Firing rate during whole NREM sleep, N(PCs)=' num2str(numPCs) ', N(nonPCs)=' num2str(length(FR.Pre.REM.nonPCs))],...
    'FontSize',16);
if savfig
    saveas(ffr_nrem,[dropbox pathfig 'FR_NREM.fig']);
    saveFigure(ffr_nrem,'FR_NREM',[dropbox pathfig]);
end

% FR during last hour of NREM
Pl = {FR.Pre.NREM_lastH.PCs; FR.Post.NREM_lastH.PCs; []; FR.Pre.NREM_lastH.nonPCs; FR.Post.NREM_lastH.nonPCs};
Cols = {[0.7 0.7 0.9], [0.2 0.2 0.9],[], [0.45 0.6 0.3], [0.3 0.4 0.2]};
ffr_lh = figure('units', 'normalized', 'outerposition', [0 0 0.5 0.95]);
MakeViolinPlot_DB(Pl,Cols,1:5,[],0);
hold on
for i=1:2
    if ~isempty(Pl{i})
        handlesplot=plotSpread(Pl{i}(idx_SZ),'distributionColors',[1 0 0],'xValues',i,'spreadWidth',0.8);
        set(handlesplot{1},'MarkerSize',15);
    end
end
hold off
set(gca,'XTick', [1:5], 'XTickLabel', {'PreREM','PostREM', '', 'PreREM', 'PostREM'}, 'FontSize', 19, 'FontWeight', 'bold');
yl = ylim;
%//
% Add groups - code stolen from <stackoverflow.com/questions/33165830/double-ticklabel-in-matlab>
groupX = [1.5 4.5]; %// central value of each group
groupY = yl(1) - yl(2)*0.08; %// vertical position of texts. Adjust as needed
deltaY = .03; %// controls vertical compression of axis. Adjust as needed
groupNames = {'PlaceCells', 'non Place Cells'};
for g = 1:numel(groupX)
    h = text(groupX(g), groupY, groupNames{g}, 'Fontsize', 23, 'Fontweight', 'bold');
    %// create text for group with appropriate font size and weight
    pos = get(h, 'Position');
    ext = get(h, 'Extent');
    pos(1) = pos(1) - ext(3)/2; %// horizontally correct position to make it centered
    set(h, 'Position', pos); %// set corrected position for text
end
pos = get(gca, 'position');
pos(2) = pos(2) + deltaY; %// vertically compress axis to make room for texts
set(gca, 'Position', pos); %/ set corrected position for axis
%//
ylabel('Firing rate (Hz)')
title(['Firing rate during second hour of NREM sleep, N(PCs)=' num2str(numPCs) ', N(nonPCs)=' num2str(length(FR.Pre.REM.nonPCs))],...
    'FontSize',16);
if savfig
    saveas(ffr_lh,[dropbox pathfig 'FR_NREM_lastH.fig']);
    saveFigure(ffr_lh,'FR_NREM_lastH',[dropbox pathfig]);
end

% Number of spikes during NREM
Pl = {numSp.Pre.NREM.PCs; numSp.Post.NREM.PCs; []; numSp.Pre.NREM.nonPCs; numSp.Post.NREM.nonPCs};
Cols = {[0.7 0.7 0.9], [0.2 0.2 0.9],[], [0.45 0.6 0.3], [0.3 0.4 0.2]};
f_nrem = figure('units', 'normalized', 'outerposition', [0 0 0.5 0.95]);
MakeViolinPlot_DB(Pl,Cols,1:5,[],0);
ax = gca;
ax.YAxis.Exponent = 0;
hold on
for i=1:2
    if ~isempty(Pl{i})
        handlesplot=plotSpread(Pl{i}(idx_SZ),'distributionColors',[1 0 0],'xValues',i,'spreadWidth',0.8);
        set(handlesplot{1},'MarkerSize',15);
    end
end
hold off
set(gca,'XTick', [1:5], 'XTickLabel', {'PreREM','PostREM', '', 'PreREM', 'PostREM'}, 'FontSize', 19, 'FontWeight', 'bold');
yl = ylim;
%//
% Add groups - code stolen from <stackoverflow.com/questions/33165830/double-ticklabel-in-matlab>
groupX = [1.5 4.5]; %// central value of each group
groupY = yl(1) - yl(2)*0.08; %// vertical position of texts. Adjust as needed
deltaY = .03; %// controls vertical compression of axis. Adjust as needed
groupNames = {'PlaceCells', 'non Place Cells'};
for g = 1:numel(groupX)
    h = text(groupX(g), groupY, groupNames{g}, 'Fontsize', 23, 'Fontweight', 'bold');
    %// create text for group with appropriate font size and weight
    pos = get(h, 'Position');
    ext = get(h, 'Extent');
    pos(1) = pos(1) - ext(3)/2; %// horizontally correct position to make it centered
    set(h, 'Position', pos); %// set corrected position for text
end
pos = get(gca, 'position');
pos(2) = pos(2) + deltaY; %// vertically compress axis to make room for texts
set(gca, 'Position', pos); %/ set corrected position for axis
%//
ylabel('Number of spikes')
title(['Number of spikes during second hour of NREM sleep, N(PCs)=' num2str(numPCs) ', N(nonPCs)=' num2str(length(FR.Pre.REM.nonPCs))],...
    'FontSize',16);
if savfig
    saveas(f_nrem,[dropbox pathfig 'numSp_NREM_lastHour.fig']);
    saveFigure(f_nrem,'numSp_NREM_lastHour',[dropbox pathfig]);
end

% Number of stimulations in second hour of NREM
Pl = {numStim.Pre.NREM_last.PCs, numStim.Pre.NREM_last.PCs};
Cols = {[0.7 0.7 0.9], [0.2 0.2 0.9]};
tits = {['PreNREM 2nd hour: all PCs (N=' num2str(sum(~isnan(numStim.Pre.NREM_last.PCs(:,1)))) ')'],...
    ['PostNREM 2nd hour: all PCs (N=' num2str(sum(~isnan(numStim.Post.NREM_last.PCs(:,1)))) ')']};
f3 = figure('units', 'normalized', 'outerposition', [0 0 0.27 0.75]);
ax = arrayfun(@(i) subplot(2,1,i, 'NextPlot', 'add', 'Box', 'off'), [1:2]);
for i = 1:length(ax)
    axes(ax(i));
    shadedErrorBar(ISI,nanmean(Pl{i},1),nanstd(Pl{i},1),...
        {'-ok','markerfacecolor',[0.7,0.4,0.3], 'Color', Cols{i}});
    xlim([ISI(1) ISI(end)]);
    set(gca, 'FontSize', 16, 'FontWeight', 'bold','XTick', ISI, 'XTickLabel',...
        {'100 ms', '500 ms', '1 s', '2s', '3s', '5s'});
    set(gca,'xscale','log');
    xlabel('log ISI');
    ylabel('# Stim');
    title(tits{i}, 'FontSize', 14);
end
yl = ylim;
if savfig
    saveas(f3,[dropbox pathfig 'numStim_NREM_lastH.fig']);
    saveFigure(f3,'numStim_NREM_lastH',[dropbox pathfig]);
end
% Number of stimulations in REM - separate figure for only SZ cells
Pl = {numStim.Pre.NREM_last.PCs(idx_SZ,:), numStim.Pre.NREM_last.PCs(idx_SZ,:)};
Cols = {[0.9 0.7 0.7], [0.9 0.2 0.2]};
tits = {['PreNREM 2nd hour: shock zone PCs (N=' num2str(sum(~isnan(numStim.Pre.NREM_last.PCs(idx_SZ,1)))) ')'],...
    ['PostNREM 2nd hour: shock zone PCs (N=' num2str(sum(~isnan(numStim.Post.NREM_last.PCs(idx_SZ,1)))) ')']};
f4 = figure('units', 'normalized', 'outerposition', [0 0 0.27 0.75]);
ax = arrayfun(@(i) subplot(2,1,i, 'NextPlot', 'add', 'Box', 'off'), [1:2]);
for i = 1:length(ax)
    axes(ax(i));
    shadedErrorBar(ISI,nanmean(Pl{i},1),nanstd(Pl{i},1),...
        {'-ok','markerfacecolor',[0.7,0.4,0.3], 'Color', Cols{i}});
    xlim([ISI(1) ISI(end)]);
    set(gca, 'FontSize', 16, 'FontWeight', 'bold','XTick', ISI, 'XTickLabel',...
        {'100 ms', '500 ms', '1 s', '2s', '3s', '5s'});
    set(gca,'xscale','log');
    xlabel('log ISI');
    ylabel('# Stim');
    title(tits{i}, 'FontSize', 14);
    ylim(yl);
end
if savfig
    saveas(f4,[dropbox pathfig 'numStim_NREM_lastH_SZ.fig']);
    saveFigure(f4,'numStim_NREM_lastH_SZ',[dropbox pathfig]);
end

%%%%%%%%%%% NREM SLEEP %%%%%%%%%%%%%%%%%%%%%%%%%