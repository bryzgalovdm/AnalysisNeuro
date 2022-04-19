function FRduringSleep(mice, experiment, varargin)

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
SleepTimeToRestrict = 2*60*60*1e4; % Sleep time to restrict (in tsd units)
savfig = false; % Do you want to save the figures?

%% Parse parameters
for i=1:2:length(varargin)
    switch(lower(varargin{i}))
        case 'sleeprestrict'
            SleepTimeToRestrict = varargin{i+1};
            if ~isa(SleepTimeToRestrict, 'numeric') && SleepTimeToRestrict <= 0
                error('Incorrect value for property ''SleepRestrict'' (type ''help ExplainedVariance_master_DB'' for details).');
            end
        case 'issavefig'
            savfig = varargin{i+1};
            if savfig ~= 1 && savfig ~= 0
                error('Incorrect value for property ''IsSaveFig'' (type ''help ExplainedVariance_master_DB'' for details).');
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

%% Get folders
Dir = PathForExperimentsERC(fetchpaths);
Dir = RestrictPathForExperiment(Dir,'nMice',mice);
numsessions = CountNumSesionsERC(Dir);

%% Preallocation
% Data
spikes = cell(numsessions, 1);
behav = cell(numsessions, 1);
sleep = cell(numsessions, 1);
% Neurons sets
PCs = cell(numsessions, 1);
nonPCs = cell(numsessions, 1);
% Epochs
NREMPre = cell(numsessions, 1);
REMPre = cell(numsessions, 1);
NREMPost = cell(numsessions, 1);
REMPost = cell(numsessions, 1);
% Results - FR
FRNREMPre = cell(numsessions, 1);
FRREMPre = cell(numsessions, 1);
FRNREMPost = cell(numsessions, 1);
FRREMPost = cell(numsessions, 1);

%% Load data - here I handle the exceptions too

cnt=1;
for imouse=1:length(Dir.path)
    for isession = 1:length(Dir.path{imouse})
        spikes{cnt} = load([Dir.path{imouse}{isession} 'SpikeData.mat'],'S','PlaceCells', 'BasicNeuronInfo');
        behav{cnt} = load([Dir.path{imouse}{isession} 'behavResources.mat'],'SessionEpoch','Vtsd','AlignedXtsd','AlignedYtsd','FreezeAccEpoch');
        try
            sleep{cnt} = load([Dir.path{imouse}{isession} 'SleepScoring_OBGamma.mat'],'SWSEpoch','REMEpoch');
        catch
            sleep{cnt} = load([Dir.path{imouse}{isession} 'SleepScoring_Accelero.mat'],'SWSEpoch','REMEpoch');
        end
        cnt=cnt+1;
    end
end

%% Split spike arrays into PCs and nonPCs
for i = 1:numsessions
    if ~isempty(spikes{i}.PlaceCells.idx)
        PCs{i} = spikes{i}.PlaceCells.idx;
        nonPCs{i} = setdiff(find(spikes{i}.BasicNeuronInfo.neuroclass > 0),spikes{i}.BasicNeuronInfo.idx_MUA);
    end
end

%% Create epochs

for i=1:numsessions
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


%% Calculate firing rate in PCs

for i = 1:numsessions
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
        

        
        % PlaceCells
        for j = 1:length(PCs{i})
            FRNREMPre{i}.PCs(j) = length(Restrict(spikes{i}.S{PCs{i}(j)}, NREMPre{i}))/sum(End(NREMPre{i}, 's') - Start(NREMPre{i}, 's'));
            FRREMPre{i}.PCs(j) = length(Restrict(spikes{i}.S{PCs{i}(j)}, REMPre{i}))/sum(End(REMPre{i}, 's') - Start(REMPre{i}, 's'));
            FRNREMPost{i}.PCs(j) = length(Restrict(spikes{i}.S{PCs{i}(j)}, NREMPost{i}))/sum(End(NREMPost{i}, 's') - Start(NREMPost{i}, 's'));
            FRREMPost{i}.PCs(j) = length(Restrict(spikes{i}.S{PCs{i}(j)}, REMPost{i}))/sum(End(REMPost{i}, 's') - Start(REMPost{i}, 's'));

        end
        
        % non-PlaceCells
        for j = 1:length(nonPCs{i})
            FRNREMPre{i}.nonPCs(j) = length(Restrict(spikes{i}.S{nonPCs{i}(j)}, NREMPre{i}))/sum(End(NREMPre{i}, 's') - Start(NREMPre{i}, 's'));
            FRREMPre{i}.nonPCs(j) = length(Restrict(spikes{i}.S{nonPCs{i}(j)}, REMPre{i}))/sum(End(REMPre{i}, 's') - Start(REMPre{i}, 's'));
            FRNREMPost{i}.nonPCs(j) = length(Restrict(spikes{i}.S{nonPCs{i}(j)}, NREMPost{i}))/sum(End(NREMPost{i}, 's') - Start(NREMPost{i}, 's'));
            FRREMPost{i}.nonPCs(j) = length(Restrict(spikes{i}.S{nonPCs{i}(j)}, REMPost{i}))/sum(End(REMPost{i}, 's') - Start(REMPost{i}, 's'));
        end
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


if numsessions > 1
    for i = 2:numsessions
        if ~isempty(PCs{i})
            FR.Pre.NREM.PCs = [FR.Pre.NREM.PCs; FRNREMPre{i}.PCs];
            FR.Pre.NREM.nonPCs = [FR.Pre.NREM.nonPCs; FRNREMPre{i}.nonPCs];
            FR.Pre.REM.PCs = [FR.Pre.REM.PCs; FRREMPre{i}.PCs];
            FR.Pre.REM.nonPCs = [FR.Pre.REM.nonPCs; FRREMPre{i}.nonPCs];
            
            FR.Post.NREM.PCs = [FR.Post.NREM.PCs; FRNREMPost{i}.PCs];
            FR.Post.NREM.nonPCs = [FR.Post.NREM.nonPCs; FRNREMPost{i}.nonPCs];
            FR.Post.REM.PCs = [FR.Post.REM.PCs; FRREMPost{i}.PCs];
            FR.Post.REM.nonPCs = [FR.Post.REM.nonPCs; FRREMPost{i}.nonPCs];
        end
        
    end
end

%% Get shock zone PCs idxs for plotting
[numPCs,numPCs_M] = CountPlaceCells(Dir, 'Verbose', false);
idx_SZ = (arrayfun(@(x)find(spikes{1}.PlaceCells.idx==x,1),spikes{1}.PlaceCells.SZ))';
for i=2:numsessions
    idx_SZ = [idx_SZ; (arrayfun(@(x)find(spikes{i}.PlaceCells.idx==x,1),spikes{i}.PlaceCells.SZ)+sum(numPCs_M(1:i-1)))'];
end

%% Illustrate the data

%%%%%%%%%%% REM SLEEP %%%%%%%%%%%%%%%%%%%%%%%%%
% FR during REM
Pl_PC = {FR.Pre.REM.PCs; FR.Post.REM.PCs};
Pl_nPC = {FR.Pre.REM.nonPCs; FR.Post.REM.nonPCs};
% Cols = {[.9856, .7372, .2537], [1 0 1]};
Cols = {[.9856, .7372, .2537], [.16 .95 1]};
% Place Cells
ffr1 = figure('units', 'normalized', 'outerposition', [0 0 0.45 0.8]);
h = MakeViolinPlot_DB(Pl_PC,Cols,1:2,[],0);
ylim([0 5]);
for k = 1:length(h)
    h{k}.violinAlpha = 0.7;
end
hold on
for i=1:2
    if ~isempty(Pl_PC{i})
        handlesplot=plotSpread(Pl_PC{i}(idx_SZ),'distributionColors',[1 0 0],'xValues',i,'spreadWidth',0.8);
        set(handlesplot{1},'MarkerSize',15);
    end
end
hold off
set(gca,'XTick', [1:2], 'XTickLabel', {'PreREM','PostREM'});
yl = ylim;
%//
% Add groups - code stolen from <stackoverflow.com/questions/33165830/double-ticklabel-in-matlab>
groupX = 1.5; %// central value of each group
groupY = yl(1) - yl(2)*0.08; %// vertical position of texts. Adjust as needed
deltaY = .03; %// controls vertical compression of axis. Adjust as needed
groupNames = {'PlaceCells'};
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
title(['Firing rate during REM sleep, N(PCs)=' num2str(numPCs)]);
makepretty_DB
if savfig
    saveas(ffr1,[dropbox pathfig 'PAG_FR_REM_PC.fig']);
    saveFigure(ffr1,'PAG_FR_REM_PC',[dropbox pathfig]);
end
% non Place cells
ffr2 = figure('units', 'normalized', 'outerposition', [0 0 0.45 0.8]);
h = MakeViolinPlot_DB(Pl_nPC,Cols,1:2,[],0);
ylim([0 30]);
for k = 1:length(h)
    h{k}.violinAlpha = 0.7;
end
set(gca,'XTick', [1:2], 'XTickLabel', {'PreREM','PostREM'});
yl = ylim;
%//
% Add groups - code stolen from <stackoverflow.com/questions/33165830/double-ticklabel-in-matlab>
groupX = 1.5; %// central value of each group
groupY = yl(1) - yl(2)*0.08; %// vertical position of texts. Adjust as needed
deltaY = .03; %// controls vertical compression of axis. Adjust as needed
groupNames = {'non PlaceCells'};
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
title(['Firing rate during REM sleep, N(nonPCs)=' num2str(length(FR.Pre.REM.nonPCs))]);
ylim([0 30]);
makepretty_DB
if savfig
    saveas(ffr2,[dropbox pathfig 'PAG_FR_REM_nPC.fig']);
    saveFigure(ffr2,'PAG_FR_REM_nPC',[dropbox pathfig]);
end
%%%%%%%%%%% REM SLEEP %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% NREM SLEEP %%%%%%%%%%%%%%%%%%%%%%%%%
% FR during NREM
Pl_PC_NREM = {FR.Pre.NREM.PCs; FR.Post.NREM.PCs;};
Pl_nPC_NREM = {FR.Pre.NREM.nonPCs; FR.Post.NREM.nonPCs};
% Place Cells
ffr3 = figure('units', 'normalized', 'outerposition', [0 0 0.45 0.8]);
h = MakeViolinPlot_DB(Pl_PC_NREM,Cols,1:2,[],0);
ylim([0 5]);
for k = 1:length(h)
    h{k}.violinAlpha = 0.7;
end
hold on
for i=1:2
    if ~isempty(Pl_PC_NREM{i})
        handlesplot=plotSpread(Pl_PC_NREM{i}(idx_SZ),'distributionColors',[1 0 0],'xValues',i,'spreadWidth',0.8);
        set(handlesplot{1},'MarkerSize',15);
    end
end
hold off
set(gca,'XTick', [1:2], 'XTickLabel', {'PreNREM','PostNREM'});
yl = ylim;
%//
% Add groups - code stolen from <stackoverflow.com/questions/33165830/double-ticklabel-in-matlab>
groupX = 1.5; %// central value of each group
groupY = yl(1) - yl(2)*0.08; %// vertical position of texts. Adjust as needed
deltaY = .03; %// controls vertical compression of axis. Adjust as needed
groupNames = {'PlaceCells'};
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
title(['Firing rate during NREM sleep, N(PCs)=' num2str(numPCs)]);
makepretty_DB
if savfig
    saveas(ffr3,[dropbox pathfig 'PAG_FR_NREM_PC.fig']);
    saveFigure(ffr3,'PAG_FR_NREM_PC',[dropbox pathfig]);
end
% non Place cells
ffr4 = figure('units', 'normalized', 'outerposition', [0 0 0.45 0.8]);
h = MakeViolinPlot_DB(Pl_nPC_NREM,Cols,1:2,[],0);
ylim([0 30]);
for k = 1:length(h)
    h{k}.violinAlpha = 0.7;
end
set(gca,'XTick', [1:2], 'XTickLabel', {'PreNREM','PostNREM'});
yl = ylim;
%//
% Add groups - code stolen from <stackoverflow.com/questions/33165830/double-ticklabel-in-matlab>
groupX = 1.5; %// central value of each group
groupY = yl(1) - yl(2)*0.08; %// vertical position of texts. Adjust as needed
deltaY = .03; %// controls vertical compression of axis. Adjust as needed
groupNames = {'non PlaceCells'};
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
title(['Firing rate during NREM sleep, N(nonPCs)=' num2str(length(FR.Pre.REM.nonPCs))]);
makepretty_DB
if savfig
    saveas(ffr4,[dropbox pathfig 'PAG_FR_NREM_nPC.fig']);
    saveFigure(ffr4,'PAG_FR_NREM_nPC',[dropbox pathfig]);
end
%%%%%%%%%%% NREM SLEEP %%%%%%%%%%%%%%%%%%%%%%%%%

end