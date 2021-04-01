function [rho_overlap, rho_distant] = PairwiseCorrelationsAcrossMice_WM1994_DB(nmice, experiment, varargin)
% PairwiseCorrelationsAcrossMice_WM1994_DB
%
% This function performs pairwise correlations on all overlapping
% and non-overlapping place cells of each particular mouse and plots the
% results
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
%   RepeatOriginalPaper if true, tt repeats the analysis of Wilson&McNaughton, 1994, Science
%                       (bins spike trains with 1ms bin and calculate average cross-correlation
%                       100 ms around). If false, calculates Pearson's 
%                       correlation coefficient on binned spiked traines
%                       (default = false) (optional)
%   RestrictSleep       if not empty, restrict sleep sessions to this time
%                       (in tsd units) (default = 20*60*1e4) (optional)
%   SaveFigure          if true, saves figures in dropbox (default = false)(optional)
%   PlotExamples        if true, plots cross-correlation histograms for
%                       every pair (optional)
% 
% OUTPUT
% 
%  None
% 
% EXAMPLE
% 
%   PairwiseCorrelationsAcrossMice_WM1994_DB([797 798 1117], 'PAG');
%   PairwiseCorrelationsAcrossMice_WM1994_DB([797 798 1117], 'MFB', 'SaveFigure', true, 'BinSize', 0.05*1e4);
%
% By Dima Bryzgalov, MOBS team, Paris,
% 01/05/2020; reworked on 05/03/2021
% github.com/bryzgalovdm

%% Parameters
% Behavioral states
if strcmp(experiment, 'PAG')
    states = {'PreSleep', 'Task', 'CondMov', 'CondFreeze', 'PostSleep', 'PostTests'};
elseif strcmp(experiment, 'MFB')
    states = {'PreSleep', 'Task', 'CondMov', 'PostSleep', 'PostTests'};
end
% Colors for plots
cols = {[0.2 0.2 0.2], [0.8 0.8 0.8]};

%% Optional parameters

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
savfig = false;
% Do you want to plot each single cross correlation of overlapping pairs?
plotexamples = false;
% Do you want to restrict your sleep time? If not put []
SleepTimeToRestrict = 20*60*1e4;

% Parse parameters
for i=1:2:length(varargin)
    switch(lower(varargin{i}))
        case 'speedthresh'
            speed_thresh = varargin{i+1};
            if ~isa(speed_thresh, 'numeric') && speed_thresh <= 0
                error('Incorrect value for property ''SpeedThresh'' (type ''help PairwiseCorrelationsAcrossMice_WM1994_DB'' for details).');
            end
         case 'binsize'
            binsize = varargin{i+1};
            if ~isa(binsize, 'numeric')
                error('Incorrect value for property ''BinSize'' (type ''help PairwiseCorrelationsAcrossMice_WM1994_DB'' for details).');
            end
         case 'overlapfactor'
            overlapFactor = varargin{i+1};
            if ~isa(overlapFactor, 'numeric')
                error('Incorrect value for property ''overlapFactor'' (type ''help PairwiseCorrelationsAcrossMice_WM1994_DB'' for details).');
            end
         case 'savefigure'
            savfig = varargin{i+1};
            if savfig ~= 1 && savfig ~= 0
                error('Incorrect value for property ''SaveFigure'' (type ''help PairwiseCorrelationsAcrossMice_WM1994_DB'' for details).');
            end
         case 'plotexamples'
            plotexamples = varargin{i+1};
            if plotexamples ~= 1 && plotexamples ~= 0
                error('Incorrect value for property ''PlotExamples'' (type ''help PairwiseCorrelationsAcrossMice_WM1994_DB'' for details).');
            end
         case 'restrictsleep'
            SleepTimeToRestrict = varargin{i+1};
            if ~isa(SleepTimeToRestrict, 'numeric')
                error('Incorrect value for property ''RestrictSleep'' (type ''help PairwiseCorrelationsAcrossMice_WM1994_DB'' for details).');
            end
         case 'repeatoriginalpaper'
            do_as_in_paper = varargin{i+1};
            if do_as_in_paper ~= 1 && do_as_in_paper ~= 0
                error('Incorrect value for property ''RepeatOriginalPaper'' (type ''help PairwiseCorrelationsAcrossMice_WM1994_DB'' for details).');
            end
    end
end

%% Manage experiment
if strcmp(experiment, 'PAG')
    fetchpaths = 'UMazePAG';
elseif strcmp(experiment, 'MFB')
    fetchpaths = 'StimMFBWake';
end

%% Load Data
% Get paths
Dir = PathForExperimentsERC(fetchpaths);
Dir = RestrictPathForExperiment(Dir,'nMice',nmice);

% Allocate memory
spikes = cell(1,length(Dir.path));
behav = cell(1,length(Dir.path));
sleep = cell(1,length(Dir.path));
Epochs = cell(1, length(Dir.path));
UMazeEpoch = cell(1, length(Dir.path));
overlappairs = cell(1, length(Dir.path));
distantpairs = cell(1, length(Dir.path));
S_PC = cell(1, length(Dir.path));
QEpoched = cell(1, length(Dir.path));
rho_overlap = cell(1, length(states));
rho_distant = cell(1, length(states));
for iepoch = 1:length(states)
    rho_overlap{iepoch} = nan(1, 1e5);
    rho_distant{iepoch} = nan(1, 1e5);
end

for i=1:length(Dir.path)
    spikes{i} = load([Dir.path{i}{1} 'SpikeData.mat'],'S','PlaceCells');
    behav{i} = load([Dir.path{i}{1} 'behavResources.mat'],'SessionEpoch','CleanVtsd','CleanXtsd','CleanYtsd','FreezeAccEpoch');
    try
        sleep{i} = load([Dir.path{i}{1} 'SleepScoring_OBGamma.mat'],'SWSEpoch','REMEpoch','Sleep');
    catch
        sleep{i} = load([Dir.path{i}{1} 'SleepScoring_Accelero.mat'],'SWSEpoch','REMEpoch','Sleep'); % REM and Sleep are not used
    end
end

%% Create Epochs
for i=1:length(Dir.path)
    
    [~, UMazeEpoch{i}, ConditioningEpoch, TaskEpoch, AfterConditioningEpoch] = ReturnMnemozyneEpochs(behav{i}.SessionEpoch,...
        'Speed', behav{i}.CleanVtsd, 'SpeedThresh', speed_thresh);
    
    % Get freezing epoch during cond
    if strcmp(experiment, 'PAG')
        [~,~,tempEpoch] = ReturnMnemozyneEpochs(behav{i}.SessionEpoch);
        ConditioningFreezingEpoch = and(behav{i}.FreezeAccEpoch, tempEpoch);
    end
    
    % SleepEpochs
    if isempty(SleepTimeToRestrict)
        PreSleepFinal = and(behav{i}.SessionEpoch.PreSleep, sleep{i}.SWSEpoch);
        PostSleepFinal = and(behav{i}.SessionEpoch.PostSleep, sleep{i}.SWSEpoch);
    else
        PreSleepFinal = RestrictToTime(and(behav{i}.SessionEpoch.PreSleep, sleep{i}.SWSEpoch),SleepTimeToRestrict);
        PostSleepFinal = RestrictToTime(and(behav{i}.SessionEpoch.PostSleep, sleep{i}.SWSEpoch),SleepTimeToRestrict);
    end
    
    % Pack everything in one cell
    if strcmp(experiment, 'PAG')
        Epochs{i} = {PreSleepFinal, TaskEpoch, ConditioningEpoch, ConditioningFreezingEpoch, PostSleepFinal, AfterConditioningEpoch};
    elseif strcmp(experiment, 'MFB')
        Epochs{i} = {PreSleepFinal, TaskEpoch, ConditioningEpoch, PostSleepFinal, AfterConditioningEpoch};
    end
end
 
%% Create binned FR vectors
for i=1:length(Dir.path)
    QEpoched{i} = cell(1, length(states));
    
    if do_as_in_paper
        Q = MakeQfromS(spikes{i}.S, 10); % 1 ms
        
        for iepoch = 1:length(states)
            QEpoched{i}{iepoch} = full(Data(Restrict(Q, Epochs{i}{iepoch})));
        end
        clear Q
    else
        Q = MakeQfromS(spikes{i}.S, binsize);
        
        for iepoch = 1:length(states)
            QEpoched{i}{iepoch} = zscore(full(Data(Restrict(Q, Epochs{i}{iepoch}))));
        end
        clear Q
    end
end

%% Calculate overlap for each place field
for i=1:length(Dir.path)
    % Get the stats for each cell
    if isfield(spikes{i}.PlaceCells,'idx')
        if length(spikes{i}.PlaceCells.idx) > PCnum_thresh %%% Take only mice with number of place cells > PCnum_thresh
            
            S_PC{i} = spikes{i}.S(spikes{i}.PlaceCells.idx);
            
            [overlappairs{i}, distantpairs{i}] = FindOverlappingPlaceFields(S_PC{i}, behav{i}.CleanXtsd, behav{i}.CleanYtsd,...
                UMazeEpoch{i}, overlapFactor);
            
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
    for imouse = 1:length(Dir.path)
        PlotCorrelationExamples(Dir.name{imouse}, spikes{imouse}, S_PC{imouse}, overlappairs{imouse},...
            Epochs{imouse}{1}, Epochs{imouse}{2}, Epochs{imouse}{5});
    end
end

%% Calculate cross-correlation between spike trains
cnt_o = 1;
cnt_d = 1;
for i=1:length(Dir.path)

    if do_as_in_paper
        % Perform calculations - overlapping cells
        if ~isempty(overlappairs{i})
            for j=1:length(overlappairs{i})
                pair = overlappairs{i}{j};
                for iepoch = 1:length(states)
                    rho_overlap{iepoch}(cnt_o) = mean(crosscorr(QEpoched{i}{iepoch}(:,spikes{i}.PlaceCells.idx(pair(1))),...
                        QEpoched{i}{iepoch}(:,spikes{i}.PlaceCells.idx(pair(2))), 'NumLags', 100));
                end
                cnt_o=cnt_o+1;
            end
        end
        clear pair
        % Perform calculations - non-overlapping cells
        if ~isempty(distantpairs{i})
            for j=1:length(distantpairs{i})
                pair = distantpairs{i}{j};
                for iepoch = 1:length(states)
                    rho_distant{iepoch}(cnt_d) = mean(crosscorr(QEpoched{i}{iepoch}(:,spikes{i}.PlaceCells.idx(pair(1))),...
                        QEpoched{i}{iepoch}(:,spikes{i}.PlaceCells.idx(pair(2))), 'NumLags', 100));
                end
                cnt_d=cnt_d+1;
            end
        end
        
    else % do_as_in_paper
        % Perform calculations - overlapping cells
        if ~isempty(overlappairs{i})
            for j=1:length(overlappairs{i})
                pair = overlappairs{i}{j};
                for iepoch = 1:length(states)
                    temp = corrcoef(QEpoched{i}{iepoch}(:,spikes{i}.PlaceCells.idx(pair(1))), ...
                        QEpoched{i}{iepoch}(:,spikes{i}.PlaceCells.idx(pair(2))));
                        rho_overlap{iepoch}(cnt_o) = temp(1,2); clear temp
                end
                cnt_o = cnt_o + 1;
            end
        end
        clear pair
        
        % Perform calculations - non-overlapping cells
        if ~isempty(distantpairs{i})
            for j=1:length(distantpairs{i})
                pair = distantpairs{i}{j};
                for iepoch = 1:length(states)
                    temp = corrcoef(QEpoched{i}{iepoch}(:,spikes{i}.PlaceCells.idx(pair(1))), ...
                        QEpoched{i}{iepoch}(:,spikes{i}.PlaceCells.idx(pair(2))));
                        rho_distant{iepoch}(cnt_d) = temp(1,2); clear temp
                end
                
                cnt_d = cnt_d + 1;
            end
        end
    end
end
% Remove missing values (PostSleep is a reference)
id_postsleep = find(contains(states, 'PostSleep'));
idx_toremove_overlap = find(isnan(rho_overlap{id_postsleep}),1);
idx_toremove_distant = find(isnan(rho_distant{id_postsleep}),1);
for iepoch = 1:length(states)
    rho_overlap{iepoch}(idx_toremove_overlap:end) = [];
    rho_distant{iepoch}(idx_toremove_distant:end) = [];
end

id_presleep = find(contains(states, 'PreSleep'));
id_task = find(contains(states, 'Task'));
id_condmov = find(contains(states, 'CondMov'));
id_condfre = find(contains(states, 'CondFreeze'));
id_post = find(contains(states, 'PostTests'));

%% Figures
% Short figure 1
foldertosave = ChooseFolderForFigures_DB('ReactReplay');

f1 = figure('units', 'normalized', 'outerposition', [0 0 0.4 0.7]);
cols_ssd = cols;
while length(cols_ssd) < 6
    cols_ssd{length(cols_ssd)+1} = cols{1};
    cols_ssd{length(cols_ssd)+1} = cols{2};
end

short_stacked_data = {rho_overlap{id_presleep}, rho_distant{id_presleep},...
    rho_overlap{id_task}, rho_distant{id_task}, rho_overlap{id_postsleep}, rho_distant{id_postsleep}};
ssp = MakeBoxPlot_DB(short_stacked_data, cols_ssd, [1 2 4 5 7 8], [], 0);
set(gca,'Xtick',1.5:3:7.5,'XtickLabel',{'PreSleep', 'Task', 'PostSleep'})
ylabel('Corellation')
title('Pairwise correlations');
ylim([-.1 .3])
legend([ssp{1}.handles.box, ssp{2}.handles.box], {'Overlapping PCs', 'Non-overlapping PCs'});
makepretty
% Save the figure
if savfig
    if exist([foldertosave filesep 'Pairwise'], 'dir') ~= 7
        mkdir([foldertosave filesep 'Pairwise']);
    end
    saveas(f1,[foldertosave filesep 'Pairwise' filesep 'Pairwise_Small_' experiment '.fig']);
    saveFigure(f1,['Pairwise_Small_' experiment], [foldertosave filesep 'Pairwise']);
end

% Final figures for Karim
% Reorganize data a bit (PreSleep-PostSleep-Maze-CondMov-CondFreeze-PostTests)
if strcmp(experiment, 'PAG')
    largefig_data{1} = {rho_overlap{id_presleep}, rho_overlap{id_postsleep}, rho_overlap{id_task},...
        rho_overlap{id_condmov}, rho_overlap{id_condfre}, rho_overlap{id_post}};
    largefig_data{2} = {rho_distant{id_presleep}, rho_distant{id_postsleep}, rho_distant{id_task},...
        rho_distant{id_condmov}, rho_distant{id_condfre}, rho_distant{id_post}};
elseif strcmp(experiment, 'MFB')
    largefig_data{1} = {rho_overlap{id_presleep}, rho_overlap{id_postsleep}, rho_overlap{id_task},...
        rho_overlap{id_condmov}, rho_overlap{id_post}};
    largefig_data{2} = {rho_distant{id_presleep}, rho_distant{id_postsleep}, rho_distant{id_task},...
        rho_distant{id_condmov}, rho_distant{id_post}};
end
pairs_toCompare = {[1 2], [1 3], [2 3]};

% Big figure
f2 = figure('units', 'normalized', 'outerposition', [0.1 0.1 0.83 0.7]);
ax{1} = axes('position', [0.05 0.05 0.43 0.9]);
ax{2} = axes('position', [0.54 0.05 0.43 0.9]);
titles = {'Overlapping place cells', 'Non-Overlapping place cells'};
for icolor = 1:length(cols)
    cols_lf{icolor} = repelem(cols(icolor), length(states));
end
% Plot
for iaxis = 1:length(ax)
    axes(ax{iaxis});
    MakeBoxPlot_DB(largefig_data{iaxis}, cols_lf{iaxis}, 1:length(states), [], 0);
    if strcmp(experiment, 'PAG')
        set(gca,'Xtick', 1:length(states), 'XtickLabel',{'PreSleep', 'PostSleep', 'Task', 'StimRun',...
            'StimFreeze', 'PostTest'})
    elseif strcmp(experiment, 'MFB')
        set(gca,'Xtick', 1:length(states), 'XtickLabel',{'PreSleep', 'PostSleep', 'Task', 'StimRun', 'PostTest'})
    end
    ylabel('Correlation', 'FontSize', 18)
    title(titles{iaxis});
    p = DoWilcoxonOnArray(largefig_data{iaxis}, pairs_toCompare);
    for ip = 1:length(p)
        if p(ip) <= 0.05
            sigstar_DB(pairs_toCompare(ip),p(ip),0,'LineWigth',16,'StarSize',24);
        end
    end
    makepretty
end
ylims_out = SelectYlim(f2);
for iaxis = 1:length(ax)
    axes(ax{iaxis});
    ylim(ylims_out);
end
if savfig
    if exist([foldertosave filesep 'Pairwise'], 'dir') ~= 7
        mkdir([foldertosave filesep 'Pairwise']);
    end
    saveas(f2,[foldertosave filesep 'Pairwise' filesep 'Pairwise_final_big_' experiment '.fig']);
    saveFigure(f2,['Pairwise_final_big_' experiment], [foldertosave filesep 'Pairwise']);
end

% Short figure 2
shortfig_data{1} = {rho_overlap{id_presleep}, rho_overlap{id_task}, rho_overlap{id_postsleep}};
shortfig_data{2} = {rho_distant{id_presleep}, rho_distant{id_task}, rho_distant{id_postsleep}};

f3 = figure('units', 'normalized', 'outerposition', [0.1 0.15 0.5 0.7]);
ax{1} = axes('position', [0.05 0.05 0.43 0.9]);
ax{2} = axes('position', [0.54 0.05 0.43 0.9]);
for icolor = 1:length(cols)
    cols_sf{icolor} = repelem(cols(icolor), 3);
end
for iaxis = 1:length(ax)
    axes(ax{iaxis});
    MakeBoxPlot_DB(shortfig_data{iaxis}, cols_sf{iaxis}, 1:3, [], 0);
    set(gca,'Xtick', 1:3, 'XtickLabel',{'PreSleep', 'Task', 'PostSleep'})
    ylabel('Correlation', 'FontSize', 18)
    title(titles{iaxis});
    p = DoWilcoxonOnArray(shortfig_data{iaxis}, pairs_toCompare);
    for ip = 1:length(p)
        if p(ip) <= 0.05
            sigstar_DB(pairs_toCompare(ip),p(ip),0,'LineWigth',16,'StarSize',24);
        end
    end
    makepretty
end
ylims_out = SelectYlim(f3);
for iaxis = 1:length(ax)
    axes(ax{iaxis});
    ylim(ylims_out);
end
if savfig
    if exist([foldertosave filesep 'Pairwise'], 'dir') ~= 7
        mkdir([foldertosave filesep 'Pairwise']);
    end
    saveas(f3,[foldertosave filesep 'Pairwise' filesep 'Pairwise_final_small_' experiment '.fig']);
    saveFigure(f3,['Pairwise_final_small_' experiment], [foldertosave filesep 'Pairwise']);
end

% Only overlapping cells
f4 = figure('units', 'normalized', 'outerposition', [0.3 0.2 0.285 0.7]);
MakeBoxPlot_DB(shortfig_data{1}, cols_sf{1}, 1:3, [], 0);
set(gca,'Xtick', 1:3, 'XtickLabel',{'PreSleep', 'Task', 'PostSleep'})
ylabel('Correlation', 'FontSize', 18)
p = DoWilcoxonOnArray(shortfig_data{1}, pairs_toCompare);
for ip = 1:length(p)
    if p(ip) <= 0.05
        sigstar_DB(pairs_toCompare(ip),p(ip),0,'LineWigth',16,'StarSize',24);
    end
end
makepretty
if savfig
    if exist([foldertosave filesep 'Pairwise'], 'dir') ~= 7
        mkdir([foldertosave filesep 'Pairwise']);
    end
    saveas(f4,[foldertosave filesep 'Pairwise' filesep 'Pairwise_fullfinal_pooled_' experiment '.fig']);
    saveFigure(f4,['Pairwise_fullfinal_pooled_'  experiment], [foldertosave filesep 'Pairwise']);
end


end







%% Auxiliary functions

function PlotCorrelationExamples(MouseName, spikes, S_PC, overlappairs, PreSleepEpoch, TaskEpoch, PostSleepEpoch, varargin)
    
% Optional parameters
binsize = 4;
windowsize = 400;
% Parse parameters
for k=1:2:length(varargin)
    switch(lower(varargin{k}))
         case 'binsize'
            binsize = varargin{k+1};
            if ~isa(binsize, 'numeric')
                error('Incorrect value for property ''BinSize'' (type ''help PlotCorrelationExamples'' for details).');
            end
         case 'windowsize'
            windowsize = varargin{k+1};
            if ~isa(windowsize, 'numeric')
                error('Incorrect value for property ''WindowSize'' (type ''help PlotCorrelationExamples'' for details).');
            end
    end
end

nbins = windowsize*2/binsize;

% Do the job and plot the figure

if ~isempty(overlappairs)
    CPre = cell(length(overlappairs),1);
    CTask = cell(length(overlappairs),1);
    CPost = cell(length(overlappairs),1);
    for ipair=1:length(overlappairs)
        [CPre{ipair},B] = CrossCorrDB(Data(Restrict(S_PC{overlappairs{ipair}(1)},PreSleepEpoch)),...
            Data(Restrict(S_PC{overlappairs{ipair}(2)},PreSleepEpoch)),binsize,nbins);
        CTask{ipair} = CrossCorrDB(Data(Restrict(S_PC{overlappairs{ipair}(1)},TaskEpoch)),...
            Data(Restrict(S_PC{overlappairs{ipair}(2)},TaskEpoch)),binsize,nbins);
        CPost{ipair} = CrossCorrDB(Data(Restrict(S_PC{overlappairs{ipair}(1)},PostSleepEpoch)),...
            Data(Restrict(S_PC{overlappairs{ipair}(2)},PostSleepEpoch)),binsize,nbins);
    end
    
    % Figure
    if length(overlappairs) < 15
        fff = figure('units', 'normalized', 'outerposition', [0 0.55 0.7 0.4], 'Name', MouseName);
        title(MouseName);
        tabs = arrayfun(@(x) uitab('Title', ['Cl' num2str(spikes.PlaceCells.idx(overlappairs{x}(1)))...
            ' vs Cl' num2str(spikes.PlaceCells.idx(overlappairs{x}(2)))]), 1:length(overlappairs));
        ax = arrayfun(@(tab) axes(tab), tabs);
        arrayfun(@(k) title(ax(k), ['Cl' num2str(spikes.PlaceCells.idx(overlappairs{k}(1)))...
            ' vs Cl' num2str(spikes.PlaceCells.idx(overlappairs{k}(2)))], 'FontSize', 16), 1:length(overlappairs));
        arrayfun(@(k) set(ax(k), 'NextPlot','add', 'FontSize', 16, 'FontWeight', 'bold'), 1:length(overlappairs));
        arrayfun(@(k) plot(ax(k), B, CPre{k}, 'k', 'LineWidth', 2),...
            1:length(overlappairs));
        arrayfun(@(k) plot(ax(k), B, CTask{k}, 'm', 'LineWidth', 2),...
            1:length(overlappairs));
        arrayfun(@(k) plot(ax(k), B, CPost{k}, 'b', 'LineWidth', 2),...
            1:length(overlappairs));
        arrayfun(@(k) line(ax(k),[0 0], ylim, 'Color', 'r', 'LineWidth',3), 1:length(overlappairs));
        arrayfun(@(k) ylabel(ax(k), 'Count'), 1:length(overlappairs));
        arrayfun(@(k) xlabel(ax(k), 'Time (ms)'), 1:length(overlappairs));
        arrayfun(@(k) legend(ax(k), {'Pre', 'Task', 'Post'}), 1:length(overlappairs));
    else
        for numfig = 1:ceil(length(overlappairs)/15)
            try
                rr = ((numfig-1)*15)+1:((numfig-1)*15)+15;
                temp_o = overlappairs(rr);
            catch
                rr = ((numfig-1)*15)+1:length(overlappairs);
                temp_o = overlappairs(rr);
            end
            fff1(numfig) = figure('units', 'normalized', 'outerposition', [0 0.55 0.7 0.4], 'Name', MouseName);
            title(MouseName);
            tabs = arrayfun(@(x) uitab('Title', ['Cl' num2str(spikes.PlaceCells.idx(temp_o{x}(1)))...
                ' vs Cl' num2str(spikes.PlaceCells.idx(temp_o{x}(2)))]), 1:length(temp_o));
            ax = arrayfun(@(tab) axes(tab), tabs);
            arrayfun(@(k) title(ax(k), ['Cl' num2str(spikes.PlaceCells.idx(temp_o{k}(1)))...
                ' vs Cl' num2str(spikes.PlaceCells.idx(temp_o{k}(2)))], 'FontSize', 16), 1:length(temp_o));
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