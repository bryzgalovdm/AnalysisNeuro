function [EV_4corr, Diff_occup, Diff_entrytime] = EV_Behav_corr_DB_new(mice, expe, EV_epoch, varargin)
%
% The function calculates Spearman rank correlation between explained
% variance and two behavioral parameters:
%   - Shock zone occupancy
%   - Latency to enter the shock zone
%
% INPUT
%
%   mice                array with mice numbers that will harvested from ERC
%                       PathForExperiments. Each mouse should contain
%                       PlaceCells structure in its SpikeData.mat
%   expe                type of experiment to analyse ('PAG' or 'MFB')
%   EV_Epoch            name of the period reactivation of which you want
%                       to assess with the EV (see the full list of options below)
%   SleeType            type of sleep during which you want to assess
%                       reactivations - NREM or REM (optional - default='NREM')
%   IsSave              if true, saves figures in dropbox (default = false)(optional)
%
% OUTPUT
%
%   EV_4corr            array of EV values used for analysis
%   Diff_occup          array of SZ occupancy values used for analysis
%   Diff_entrytime      array of SZ entry latency values used for analysis
%
% LIST OF EPOCHS PERMITTED
%
%   PAG expe           'Explo', 'CondMov', 'CondFreeze', 'FullTask', 'RipplesEpoch', 'PostTests'
%   MFB expe           'Explo', 'CondMov', 'FullTask', 'RipplesEpoch', 'PostTests'
%   Novel expe         'Explo', 'CondMov', 'FullTask', 'RipplesEpoch'
%
% EXAMPLE
%
%   [EV_4corr, Diff_occup, Diff_entrytime] = EV_Behav_corr_DB_new(mice, expe, EV_epoch);
%   [EV_4corr, Diff_occup, Diff_entrytime] = EV_Behav_corr_DB_new(mice, expe, EV_epoch, 'SleepType', 'REM');
%
% By Dima Bryzgalov, MOBS team, Paris,
% 15/09/2021
% github.com/bryzgalovdm

%% Parameters
type_sleep = 'NREM';
sav = false;

%% Optional Arguments
for i=1:2:length(varargin)    
    switch(lower(varargin{i})) 
        case 'issave'
            sav = varargin{i+1};
            if sav ~= 1 && sav ~= 0
                error('Incorrect value for property ''IsSave'' (type ''help EV_Behav_corr_DB_new'' for details).');
            end
        case 'sleptype'
            type_sleep = varargin{i+1};
            if isa(type_sleep, 'char')
                error('Incorrect value for property ''SleepType'' (type ''help EV_Behav_corr_DB_new'' for details).');
            end
    end
end

%% Manage EV_epoch
if strcmp(expe, 'PAG')
    states = {'Explo', 'CondMov', 'CondFreeze', 'FullTask', 'RipplesEpoch', 'PostTests'};
    if sum(strcmp(states, EV_epoch)) < 1
        error('Epoch name does not match any name within EV function. Please check help');
    end
elseif strcmp(expe, 'MFB')
    states = {'Explo', 'CondMov', 'FullTask', 'RipplesEpoch', 'PostTests'};
    if sum(strcmp(states, EV_epoch)) < 1
        error('Epoch name does not match any name within EV function. Please check help');
    end
elseif strcmp(expe, 'Novel')
    states = {'Explo', 'CondMov', 'FullTask', 'RipplesEpoch'};
    if sum(strcmp(states, EV_epoch)) < 1
        error('Epoch name does not match any name within EV function. Please check help');
    end
end

%% Calculate EV
[EV, ~, states_sleep, states_wake, num_mice] = ExplainedVariance_master_DB(mice, expe, 'PlotResults', false, 'IsII', false);

%% Calculate behavioral results
% Manage experiment
if strcmp(expe, 'PAG')
    fetchpaths = 'UMazePAG';
elseif strcmp(expe, 'MFB')
    fetchpaths = 'StimMFBWake';
elseif strcmp(expe, 'Novel')
    fetchpaths = 'Novel';
end

% Get paths of each individual mouse
mice_EV = mice(num_mice);
Dir = PathForExperimentsERC(fetchpaths);
Dir = RestrictPathForExperiment(Dir,'nMice',mice_EV);

% Load the data
cnt=1;
for imouse=1:length(Dir.path)
    for isession = 1:length(Dir.path{imouse})
        a{cnt} = load([Dir.path{imouse}{isession} 'behavResources.mat'], 'behavResources');
        cnt=cnt+1;
    end
end

% Find indices of PreTests and PostTest session in the structure
id_Pre = cell(1,length(a));
id_Post = cell(1,length(a));

for i=1:length(a)
    id_Pre{i} = FindSessionID_ERC(a{i}.behavResources, 'TestPre');
    id_Post{i} = FindSessionID_ERC(a{i}.behavResources, 'TestPost');
end

% Calculate average occupancy
occup_shock = nan(length(Dir.path), 4, 2); % 4 tests, Pre and Post
for i=1:length(a)
    for k=1:length(id_Pre{i})
        temp = CalculateZoneOccupancy(a{i}.behavResources(id_Pre{i}(k)));
        occup_shock(i,k,1) = temp(1);
        temp = CalculateZoneOccupancy(a{i}.behavResources(id_Post{i}(k)));
        occup_shock(i,k,2) = temp(1);
    end
end

occup_shock_mean = nan(length(Dir.path), 2);
for izone = 1:2 % 1 codes for preTest, 2 for postTest
    occup_shock_mean(:,izone) = mean(squeeze(occup_shock(:,:,izone)),2);
end

% Prepare the 'first enter to shock zone' array
EntryTime_shock = nan(length(Dir.path), 4, 2); % 4 tests, Pre and Post
for i = 1:length(a)
    for k=1:length(id_Pre{i})
        temp = CalculateFirstEntryZoneTime(a{i}.behavResources(id_Pre{i}(k)), 240);
        EntryTime_shock(i,k,1) = temp(1);
        temp = CalculateFirstEntryZoneTime(a{i}.behavResources(id_Post{i}(k)), 240);
        EntryTime_shock(i,k,2) = temp(1);
    end  
end
    
EntryTime_shock_mean = nan(length(Dir.path), 2);
for izone = 1:2 % 1 codes for preTest, 2 for postTest
    EntryTime_shock_mean(:,izone) = mean(squeeze(EntryTime_shock(:,:,izone)),2);
end

Diff_occup = diff(occup_shock_mean,1,2)*100;
Diff_entrytime = diff(EntryTime_shock_mean,1,2);

%% Choose epochs
num_epoch_wake = strcmp(states_wake, EV_epoch);
num_epoch_sleep = strcmp(states_sleep, type_sleep);
EV_4corr = nonzeros(EV{num_epoch_sleep}{num_epoch_wake});
EV_4corr(isnan(EV_4corr)) = [];

%% Calculate correlations
% Occupancy
[rho, p] = corr(EV_4corr, Diff_occup, 'Type', 'Spearman'); 
r(1) = rho;
pvalue(1) = p;

% Entries
[rho, p] = corr(EV_4corr, Diff_entrytime, 'Type', 'Spearman'); 
r(2) = rho;
pvalue(2) = p;

%% Plot figures
f= figure('units', 'normalized', 'outerposition', [0 0 1 0.65]);
subplot(121)
scatter(EV_4corr*100, Diff_occup, 150, 'k', 'filled')
t1(izone) = text(.4, .9, ['rho = ' num2str(round(r(1), 3))], 'sc', 'FontSize', 13);
t2(izone) = text(.4, .8, ['p = ' num2str(round(pvalue(1), 3))], 'sc', 'FontSize', 13);
l1 = lsline;
l1.Color = 'k';
l1.LineWidth = 1.5;
xlabel('EV');
ylabel('Post - Pre SZ occupancies');
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Helvetica');
set(gca, 'LineWidth', 1);
makepretty_DB

subplot(122)
scatter(EV_4corr*100, Diff_entrytime, 150, 'k', 'filled')
t1(izone) = text(.4, .9, ['rho = ' num2str(round(r(2), 3))], 'sc', 'FontSize', 13);
t2(izone) = text(.4, .8, ['p = ' num2str(round(pvalue(2), 3))], 'sc', 'FontSize', 13);
l1 = lsline;
l1.Color = 'k';
l1.LineWidth = 1.5;
xlabel('EV');
ylabel('Post - Pre entry latencies');
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Helvetica');
set(gca, 'LineWidth', 1);
makepretty_DB

if sav
    foldertosave = ChooseFolderForFigures_DB('ReactReplay');
    saveas(f,[foldertosave filesep 'EV' filesep 'Correlations' filesep 'CorrEV_' type_sleep '_' EV_epoch '.fig']);
    saveFigure(f,['CorrEV_' type_sleep '_' EV_epoch], [foldertosave filesep 'EV' filesep 'Correlations']);
end

end