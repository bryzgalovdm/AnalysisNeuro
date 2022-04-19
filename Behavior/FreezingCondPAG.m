function FreezingCondPAG(mice)
% 
% This function plots two figures of freezing in shock and safe
% compartments: one for learning sessions and one for post-tests
% 
% INPUT
%   
%   mice                    mice to analyse
% 
% OUTPUT
% 
% 
% EXAMPLE
% 
%   FreezingCondPAG([797 798])
% 
% SEE ALSO
%
% 
% By Dima Bryzgalov, MOBS team, Paris,
% 08-09/2021
% github.com/bryzgalovdm

%% Get data
Dir = PathForExperimentsERC('UMazePAG');
Dir = RestrictPathForExperiment(Dir,'nMice',mice);
numsessions = CountNumSesionsERC(Dir);

%% Allocate
CondEpoch = cell(numsessions, 1);
TestPostEpoch = cell(numsessions, 1);
FreezeEpoch = cell(numsessions, 1);
FreezeEpochShock = cell(numsessions, 1);
FreezeEpochSafe = cell(numsessions, 1);
FreezeEpochPost = cell(numsessions, 1);
FreezeEpochPostShock = cell(numsessions, 1);
FreezeEpochPostSafe = cell(numsessions, 1);
OverallFreeze = nan(numsessions, 1);
ShockFreeze = nan(numsessions, 1);
SafeFreeze = nan(numsessions, 1);
OverallFreezePost = nan(numsessions, 1);
ShockFreezePost = nan(numsessions, 1);
SafeFreezePost = nan(numsessions, 1);

%% Load dara
cnt=1;
for imouse=1:length(Dir.path)
    for isession = 1:length(Dir.path{imouse})
        b{cnt} = load([Dir.path{imouse}{isession} 'behavResources.mat'],'SessionEpoch', 'FreezeAccEpoch', 'ZoneEpoch');
        try
            sleep{cnt} = load([Dir.path{imouse}{isession} 'SleepScoring_OBGamma.mat'],'SWSEpoch','REMEpoch','Wake', 'TotalNoiseEpoch');
        catch
            sleep{cnt} = load([Dir.path{imouse}{isession} 'SleepScoring_Accelero.mat'],'SWSEpoch','REMEpoch','Wake', 'TotalNoiseEpoch');
        end
        cnt=cnt+1;
    end
end



%% Prepare epochs
for isession = 1:numsessions
    [~, ~, ~, CondEpoch{isession}, ~, TestPostEpoch{isession}] = ReturnMnemozyneEpochs(b{isession}.SessionEpoch);
    FreezeEpoch{isession} = and(sleep{isession}.Wake, and(b{isession}.FreezeAccEpoch, CondEpoch{isession}));
    FreezeEpochShock{isession} = and(FreezeEpoch{isession}, b{isession}.ZoneEpoch.Shock);
    FreezeEpochSafe{isession} = and(FreezeEpoch{isession}, b{isession}.ZoneEpoch.NoShock);
    
    FreezeEpochPost{isession} = and(sleep{isession}.Wake, and(b{isession}.FreezeAccEpoch, TestPostEpoch{isession}));
    FreezeEpochPostShock{isession} = and(FreezeEpochPost{isession}, b{isession}.ZoneEpoch.Shock);
    FreezeEpochPostSafe{isession} = and(FreezeEpochPost{isession}, b{isession}.ZoneEpoch.NoShock);
end

%% Calculate percentage of freezing
for isession = 1:numsessions
    OverallFreeze(isession) = tot_length(FreezeEpoch{isession})/tot_length(CondEpoch{isession});
    ShockFreeze(isession) = tot_length(FreezeEpochShock{isession})/tot_length(CondEpoch{isession});
    SafeFreeze(isession) = tot_length(FreezeEpochSafe{isession})/tot_length(CondEpoch{isession});
    
    OverallFreezePost(isession) = tot_length(FreezeEpochPost{isession})/tot_length(TestPostEpoch{isession});
    ShockFreezePost(isession) = tot_length(FreezeEpochPostShock{isession})/tot_length(TestPostEpoch{isession});
    SafeFreezePost(isession) = tot_length(FreezeEpochPostSafe{isession})/tot_length(TestPostEpoch{isession});
end

%% Plot
fh = figure('units', 'normalized', 'outerposition', [0.1 0.2 0.4 0.7]);

% Occupancy
[~,h_occ] = PlotErrorBarN_DB([OverallFreeze*100 ShockFreeze*100 SafeFreeze*100],...
    'barcolors', [0 0 0], 'barwidth', 0.4, 'newfig', 0, 'showpoints', 1);
h_occ.FaceColor = 'flat';
h_occ.FaceAlpha = .6;
h_occ.CData(2,:) = [.9 0 0];
h_occ.CData(3,:) = [0 0 0.9];
set(gca,'Xtick',[1:3],'XtickLabel',{'Full Maze', 'Shock', 'Safe'});
ylabel('% time from the whole conditioning');
title('Freezing during condtioning')
makepretty


fh1 = figure('units', 'normalized', 'outerposition', [0.1 0.2 0.4 0.7]);

% Occupancy
[~,h_occ] = PlotErrorBarN_DB([OverallFreezePost*100 ShockFreezePost*100 SafeFreezePost*100],...
    'barcolors', [0 0 0], 'barwidth', 0.4, 'newfig', 0, 'showpoints', 1);
h_occ.FaceColor = 'flat';
h_occ.FaceAlpha = .6;
h_occ.CData(2,:) = [.9 0 0];
h_occ.CData(3,:) = [0 0 0.9];
set(gca,'Xtick',[1:3],'XtickLabel',{'Full Maze', 'Shock', 'Safe'});
ylabel('% time from the whole conditioning');
title('Freezing during PostTests (16 min)')
makepretty

