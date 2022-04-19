function fh = BehaviorERC(nMice, experiment, varargin)
%
% This function plots 4 behavioral metrics of the Shock Zone in the UMaze
% for Pre- and Post-Tests:
%  - Occupancy of the zone
%  - Number of entries into the zone
%  - Time to enter the zone for the first time
%  - Average speed in the zone
%
% INPUT
%
%     nMice      numbers of mice to include in the analysis
%     experiment type of ERC experiment ('PAG', 'MFB' or 'Novel')
%     IsSave     (optional) true or 1 if to save the figure (default=false)
%     SownPoints (optional) 1 if you want to show individual points (default=0)
% 
%  OUTPUT
%
%     nothing
%
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 05/11/2020 based on the script I wrote in 2018
% github.com/bryzgalovdm

%% Parameters
FigName = 'UMazePAG_Behavior';
IsSave = false;
SP = 1;
numtest = 4;
% nMice = [797 798 828 861 882 905 906 911 912 977 994 1117 1124 1161 1162 1168];

%% Optional Parameters
for i=1:2:length(varargin)    
    switch(lower(varargin{i})) 
        case 'issave'
            IsSave = varargin{i+1};
            if IsSave ~= 1 && IsSave ~= 0
                error('Incorrect value for property ''IsSave'' (type ''help PreTestCharacteristics'' for details).');
            end
        case 'showpoints'
            SP = varargin{i+1};
            if SP ~= 1 && SP ~= 0
                error('Incorrect value for property ''ShowPoints'' (type ''help PreTestCharacteristics'' for details).');
            end
        case 'numtest'
            numtest =  varargin{i+1};
            if ~isa(numtest, 'double') && numtest < 1
                error('Incorrect value for property ''NumTest'' (type ''help ExampleTrajectory_DB'' for details).');
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

%% Get data
Dir = PathForExperimentsERC(fetchpaths);
Dir = RestrictPathForExperiment(Dir,'nMice', nMice);
numsessions = CountNumSesionsERC(Dir);

a = cell(numsessions, 1);
cnt = 1;
for imouse = 1:length(Dir.path)
    for isession = 1:length(Dir.path{imouse})
        a{cnt} = load([Dir.path{imouse}{isession} '/behavResources.mat'], 'behavResources');
        cnt = cnt + 1;
    end
end

%% Find indices of PreTests and PostTest session in the structure
id_Pre = cell(1,length(a));
id_Post = cell(1,length(a));

for i=1:length(a)
    id_Pre{i} = FindSessionID_ERC(a{i}.behavResources, 'TestPre');
    id_Post{i} = FindSessionID_ERC(a{i}.behavResources, 'TestPost');
end

%% Calculate average occupancy
% Calculate occupancy de novo
occup_shock = nan(numsessions, numtest, 2); % 4 tests, Pre and Post
occup_safe = nan(numsessions, numtest, 2); % 4 tests, Pre and Post

for i=1:length(a)
    for k=1:length(id_Pre{i})
        temp = CalculateZoneOccupancy(a{i}.behavResources(id_Pre{i}(k)));
        occup_shock(i,k,1) = temp(1);
        occup_safe(i,k,1) = temp(2);
        temp = CalculateZoneOccupancy(a{i}.behavResources(id_Post{i}(k)));
        occup_shock(i,k,2) = temp(1);
        occup_safe(i,k,2) = temp(2);
    end
end

occup_shock_mean = nan(numsessions, 2);
occup_shock_std = nan(numsessions, 2);
occup_safe_mean = nan(numsessions, 2);
occup_safe_std = nan(numsessions, 2);

for izone = 1:2 % 1 codes for preTest, 2 for postTest
    occup_shock_mean(:,izone) = mean(squeeze(occup_shock(:,:,izone)),2);
    occup_shock_std(:,izone) = std(squeeze(occup_shock(:,:,izone)),0,2);
    occup_safe_mean(:,izone) = mean(squeeze(occup_safe(:,:,izone)),2);
    occup_safe_std(:,izone) = std(squeeze(occup_safe(:,:,izone)),0,2);
end

%% Prepare the 'first enter to shock zone' array

EntryTime_shock = nan(numsessions, numtest, 2); % 4 tests, Pre and Post
EntryTime_safe = nan(numsessions, numtest, 2); % 4 tests, Pre and Post

for i = 1:length(a)
    for k=1:length(id_Pre{i})
        temp = CalculateFirstEntryZoneTime(a{i}.behavResources(id_Pre{i}(k)), 240);
        EntryTime_shock(i,k,1) = temp(1);
        EntryTime_safe(i,k,1) = temp(2);
        temp = CalculateFirstEntryZoneTime(a{i}.behavResources(id_Post{i}(k)), 240);
        EntryTime_shock(i,k,2) = temp(1);
        EntryTime_safe(i,k,2) = temp(2);
    end  
end
    
EntryTime_shock_mean = nan(numsessions, 2);
EntryTime_shock_std = nan(numsessions, 2);
EntryTime_safe_mean = nan(numsessions, 2);
EntryTime_safe_std = nan(numsessions, 2);
for izone = 1:2 % 1 codes for preTest, 2 for postTest
    EntryTime_shock_mean(:,izone) = mean(squeeze(EntryTime_shock(:,:,izone)),2);
    EntryTime_shock_std(:,izone) = std(squeeze(EntryTime_shock(:,:,izone)),0,2);
    EntryTime_safe_mean(:,izone) = mean(squeeze(EntryTime_safe(:,:,izone)),2);
    EntryTime_safe_std(:,izone) = std(squeeze(EntryTime_safe(:,:,izone)),0,2);
end

%% Calculate number of entries into the shock zone

NumEntries_shock = nan(numsessions, numtest, 2); % 4 tests, Pre and Post
NumEntries_safe = nan(numsessions, numtest, 2); % 4 tests, Pre and Post

for i = 1:length(a)
    for k=1:length(id_Pre{i})
        temp = CalculateNumEntriesZone(a{i}.behavResources(id_Pre{i}(k)));
        NumEntries_shock(i,k,1) = temp(1);
        NumEntries_safe(i,k,1) = temp(2);
        temp = CalculateNumEntriesZone(a{i}.behavResources(id_Post{i}(k)));
        NumEntries_shock(i,k,2) = temp(1);
        NumEntries_safe(i,k,2) = temp(2);
    end  
end

NumEntries_shock_mean = nan(numsessions, 2);
NumEntries_shock_std = nan(numsessions, 2);
NumEntries_safe_mean = nan(numsessions, 2);
NumEntries_safe_std = nan(numsessions, 2);
for izone = 1:2 % 1 codes for preTest, 2 for postTest
    NumEntries_shock_mean(:,izone) = mean(squeeze(NumEntries_shock(:,:,izone)),2);
    NumEntries_shock_std(:,izone) = std(squeeze(NumEntries_shock(:,:,izone)),0,2);
    NumEntries_safe_mean(:,izone) = mean(squeeze(NumEntries_safe(:,:,izone)),2);
    NumEntries_safe_std(:,izone) = std(squeeze(NumEntries_safe(:,:,izone)),0,2);
end

%% Calculate speed in the safe zone and in the noshock + shock vs everything else
% I skip the last point in ZoneIndices because length(Xtsd)=length(Vtsd)+1
% - UPD 18/07/2018 - Could do length(Start(ZoneEpoch))
Speed_shock = nan(numsessions, numtest, 2); % 4 tests, Pre and Post
Speed_safe = nan(numsessions, numtest, 2); % 4 tests, Pre and Post

for i = 1:length(a)
    for k=1:length(id_Pre{i})
        temp = CalculateSpeedZone(a{i}.behavResources(id_Pre{i}(k)));
        Speed_shock(i,k,1) = temp(1);
        Speed_safe(i,k,1) = temp(2);
        temp = CalculateSpeedZone(a{i}.behavResources(id_Post{i}(k)));
        Speed_shock(i,k,2) = temp(1);
        Speed_safe(i,k,2) = temp(2);
    end 
end

Speed_shock_mean = nan(numsessions, 2);
Speed_shock_std = nan(numsessions, 2);
Speed_safe_mean = nan(numsessions, 2);
Speed_safe_std = nan(numsessions, 2);
for izone = 1:2 % 1 codes for preTest, 2 for postTest
    Speed_shock_mean(:,izone) = mean(squeeze(Speed_shock(:,:,izone)),2);
    Speed_shock_std(:,izone) = std(squeeze(Speed_shock(:,:,izone)),0,2);
    Speed_safe_mean(:,izone) = mean(squeeze(Speed_safe(:,:,izone)),2);
    Speed_safe_std(:,izone) = std(squeeze(Speed_safe(:,:,izone)),0,2);
end

%% Plot
if strcmp(experiment, 'PAG')
    cols = {[.9 0 0], [0 0 0]};
elseif strcmp(experiment, 'MFB')
    cols = {[0 .9 0], [0 0 0]};
elseif strcmp(experiment, 'Novel')
    cols = {[0.793 .902 0.184], [0 0 0]};
elseif strcmp(experiment, 'Known')
    cols = {[0.793 .902 0.184], [0 0 0]};
end
% Axes
fh = figure('units', 'normalized', 'outerposition', [0.1 0.2 0.53 0.85]);
Occupancy_Axes = axes('position', [0.09 0.55 0.41 0.41]);
First_Axes = axes('position', [0.57 0.55 0.41 0.41]);
NumEntr_Axes = axes('position', [0.09 0.05 0.41 0.41]);
Speed_Axes = axes('position', [0.57 0.05 0.41 0.41]);

% Occupancy
axes(Occupancy_Axes);
[~,h_occ] = PlotErrorBarN_DB([occup_shock_mean(:,1)*100 occup_shock_mean(:,2)*100],...
    'barcolors', cols{2}, 'barwidth', 0.4, 'newfig', 0, 'showpoints', SP);
h_occ.FaceColor = 'flat';
h_occ.FaceAlpha = .6;
h_occ.CData(2,:) = cols{1};
set(gca,'Xtick',[1:2],'XtickLabel',{'PreTests', 'PostTests'});
makepretty
line(xlim,[21.5 21.5],'Color','k','LineStyle','--','LineWidth',5);
ylim([0 70])
% text(1.85,23.2,'Random Occupancy', 'FontWeight','bold','FontSize',13);
ylabel('% time');
title('Shock zone occupancy')
makepretty

axes(NumEntr_Axes);
[~,h_nent] = PlotErrorBarN_DB([NumEntries_shock(:,1) NumEntries_shock(:,2)],...
    'barcolors', cols{2}, 'barwidth', 0.4, 'newfig', 0, 'showpoints',SP);
h_nent.FaceColor = 'flat';
h_nent.FaceAlpha = .6;
h_nent.CData(2,:) = cols{1};
set(gca,'Xtick',[1:2],'XtickLabel',{'PreTests', 'PostTests'});
ylim([0 7])
ylabel('Number of entries');
title('Entries to shock zone');
makepretty

axes(First_Axes);
[~,h_first] = PlotErrorBarN_DB([EntryTime_shock_mean(:,1) EntryTime_shock_mean(:,2)],...
    'barcolors', cols{2}, 'barwidth', 0.4, 'newfig', 0, 'showpoints',SP);
h_first.FaceColor = 'flat';
h_first.FaceAlpha = .6;
h_first.CData(2,:) = cols{1};
set(gca,'Xtick',[1:2],'XtickLabel',{'PreTests', 'PostTests'});
ylim([0 250])
ylabel('Time (s)');
title('Latency to enter shock zone')
makepretty

axes(Speed_Axes);
[~,h_speed] = PlotErrorBarN_DB([Speed_shock_mean(:,1) Speed_shock_mean(:,2)],...
    'barcolors', cols{2}, 'barwidth', 0.4, 'newfig', 0, 'showpoints',SP);
h_speed.FaceColor = 'flat';
h_speed.FaceAlpha = .6;
h_speed.CData(2,:) = cols{1};
set(gca,'Xtick',[1:2],'XtickLabel',{'PreTests', 'PostTests'});
ylim([0 6])
ylabel('Speed (cm/s)');
title('Speed in safe zone')
makepretty

%% Save
if IsSave
    foldertosave = ChooseFolderForFigures_DB('Behavior');
    saveas(fh,[foldertosave '/' FigName '_' experiment '.fig']);
    saveFigure(fh, FigName, foldertosave);
end

%% %% Plot difference shock-safe
% Axes
f1 = figure('units', 'normalized', 'outerposition', [0.1 0.2 0.53 0.85]);
Occupancy_Axes = axes('position', [0.09 0.55 0.41 0.41]);
First_Axes = axes('position', [0.57 0.55 0.41 0.41]);
NumEntr_Axes = axes('position', [0.09 0.05 0.41 0.41]);
Speed_Axes = axes('position', [0.57 0.05 0.41 0.41]);

% Occupancy
axes(Occupancy_Axes);
[~,h_occ] = PlotErrorBarN_DB([occup_shock_mean(:,1)*100-occup_safe_mean(:,1)*100 ...
    occup_shock_mean(:,2)*100-occup_safe_mean(:,2)*100],...
    'barcolors', cols{2}, 'barwidth', 0.4, 'newfig', 0, 'showpoints', 0);
h_occ.FaceColor = 'flat';
h_occ.FaceAlpha = .6;
h_occ.CData(2,:) = cols{1};
set(gca,'Xtick',[1:2],'XtickLabel',{'PreTests', 'PostTests'});
makepretty
ylim([-30 50])
ylabel('Shock-Safe: % time');
title('Zone occupancy')
makepretty

axes(NumEntr_Axes);
[~,h_nent] = PlotErrorBarN_DB([NumEntries_shock(:,1)-NumEntries_safe(:,1) ...
    NumEntries_shock(:,2)-NumEntries_safe(:,2)],...
    'barcolors', cols{2}, 'barwidth', 0.4, 'newfig', 0, 'showpoints',0);
h_nent.FaceColor = 'flat';
h_nent.FaceAlpha = .6;
h_nent.CData(2,:) = cols{1};
set(gca,'Xtick',[1:2],'XtickLabel',{'PreTests', 'PostTests'});
ylim([-1.2 .8])
ylabel('Shock-Safe: # of entries');
title('Entries to the zone');
makepretty

axes(First_Axes);
[~,h_first] = PlotErrorBarN_DB([EntryTime_shock_mean(:,1)-EntryTime_safe_mean(:,1) ...
    EntryTime_shock_mean(:,2)-EntryTime_safe_mean(:,2)],...
    'barcolors', cols{2}, 'barwidth', 0.4, 'newfig', 0, 'showpoints',0);
h_first.FaceColor = 'flat';
h_first.FaceAlpha = .6;
h_first.CData(2,:) = cols{1};
set(gca,'Xtick',[1:2],'XtickLabel',{'PreTests', 'PostTests'});
ylim([-120 120])
ylabel('Shock-Safe: time (s)');
title('Latency to enter the zone zone')
makepretty

axes(Speed_Axes);
[~,h_speed] = PlotErrorBarN_DB([Speed_shock_mean(:,1)-Speed_safe_mean(:,1) ...
    Speed_shock_mean(:,2)-Speed_safe_mean(:,2)],...
    'barcolors', cols{2}, 'barwidth', 0.4, 'newfig', 0, 'showpoints',0);
h_speed.FaceColor = 'flat';
h_speed.FaceAlpha = .6;
h_speed.CData(2,:) = cols{1};
set(gca,'Xtick',[1:2],'XtickLabel',{'PreTests', 'PostTests'});
ylim([-1.2 1.5])
ylabel('Shock-Safe: speed (cm/s)');
title('Speed in the zone')
makepretty

%% Save
if IsSave
    foldertosave = ChooseFolderForFigures_DB('Behavior');
    saveas(f1,[foldertosave '/' FigName '_' experiment '_diff.fig']);
    saveFigure(f1, [FigName '_diff'], foldertosave);
end


%% Plot neutral
% Axes
f2 = figure('units', 'normalized', 'outerposition', [0.1 0.2 0.53 0.85]);
Occupancy_Axes = axes('position', [0.09 0.55 0.41 0.41]);
First_Axes = axes('position', [0.57 0.55 0.41 0.41]);
NumEntr_Axes = axes('position', [0.09 0.05 0.41 0.41]);
Speed_Axes = axes('position', [0.57 0.05 0.41 0.41]);

% Occupancy
axes(Occupancy_Axes);
[~,h_occ] = PlotErrorBarN_DB([occup_safe_mean(:,1)*100 occup_safe_mean(:,2)*100],...
    'barcolors', cols{2}, 'barwidth', 0.4, 'newfig', 0, 'showpoints', SP);
h_occ.FaceColor = 'flat';
h_occ.FaceAlpha = .6;
h_occ.CData(2,:) = cols{1};
set(gca,'Xtick',[1:2],'XtickLabel',{'PreTests', 'PostTests'});
makepretty
line(xlim,[21.5 21.5],'Color','k','LineStyle','--','LineWidth',5);
ylim([0 70])
% text(1.85,23.2,'Random Occupancy', 'FontWeight','bold','FontSize',13);
ylabel('% time');
title('Safe zone occupancy')
makepretty

axes(NumEntr_Axes);
[~,h_nent] = PlotErrorBarN_DB([NumEntries_safe(:,1) NumEntries_safe(:,2)],...
    'barcolors', cols{2}, 'barwidth', 0.4, 'newfig', 0, 'showpoints',SP);
h_nent.FaceColor = 'flat';
h_nent.FaceAlpha = .6;
h_nent.CData(2,:) = cols{1};
set(gca,'Xtick',[1:2],'XtickLabel',{'PreTests', 'PostTests'});
ylim([0 7])
ylabel('Number of entries');
title('Entries to safe zone');
makepretty

axes(First_Axes);
[~,h_first] = PlotErrorBarN_DB([EntryTime_safe_mean(:,1) EntryTime_safe_mean(:,2)],...
    'barcolors', cols{2}, 'barwidth', 0.4, 'newfig', 0, 'showpoints',SP);
h_first.FaceColor = 'flat';
h_first.FaceAlpha = .6;
h_first.CData(2,:) = cols{1};
set(gca,'Xtick',[1:2],'XtickLabel',{'PreTests', 'PostTests'});
ylim([0 250])
ylabel('Time (s)');
title('Latency to enter safe zone')
makepretty

axes(Speed_Axes);
[~,h_speed] = PlotErrorBarN_DB([Speed_safe_mean(:,1) Speed_safe_mean(:,2)],...
    'barcolors', cols{2}, 'barwidth', 0.4, 'newfig', 0, 'showpoints',SP);
h_speed.FaceColor = 'flat';
h_speed.FaceAlpha = .6;
h_speed.CData(2,:) = cols{1};
set(gca,'Xtick',[1:2],'XtickLabel',{'PreTests', 'PostTests'});
ylim([0 6])
ylabel('Speed (cm/s)');
title('Speed in safe zone')
makepretty

if IsSave
    foldertosave = ChooseFolderForFigures_DB('Behavior');
    saveas(f2,[foldertosave '/' FigName '_' experiment '_neutral.fig']);
    saveFigure(f2, [FigName '_neutral'], foldertosave);
end

% %% Save_data
% if save_res == 1
%     diff_occup = [occup_shock_mean(:,1)*100-occup_safe_mean(:,1)*100 ...
%     occup_shock_mean(:,2)*100-occup_safe_mean(:,2)*100];
%     save([foldertosave filesep 'Behav_ratios.mat'], 'ratio_occup', 'diff_occup', 'ratio_fentry', 'diff_fentry', 'ratio_nentries', 'diff_nentries', 'ratio_speed', 'diff_speed');
% end