function fh = BehaviorERC_PAGvsMFB(pagMice, mfbMice, varargin)
%
% This function plots 4 behavioral metrics of the stim Zone in the UMaze
% for Pre- and Post-Tests:
%  - Occupancy of the zone
%  - Number of entries into the zone
%  - Time to enter the zone for the first time
%  - Average speed in the zone zone
%
% INPUT
%
%     nMice      numbers of mice to include in the analysis
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
FigName = 'PAGvsMFB_Behavior';
IsSave = false;
SP = 0;

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
    end
end

%% Get data
% PAG 
Dir{1} = PathForExperimentsERC_Dima('UMazePAG');
Dir{1} = RestrictPathForExperiment(Dir,'nMice', pagMice);
% MFB
Dir{2} = PathForExperimentsERC_SL(expe);
Dir{2} = RestrictPathForExperiment(Dir,'nMice', Mice_to_analyze);

nDir = length(length(Dir{1}.path));
nDirMfb = length(length(Dir{2}.path));
for igr=1:2
    nDir(igr) = length(length(Dir{igr}.path));
    a{igrg} = cell(nDirPag,1);
    for i = 1:length(Dir.path)
        % PreTests
        a{i} = load([Dir.path{i}{1} '/behavResources.mat'], 'behavResources');
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
occup = nan(length(Dir.path), 4, 2); % 4 tests, Pre and Post

for i=1:length(a)
    for k=1:length(id_Pre{i})
        temp = CalculateZoneOccupancy(a{i}.behavResources(id_Pre{i}(k)));
        occup(i,k,1) = temp(1);
        temp = CalculateZoneOccupancy(a{i}.behavResources(id_Post{i}(k)));
        occup(i,k,2) = temp(1);
    end
end

occup_mean = nan(length(Dir.path), 2);
occup_std = nan(length(Dir.path), 2);

for izone = 1:2 % 1 codes for preTest, 2 for postTest
    occup_mean(:,izone) = mean(squeeze(occup(:,:,izone)),2);
    occup_std(:,izone) = std(squeeze(occup(:,:,izone)),0,2);
end

%% Prepare the 'first enter to shock zone' array

EntryTime = nan(length(Dir.path), 4, 2); % 4 tests, Pre and Post

for i = 1:length(a)
    for k=1:length(id_Pre{i})
        temp = CalculateFirstEntryZoneTime(a{i}.behavResources(id_Pre{i}(k)), 240);
        EntryTime(i,k,1) = temp(1);
        temp = CalculateFirstEntryZoneTime(a{i}.behavResources(id_Post{i}(k)), 240);
        EntryTime(i,k,2) = temp(1);
    end  
end
    
EntryTime_mean = nan(length(Dir.path), 2);
EntryTime_std = nan(length(Dir.path), 2);
for izone = 1:2 % 1 codes for preTest, 2 for postTest
    EntryTime_mean(:,izone) = mean(squeeze(EntryTime(:,:,izone)),2);
    EntryTime_std(:,izone) = std(squeeze(EntryTime(:,:,izone)),0,2);
end

%% Calculate number of entries into the shock zone

NumEntries = nan(length(Dir.path), 4, 2); % 4 tests, Pre and Post

for i = 1:length(a)
    for k=1:length(id_Pre{i})
        temp = CalculateNumEntriesZone(a{i}.behavResources(id_Pre{i}(k)));
        NumEntries(i,k,1) = temp(1);
        temp = CalculateNumEntriesZone(a{i}.behavResources(id_Post{i}(k)));
        NumEntries(i,k,2) = temp(1);
    end  
end

NumEntries_mean = nan(length(Dir.path), 2);
NumEntries_std = nan(length(Dir.path), 2);
for izone = 1:2 % 1 codes for preTest, 2 for postTest
    NumEntries_mean(:,izone) = mean(squeeze(NumEntries(:,:,izone)),2);
    NumEntries_std(:,izone) = std(squeeze(NumEntries(:,:,izone)),0,2);
end

%% Calculate speed in the safe zone and in the noshock + shock vs everything else
% I skip the last point in ZoneIndices because length(Xtsd)=length(Vtsd)+1
% - UPD 18/07/2018 - Could do length(Start(ZoneEpoch))
Speed = nan(length(Dir.path), 4, 2); % 4 tests, 2 zones

for i = 1:length(a)
    for k=1:length(id_Pre{i})
        temp = CalculateSpeedZone(a{i}.behavResources(id_Pre{i}(k)));
        Speed(i,k,1) = temp(1);
        temp = CalculateSpeedZone(a{i}.behavResources(id_Post{i}(k)));
        Speed(i,k,2) = temp(1);
    end 
end

Speed_mean = nan(length(Dir.path), 2);
Speed_std = nan(length(Dir.path), 2);
for izone = 1:2 % 1 codes for preTest, 2 for postTest
    Speed_mean(:,izone) = mean(squeeze(Speed(:,:,izone)),2);
    Speed_std(:,izone) = std(squeeze(Speed(:,:,izone)),0,2);
end

%% Plot
% Axes
fh = figure('units', 'normalized', 'outerposition', [0 0 0.65 0.65]);
Occupancy_Axes = axes('position', [0.07 0.55 0.41 0.41]);
NumEntr_Axes = axes('position', [0.55 0.55 0.41 0.41]);
First_Axes = axes('position', [0.07 0.05 0.41 0.41]);
Speed_Axes = axes('position', [0.55 0.05 0.41 0.41]);

% Occupancy
axes(Occupancy_Axes);
[p_occ,h_occ, her_occ] = PlotErrorBarN_DB([occup_mean(:,1)*100 occup_mean(:,2)*100],...
    'barcolors', [1 0 0], 'barwidth', 0.6, 'newfig', 0, 'showpoints', SP);
h_occ.FaceColor = 'flat';
h_occ.CData(2,:) = [0 0 1];
set(gca,'Xtick',[1:2],'XtickLabel',{'Shock', 'Safe'});
set(gca, 'FontSize', 18, 'FontWeight',  'bold','FontName','Times New Roman');
set(gca, 'LineWidth', 3);
set(h_occ, 'LineWidth', 3);
set(her_occ, 'LineWidth', 3);
line(xlim,[21.5 21.5],'Color','k','LineStyle','--','LineWidth',5);
% text(1.85,23.2,'Random Occupancy', 'FontWeight','bold','FontSize',13);
ylabel('% time');
title('Occupancy percentage', 'FontSize', 14);
% ylim([0 50])

axes(NumEntr_Axes);
[p_nent,h_nent, her_nent] = PlotErrorBarN_DB([NumEntries_mean(:,1) NumEntries_mean(:,2)],...
    'barcolors', [1 0 0], 'barwidth', 0.6, 'newfig', 0, 'showpoints',SP);
h_nent.FaceColor = 'flat';
h_nent.CData(2,:) = [0 0 1];
set(gca,'Xtick',[1:2],'XtickLabel',{'Shock', 'Safe'});
set(gca, 'FontSize', 18, 'FontWeight',  'bold','FontName','Times New Roman');
set(gca, 'LineWidth', 3);
set(h_nent, 'LineWidth', 3);
set(her_nent, 'LineWidth', 3);
ylabel('Number of entries');
title('# of entries to the Zone', 'FontSize', 14);
% ylim([0 6])

axes(First_Axes);
[p_first,h_first, her_first] = PlotErrorBarN_DB([EntryTime_mean(:,1) EntryTime_mean(:,2)],...
    'barcolors', [1 0 0], 'barwidth', 0.6, 'newfig', 0, 'showpoints',SP);
h_first.FaceColor = 'flat';
h_first.CData(2,:) = [0 0 1];
set(gca,'Xtick',[1:2],'XtickLabel',{'Shock', 'Safe'});
set(gca, 'FontSize', 18, 'FontWeight', 'bold','FontName','Times New Roman');
set(gca, 'LineWidth', 3);
set(h_first, 'LineWidth', 3);
set(her_first, 'LineWidth', 3);
ylabel('Time (s)');
title('First time to enter the zone', 'FontSize', 14);

axes(Speed_Axes);
[p_speed,h_speed, her_speed] = PlotErrorBarN_DB([Speed_mean(:,1) Speed_mean(:,2)],...
    'barcolors', [1 0 0], 'barwidth', 0.6, 'newfig', 0, 'showpoints',SP);
h_speed.FaceColor = 'flat';
h_speed.CData(2,:) = [0 0 1];
set(gca,'Xtick',[1:2],'XtickLabel',{'Shock', 'Safe'});
set(gca, 'FontSize', 18, 'FontWeight',  'bold','FontName','Times New Roman');
set(gca, 'LineWidth', 3);
set(h_speed, 'LineWidth', 3);
set(her_speed, 'LineWidth', 3);
ylabel('Speed (cm/s)');
title('Average speed in the zone', 'FontSize', 14);
% ylim([0 8])

%% Save
if IsSave
    saveas(fh,[pwd '/' FigName '.fig']);
    saveFigure(fh,pwd, dirsave);
end