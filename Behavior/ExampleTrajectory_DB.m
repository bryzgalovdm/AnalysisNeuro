
function ExampleTrajectory_DB(Dir, experiment, IsPCDriven, varargin)
%
% This functions plot trajectories of a selected mouse in pretests,
% condiioning sessions and posttests as well as bar plots of shock-safe
%  zone occupancies
%
%  INPUT
%       
%       Dir             directory with behavResources.mat
%       ISPCDriven      true if the mouse was stimulated using BCI and info
%                           about stims is stored only in StimEpoch (boolean)
%       IsSave          true if you want to save the resulting figure in
%                           Dir (boolean - default = false)
%       NumTest         number of tests to plot the trajectories (default = 4)
% 
% 
%  OUTPUT
%
%       Figure
%
%       See
%   
%       BehaviorERC, PathForExperimentERC_Dima
% 
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 27/10/2020, based on my own script coded in 2019
% github.com/bryzgalovdm


%% Default parameters
IsSave = false;
numtest=4;
IsColorTests = false;

%% Optional parameters
for i=1:2:length(varargin)    
    switch(lower(varargin{i})) 
        case 'issave'
            IsSave = varargin{i+1};
            if IsSave ~= 1 && IsSave ~= 0
                error('Incorrect value for property ''IsSave'' (type ''help ExampleTrajectory_DB'' for details).');
            end
        case 'numtest'
            numtest =  varargin{i+1};
            if ~isa(numtest, 'double')
                error('Incorrect value for property ''NumTest'' (type ''help ExampleTrajectory_DB'' for details).');
            end
        case 'colortests'
            IsColorTests =  varargin{i+1};
            if IsColorTests ~= 1 && IsColorTests ~= 0
                error('Incorrect value for property ''ColorTests'' (type ''help ExampleTrajectory_DB'' for details).');
            end    
    end
end

%% Get data
a = load([Dir '/behavResources.mat'], 'behavResources', 'TTLInfo', 'SessionEpoch');

%% Find indices of PreTests and PostTest session in the structure

id_Pre = zeros(1,numtest);
id_Cond = cell(1,numtest);
if isfield(a.SessionEpoch, 'TestPost1')
    id_Post = cell(1,numtest);
end

id_Pre = zeros(1,length(a.behavResources));
id_Cond = zeros(1,length(a.behavResources));
if isfield(a.SessionEpoch, 'TestPost1')
    id_Post = zeros(1,length(a.behavResources));
end

for k=1:length(a.behavResources)
    if ~isempty(strfind(a.behavResources(k).SessionName,'TestPre'))
        id_Pre(k) = 1;
    end
    if ~isempty(strfind(a.behavResources(k).SessionName,'Cond'))
        id_Cond(k) = 1;
    end
    if isfield(a.SessionEpoch, 'TestPost1')
        if ~isempty(strfind(a.behavResources(k).SessionName,'TestPost'))
            id_Post(k) = 1;
        end
    end
end
id_Pre=find(id_Pre);
id_Cond=find(id_Cond);
if isfield(a.SessionEpoch, 'TestPost1')
    id_Post=find(id_Post);
end

%% Calculate average occupancy
% Calculate occupancy de novo
for k=1:length(id_Pre)
    for t=1:length(a.behavResources(id_Pre(k)).Zone)
        Pre_Occup(k,t)=size(a.behavResources(id_Pre(k)).ZoneIndices{t},1)./...
            size(Data(a.behavResources(id_Pre(k)).Xtsd),1);
    end
end
for k=1:length(id_Cond)
    for t=1:length(a.behavResources(id_Cond(k)).Zone)
        Cond_Occup(k,t)=size(a.behavResources(id_Cond(k)).ZoneIndices{t},1)./...
            size(Data(a.behavResources(id_Cond(k)).Xtsd),1);
    end
end
Pre_Occup_Shock = squeeze(Pre_Occup(:,1));
Cond_Occup_Shock = squeeze(Cond_Occup(:,1));

Pre_Occup_Safe = squeeze(Pre_Occup(:,2));
Cond_Occup_Safe = squeeze(Cond_Occup(:,2));

Pre_Occup_Shock_mean = mean(Pre_Occup_Shock,2);
Pre_Occup_Shock_std = std(Pre_Occup_Shock,0,2);
Cond_Occup_Shock_mean = mean(Cond_Occup_Shock,2);
Cond_Occup_Shock_std = std(Cond_Occup_Shock,0,2);

Pre_Occup_Safe_mean = mean(Pre_Occup_Safe,2);
Pre_Occup_Safe_std = std(Pre_Occup_Safe,0,2);
Cond_Occup_Safe_mean = mean(Cond_Occup_Safe,2);
Cond_Occup_Safe_std = std(Cond_Occup_Safe,0,2);
    
%%% PostTests
if exist('id_Post', 'var')
    for k=1:length(id_Post)
        for t=1:length(a.behavResources(id_Post(k)).Zone)
            Post_Occup(k,t)=size(a.behavResources(id_Post(k)).ZoneIndices{t},1)./...
                size(Data(a.behavResources(id_Post(k)).Xtsd),1);
        end
    end
    Post_Occup_Shock = squeeze(Post_Occup(:,1));
    Post_Occup_Safe = squeeze(Post_Occup(:,2));
    Post_Occup_Shock_mean = mean(Post_Occup_Shock,2);
    Post_Occup_Shock_std = std(Post_Occup_Shock,0,2);
    Post_Occup_Safe_mean = mean(Post_Occup_Safe,2);
    Post_Occup_Safe_std = std(Post_Occup_Safe,0,2);
end



% Wilcoxon signed rank task between Pre and PostTest
% p_pre_post = signrank(Pre_Occup_Shock_mean, Post_Occup_Shock_mean);

if IsPCDriven
    tsStim = ts(Start(a.TTLInfo.StimEpoch));
end

%% Set up the figure
clrs = {[0.466 0.674 0.188], [0 0.447 0.741], [0.85 0.325 0.098], [0.494 0.184 0.556],...
    [0.635 0.078 0.184], [0.3010 .745 .933], [0.929 .694 .125], 'r'};

% Axes
fh = figure('units', 'normalized', 'outerposition', [0 0 0.9 0.9]);
PreTest_Axes = axes('position', [0.07 0.07 0.28 0.38]);
Cond_Axes = axes('position', [0.38 0.07 0.28 0.38]);
PostTest_Axes = axes('position', [0.69 0.07 0.28 0.38]);

PreTestTraj_Axes = axes('position', [0.07 0.54 0.28 0.38]);
CondTraj_Axes = axes('position', [0.38 0.54 0.28 0.38]);
PostTestTraj_Axes = axes('position', [0.69 0.54 0.28 0.38]);

maze = [0 0; 0 1; 1 1; 1 0; 0.63 0; 0.63 0.75; 0.35 0.75; 0.35 0; 0 0];
shockZone = [0 0; 0 0.43; 0.35 0.43; 0.35 0; 0 0]; 
safeZone = [1 0; 0.63 0; 0.63 0.43; 1 0.43; 1 0]; 

%%% 
if strcmp(experiment, 'PAG')
    colors = {[.8 0 0], [0 0 0.8]};
    titles = {'PreTests', 'Learning', 'PostTests'};
    xticklabels = {'Shock', 'Safe'};
elseif strcmp(experiment, 'MFB')
    colors = {[0 0.8 0], [.8 .9 0.8]};
    titles = {'PreTests', 'Learning', 'PostTests'};
    xticklabels = {'Stim', 'Neutral'};
elseif strcmp(experiment, 'Novel')    
    colors = {[.4 .8 .2], [0.793 .902 0.184]};
    titles = {'PreTests', 'LearningFake', 'PostTests'};
    xticklabels = {'Left', 'Right'};
end

%% Plot
% Trajectories
axes(PreTestTraj_Axes);
hold on
for i=1:numtest
    if IsColorTests
        plot(Data(a.behavResources(id_Pre(i)).AlignedXtsd),Data(a.behavResources(id_Pre(i)).AlignedYtsd),...
            'LineWidth',3, 'Color', clrs{i});
    else
        plot(Data(a.behavResources(id_Pre(i)).AlignedXtsd),Data(a.behavResources(id_Pre(i)).AlignedYtsd),...
            'LineWidth',3, 'Color', [0.8 0.8 0.8]);
    end
end
box off
set(gca,'XtickLabel',{},'YtickLabel',{});
plot(maze(:,1),maze(:,2),'k','LineWidth',3)
plot(shockZone(:,1),shockZone(:,2), 'color', colors{1}, 'LineWidth',3)
plot(safeZone(:,1),safeZone(:,2), 'color', colors{2}, 'LineWidth',3)
title(titles{1},'FontSize',18,'FontWeight','bold');
xlim([0 1])
ylim([0 1])
hold off
% set(gca, 'FontSize', 18, 'FontWeight',  'bold');

% Cond
axes(CondTraj_Axes);
hold on
for i=1:numtest
    if IsColorTests
        plot(Data(a.behavResources(id_Cond(i)).AlignedXtsd),Data(a.behavResources(id_Cond(i)).AlignedYtsd),...
            'LineWidth',3, 'Color', clrs{i});
    else
        plot(Data(a.behavResources(id_Cond(i)).AlignedXtsd),Data(a.behavResources(id_Cond(i)).AlignedYtsd),...
            'LineWidth',3, 'Color', [0.8 0.8 0.8]);
    end
    if strcmp(experiment, 'PAG') || strcmp(experiment, 'MFB')
        tempX = Data(a.behavResources(id_Cond(i)).AlignedXtsd);
        tempY = Data(a.behavResources(id_Cond(i)).AlignedYtsd);
        if IsPCDriven
            StimTimes_TTL = Restrict(tsStim, a.SessionEpoch.(['Cond' num2str(i)]));
            Times_MAT = ts(a.behavResources(id_Cond(i)).PosMat(:,1)*1e4);
            StimTimes_plot = Range(Restrict(Times_MAT,StimTimes_TTL, 'align', 'closest'));
            for j = 1:length(StimTimes_plot)
                plot(tempX(a.behavResources(id_Cond(i)).PosMat(:,1)==StimTimes_plot(j)/1e4),...
                    tempY(a.behavResources(id_Cond(i)).PosMat(:,1)==StimTimes_plot(j)/1e4),...
                    'p','Color','r','MarkerFaceColor','red','MarkerSize',16);
            end
        else
            plot(tempX(a.behavResources(id_Cond(i)).PosMat(:,4)==1),tempY(a.behavResources(id_Cond(i)).PosMat(:,4)==1),...
                'p','Color','r','MarkerFaceColor','red','MarkerSize',16);
        end
        clear tempX tempY
    end
end
box off
set(gca,'XtickLabel',{},'YtickLabel',{});
plot(maze(:,1),maze(:,2),'k','LineWidth',3)
plot(shockZone(:,1),shockZone(:,2), 'color', colors{1}, 'LineWidth',3)
plot(safeZone(:,1),safeZone(:,2), 'color', colors{2}, 'LineWidth',3)
title(titles{2},'FontSize',18,'FontWeight','bold');
xlim([0 1])
ylim([0 1])
hold off
% set(gca, 'FontSize', 18, 'FontWeight',  'bold');

% PostTests
axes(PostTestTraj_Axes);
if exist('id_Post', 'var') 
    hold on
    for i=1:numtest
        if IsColorTests
            plot(Data(a.behavResources(id_Post(i)).AlignedXtsd),Data(a.behavResources(id_Post(i)).AlignedYtsd),...
                'LineWidth',3, 'Color', clrs{i});
        else
            plot(Data(a.behavResources(id_Post(i)).AlignedXtsd),Data(a.behavResources(id_Post(i)).AlignedYtsd),...
                'LineWidth',3, 'Color', [0.8 0.8 0.8]);
        end
    end
    box off
    set(gca,'XtickLabel',{},'YtickLabel',{});
    plot(maze(:,1),maze(:,2),'k','LineWidth',3)
    plot(shockZone(:,1),shockZone(:,2), 'color', colors{1},'LineWidth',3)
    plot(safeZone(:,1),safeZone(:,2), 'color', colors{2}, 'LineWidth',3)
    title(titles{3},'FontSize',18,'FontWeight','bold');
    xlim([0 1])
    ylim([0 1])
    hold off
    % set(gca, 'FontSize', 18, 'FontWeight',  'bold');
else
    set(gca,'XtickLabel',{},'YtickLabel',{});
    text(.2,.4,'No PostTests!', 'FontWeight','bold','FontSize',24);
end

axes(PreTest_Axes)
[p_occ,h_occ, her_occ] = PlotErrorBarN_DB([Pre_Occup_Shock_mean*100 Pre_Occup_Safe_mean*100], 'barcolors', colors{1}, 'barwidth', 0.6, 'newfig', 0, 'showpoints',0);
h_occ.FaceColor = 'flat';
h_occ.FaceAlpha = .45;
h_occ.CData(2,:) = colors{2};
% set(gca,'Xtick',[1:2],'XtickLabel',{'PreTest', 'PostTest'});
set(gca, 'FontSize', 18, 'FontWeight',  'bold','FontName','Times New Roman');
set(gca, 'LineWidth', 5,'Xtick',[1:2],'XTickLabel', xticklabels);
set(h_occ, 'LineWidth', 3);
set(her_occ, 'LineWidth', 3);
line(xlim,[21.5 21.5],'Color','k','LineStyle','--','LineWidth',5);
% % text(1.85,23.2,'Random Occupancy', 'FontWeight','bold','FontSize',13);
% title('Percentage of the ShockZone occupancy', 'FontSize', 14);
ylabel('% occupancy','FontSize',18,'FontWeight','bold')
ylim([0 90])

axes(Cond_Axes)
[p_occ,h_occ, her_occ] = PlotErrorBarN_DB([Cond_Occup_Shock_mean*100 Cond_Occup_Safe_mean*100], 'barcolors', colors{1}, 'barwidth', 0.6, 'newfig', 0, 'showpoints',0);
h_occ.FaceColor = 'flat';
h_occ.FaceAlpha = .45;
h_occ.CData(2,:) = colors{2};
set(gca, 'FontSize', 18, 'FontWeight',  'bold','FontName','Times New Roman');
set(gca, 'LineWidth', 5,'Xtick',[1:2],'XTickLabel', xticklabels);
set(h_occ, 'LineWidth', 3);
set(her_occ, 'LineWidth', 3);
line(xlim,[21.5 21.5],'Color','k','LineStyle','--','LineWidth',5);
% text(1.85,23.2,'Random Occupancy', 'FontWeight','bold','FontSize',13);
% title('Percentage of the ShockZone occupancy', 'FontSize', 14);
ylim([0 90])

axes(PostTest_Axes)
if exist('id_Post', 'var')
    [p_occ,h_occ, her_occ] = PlotErrorBarN_DB([Post_Occup_Shock_mean*100 Post_Occup_Safe_mean*100], 'barcolors', colors{1}, 'barwidth', 0.6, 'newfig', 0, 'showpoints',0);
    h_occ.FaceColor = 'flat';
    h_occ.FaceAlpha = .45;
    h_occ.CData(2,:) = colors{2};
    set(gca, 'FontSize', 18, 'FontWeight',  'bold','FontName','Times New Roman');
    set(gca, 'LineWidth', 5,'Xtick',[1:2],'XTickLabel', xticklabels);
    set(h_occ, 'LineWidth', 3);
    set(her_occ, 'LineWidth', 3);
    line(xlim,[21.5 21.5],'Color','k','LineStyle','--','LineWidth',5);
    % text(1.85,23.2,'Random Occupancy', 'FontWeight','bold','FontSize',13);
    % title('Percentage of the ShockZone occupancy', 'FontSize', 14);
    ylim([0 90])
else
    set(gca,'XtickLabel',{},'YtickLabel',{});
    text(.2,.4,'No PostTests!', 'FontWeight','bold','FontSize',24);
end


%% Save it
if IsSave
    saveas(fh, [Dir '/Trajectories.fig']);
    saveFigure(fh,'Trajectories',Dir);
end

%% Write to xls file
% T = table(Pre_Occup_mean, Post_Occup_mean, Pre_entnum_mean, Post_entnum_mean,Pre_FirstTime_mean,Post_FirstTime_mean,...
%     Pre_VZmean_mean, Post_VZmean_mean);
% 
% filenme = [dir_out 'finalxls.xlsx'];
% writetable(T, filenme, 'Sheet',1,'Range','A1');
