function PlotSessionTrajectories(mice)
%
% This functions plot trajectories of a selected mouse in pretests,
% condiioning sessions and posttests as well as bar plots of shock-safe
%  zone occupancies
%
%  INPUT
%       
%       mice             numbers of mice from PAG exp
% 
% 
%  OUTPUT
%
%       Figure
%
%       See
%   
%       BehaviorERC, ExampleTrajectory_DB
% 
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 03/05/2021
% github.com/bryzgalovdm

%% Mazes
maze = [0 0; 0 1; 1 1; 1 0; 0.63 0; 0.63 0.75; 0.35 0.75; 0.35 0; 0 0];
shockZone = [0 0; 0 0.43; 0.35 0.43; 0.35 0; 0 0];
numtest = 4;

%% Get data
Dir = PathForExperimentsERC('UMazePAG');
Dir = RestrictPathForExperiment(Dir,'nMice', mice);

a = cell(length(Dir.path),1);
for imouse = 1:length(Dir.path)
    if strcmp(Dir.name{imouse}, 'Mouse905')
        a{imouse} = load([Dir.path{imouse}{1} '/behavResources_backup.mat'], 'behavResources');
    else
        a{imouse} = load([Dir.path{imouse}{1} '/behavResources.mat'], 'behavResources');
    end
end

%% Find necessary tests
id_Pre = cell(1,length(a));
id_Cond = cell(1,length(a));
id_Post = cell(1,length(a));

for i=1:length(a)
    id_Pre{i} = FindSessionID_ERC(a{i}.behavResources, 'TestPre');
    id_Cond{i} = FindSessionID_ERC(a{i}.behavResources, 'Cond');
    id_Post{i} = FindSessionID_ERC(a{i}.behavResources, 'TestPost');
end

%% Plot 
fh = figure('units', 'normalized', 'outerposition', [0 0 0.85 0.5]);
tabs = arrayfun(@(x) uitab('Title', Dir.name{x}), 1:length(Dir.path));
for itab = 1:length(tabs)
    curtab = tabs(itab);
    % PreTests
    ax(1) = axes('Parent', curtab, 'position', [0.03 0.07 0.3 0.85]);
    axes(ax(1));
    hold on
    for itest = 1:numtest
        plot(Data(a{itab}.behavResources(id_Pre{itab}(itest)).CleanAlignedXtsd),Data(a{itab}.behavResources(id_Pre{itab}(itest)).CleanAlignedYtsd),...
            'LineWidth',3, 'Color', [0.2 0.2 0.2]);
    end
    set(gca,'XtickLabel',{},'YtickLabel',{});
    hold on
    plot(maze(:,1),maze(:,2),'k','LineWidth',3)
    plot(shockZone(:,1),shockZone(:,2),'r','LineWidth',3)
    title('PreTests','FontSize',18,'FontWeight','bold');
    xlim([0 1])
    ylim([0 1])
    hold off
    makepretty
    % Cond
    ax(2) = axes('Parent', curtab, 'position', [0.36 0.07 0.3 0.85]);
    axes(ax(2));
    hold on
    for itest = 1:numtest
        plot(Data(a{itab}.behavResources(id_Cond{itab}(itest)).CleanAlignedXtsd),Data(a{itab}.behavResources(id_Cond{itab}(itest)).CleanAlignedYtsd),...
            'LineWidth',3, 'Color', [0.2 0.2 0.2]);
        tempX = Data(a{itab}.behavResources(id_Cond{itab}(itest)).CleanAlignedXtsd);
        tempY = Data(a{itab}.behavResources(id_Cond{itab}(itest)).CleanAlignedYtsd);
        plot(tempX(a{itab}.behavResources(id_Cond{itab}(itest)).PosMat(:,4)==1),tempY(a{itab}.behavResources(id_Cond{itab}(itest)).PosMat(:,4)==1),...
            'p','Color','r','MarkerFaceColor','red','MarkerSize',16);
        clear tempX tempY
    end
    set(gca,'XtickLabel',{},'YtickLabel',{});
    plot(maze(:,1),maze(:,2),'k','LineWidth',3)
    plot(shockZone(:,1),shockZone(:,2),'r','LineWidth',3)
    title('Conditioning','FontSize',18,'FontWeight','bold');
    xlim([0 1])
    ylim([0 1])
    hold off
    makepretty
    % PostTest
    ax(3) = axes('Parent', curtab, 'position', [0.69 0.07 0.3 0.85]);
    axes(ax(3));
    hold on
    for itest = 1:numtest
        plot(Data(a{itab}.behavResources(id_Post{itab}(itest)).CleanAlignedXtsd),Data(a{itab}.behavResources(id_Post{itab}(itest)).CleanAlignedYtsd),...
            'LineWidth',3, 'Color', [0.2 0.2 0.2]);
    end
    set(gca,'XtickLabel',{},'YtickLabel',{});
    plot(maze(:,1),maze(:,2),'k','LineWidth',3)
    plot(shockZone(:,1),shockZone(:,2),'r','LineWidth',3)
    title('PostTests','FontSize',18,'FontWeight','bold');
    xlim([0 1])
    ylim([0 1])
    hold off
    makepretty
end

end