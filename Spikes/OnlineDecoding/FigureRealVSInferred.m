function f1 = FigureRealVSInferred(folder, thresh, varargin)
% 
% This funciton will plot the figure with two subplots:
% - Real position in test set
% - Inferred position in the test set
% Inferred position are filtered acording to the evaluated loss (confidence of the network)
% 
%  INPUT
% 
% 
%         folder          folder with inferring.mat
%         thresh          threshold on evaluated loss (data LESS than <thresh> are taken)
% 
% 
%         Optional:
% 
%         SzSide          side on which to-be-stimulated arm is located ('left' or 'right')
% 
%         
%  OUTPUT
%  
%         fh              figure handle
% 
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 23/06/2020; Extended 21/08/2020
% github.com/bryzgalovdm

%% Defaults
szside = 'left';

%% Arguments management
for i=1:2:length(varargin)    
    switch(lower(varargin{i})) 
        case 'szside'
            szside = varargin{i+1};
            if ~isa(szside,'char')
                error('Incorrect value for property ''SZside'' (type ''help FigureRealVSInferred'' for details).');
            elseif ~strcmp(szside, 'right')
                if ~strcmp(szside, 'left')
                    error('Shock zone side could be either left or right (type ''help FigureRealVSInferred'' for details).');
                end
            elseif ~strcmp(szside, 'left')
                if ~strcmp(szside, 'right')
                    error('Shock zone side could be either left or right (type ''help FigureRealVSInferred'' for details).');
                end
            end
    end
end

%% Load data
cd(folder);
file = load('inferring.mat');
if isunix
    cd ../
else
    idcs = strfind(folder, '/');
    if sum(idcs == length(folder)) == 0
        cd(folder(1:idcs(end)));
    else
        cd(folder(1:idcs(end-1)));
    end
end
orig = load('nnBehavior.mat', 'behavior');

%% Do a folder
if exist('RealVsInferred_figs', 'dir') ~= 7
    mkdir('RealVsInferred_figs');
end

%% Initialize the maze
maze = [0.36 0.38; 0.36 .96; 1.04 .96; 1.04 0.38; 0.78 0.38; 0.78 0.80; 0.62 0.80; 0.62 0.38; 0.36 0.38];
switch szside
    case 'right'
        minx=.8;
        maxx=1;
        miny=0;
        maxy=.55;
        ShockZone = [.78 .38 .26 .21];
        
    case 'left'
        minx=.35;
        maxx=.6;
        miny=.38;
        maxy=.58;
        ShockZone = [.36 .38 .26 .21];
end

%% Loop over thresholds
for itresh = 1:length(thresh)
    % Remove all the points with evaluated loss lower than thresh
    inferringPlot = file.inferring(file.inferring(:,3)<thresh(itresh),:);
    inferredposX = inferringPlot(:,1);
    inferredposY = inferringPlot(:,2);
    
    timePlot = file.times(file.inferring(:,3)<thresh(itresh));
    
    posPlot = file.pos(file.inferring(:,3)<thresh(itresh),:);
    posx = posPlot(:,1);
    posy = posPlot(:,2);
    % eliminate zeros
    posx(inferredposX==0)=[];
    posy(inferredposY==0)=[];
    timePlot(inferredposY==0)=[];
    inferredposX(inferredposX==0)=[];
    inferredposY(inferredposY==0)=[];

    inferredposX(posx==0)=[];
    inferredposY(posy==0)=[];
    timePlot(posx==0)=[];
    posx(posx==0)=[];
    posy(posy==0)=[];
    

    maxt = length(posx);
    tps=timePlot-file.times(1);
    
    % restrict to shock zone position
    for i=1:length(inferredposX)
        idx=find(inferredposX>minx & inferredposX<maxx);
        idy=find(inferredposY>miny & inferredposY<maxy);
        val = intersect(idx,idy);
    end
    
    % Figure RealVsInferred
    f1 = figure('rend','painters','pos',[400+itresh*10 400+itresh*10 900 600]);
    hold on
    plot(posx(1:maxt),posy(1:maxt),'ko','markerfacecolor','k')
    line([posx(1:maxt) inferredposX(1:maxt)]', [posy(1:maxt) inferredposY(1:maxt)]','color',[0.7 0.7 0.7])
    plot(inferredposX(1:maxt), inferredposY(1:maxt),'k.')
    scatter(posx,posy,15,tps,'filled')
    plot(maze(:,1),maze(:,2),'k','LineWidth',2)
    rectangle('Position',ShockZone,'EdgeColor','r','LineWidth',2)
    title(['Thresh = ' num2str(thresh(itresh))])
    set(gca, 'XTickLabel', {}, 'YTickLabel', {});
    set(gca, 'FontSize', 16, 'FontWeight', 'bold');
    legend('Real pos', 'Inference')
    
    saveas(f1,[pwd '/RealVsInferred_figs/AllPoints_' szside '_' num2str(thresh(itresh)) '.fig']);
    saveFigure(f1, ['AllPoints_' szside '_' num2str(thresh(itresh))],[pwd '/RealVsInferred_figs/']);
    
    % Figure inferred and real
    f2 = figure('Color',[1 1 1],'rend','painters','pos',[100+itresh*10 20+itresh*10 1600 600])
    subplot(121)
    scatter(posx(val),posy(val),50,tps(val)','filled');
    title('Real position', 'FontSize', 16)
    hold on
    plot(maze(:,1),maze(:,2),'k','LineWidth',2)
    rectangle('Position',ShockZone,'EdgeColor','r','LineWidth',2)
    xlim([0.335 1.045]);
    ylim([0.36 1])
    set(gca,'visible','off')
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    set(gca, 'XTickLabel', {}, 'YTickLabel', {});
    set(gca, 'FontSize', 16, 'FontWeight', 'bold');
    subplot(122)
    scatter(inferredposX(val),inferredposY(val),50,tps(val)','filled');
    title('Online inference', 'FontSize', 16)
    colorbar
    hold on
    plot(maze(:,1),maze(:,2),'k','LineWidth',2)
    rectangle('Position',ShockZone,'EdgeColor','r','LineWidth',2)
    xlim([0.335 1.045])
    ylim([0.36 1])
    set(gca,'visible','off')
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    set(gca, 'XTickLabel', {}, 'YTickLabel', {});
    set(gca, 'FontSize', 16, 'FontWeight', 'bold');
    annotation('textbox', [0.46 0.86 0.1 0.1], 'String', ['Thresh = ' num2str(thresh(itresh))], 'FontWeight', 'bold', 'FontSize', 12,...
        'EdgeColor', 'none', 'FontSize', 16, 'FontWeight', 'bold')
    annotation('textbox', [0.45 0.03 0.1 0.1], 'String', [num2str(length(val)) ' points decoded'], 'FontWeight', 'bold', 'FontSize', 12,...
        'EdgeColor', 'none', 'FontSize', 16, 'FontWeight', 'bold')
    
    saveas(f2,[pwd '/RealVsInferred_figs/RealVsInf_' szside '_' num2str(thresh(itresh)) '.fig']);
    saveFigure(f2, ['RealVsInf_' szside '_' num2str(thresh(itresh))],[pwd '/RealVsInferred_figs/']);
end

