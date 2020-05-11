function MakeViolinPlot_DB(A,Cols,X,Legends,ShowPoints, varargin)
% Input variables
% A = your data in a cell format A = {Data1,Data2,Data3}
% Cols = The color for each set of data Cols = {[1 0 0],[0 0 1],[0 1 0]}
% If you don't care about colors just put {} and everything wille be grey
% X = the position to plot your data X = [1,2,3]
% If you don't care just put []
% Legends = the identity of your datasets for xlabels Legends = {'bla' 'bla' 'bla'}
% If you don't care about colors just put {}
% ShowPoints = 0 for no points = 1 for points
%
% Modified from MakeSpreadAndBoxPlot_SB.m

if isempty(Cols)
    for i = 1:length(A)
        Cols{i} = [0.6 0.6 0.6];
    end
end

if isempty(X)
    for i = 1:length(A)
        X(i) = i*2;
    end
end

if isempty(Legends)
    for i = 1:length(A)
        Legends{i} = num2str(i);
    end
end

for k = 1 : length(A)
    if sum(isnan(A{k}))<length(A{k})
        a=iosr.statistics.boxPlot(X(k),A{k}(:),'boxColor', Cols{k}, 'boxAlpha', 0, 'lineColor',[0.95 0.95 0.95],...
            'medianColor','k','boxWidth','auto', 'LineColor', 'k', 'LineWidth', 3, 'showOutliers',false,'showViolin',true,'violinColor', Cols{k});
        a.handles.upperWhiskers.Visible='on';a.handles.upperWhiskerTips.Visible='on';
        a.handles.lowerWhiskers.Visible='on';a.handles.lowerWhiskerTips.Visible='on';
        a.handles.medianLines.LineWidth = 5;
        hold on
        if ShowPoints
            handlesplot=plotSpread(A{k}(:),'distributionColors','k','xValues',X(k),'spreadWidth',0.8); hold on;
            set(handlesplot{1},'MarkerSize',20)
            handlesplot=plotSpread(A{k}(:),'distributionColors',Cols{k}*0.4,'xValues',X(k),'spreadWidth',0.8); hold on;
            set(handlesplot{1},'MarkerSize',10)
        end
        MinA(k) = min(A{k});
        MaxA(k) = max(A{k});
    end
    
end
xlim([min(X)-1 max(X)+1])
rg = max(MaxA)-min(MinA);
if ShowPoints
    ylim([min(MinA)-rg*0.1 max(MaxA)+rg*0.1])
end
if exist('Legends')
    set(gca,'XTick',X,'XTickLabel',Legends)
    box off
    set(gca,'FontSize',10,'Linewidth',1)
end