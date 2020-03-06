%PlotBasicSpikeData - Plot basic properties of sorted clusters.
%
% Plot the basic properties obtained from <CalcBasicNeuronInfo> function
%
%  USAGE
%
%    f = PlotBasicSpikeData(BasicNeuronInfo, saveplot)
%
%    BasicNeuronInfo       Basic information (see CalcBasicNeuronInfo.m)
%    saveplot     if 1 will save the plot; 0 - not      
%
%  OUTPUT
%
%    f      Handle to the figure
%
%       See
%   
%       CalcBasicNeuronInfo
% 
% Copyright (C) 2018 by Dmitri Bryzgalov
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

%% TODO: 4Hz modulation - and correct it - put imagesc figure of phases for each neuron

function f = PlotBasicSpikeData(BasicNeuronInfo, Quality, saveplot)

%% Prepare
% Load data


% Calculate number of MUA and SUA
num_MUA = length(BasicNeuronInfo.idx_MUA);
num_SUA = length(BasicNeuronInfo.idx_SUA);

% Calculate number in both classes
AllPyr = sum(BasicNeuronInfo.neuroclass(BasicNeuronInfo.idx_SUA)>0);
id_Pyr = find(BasicNeuronInfo.neuroclass(BasicNeuronInfo.idx_SUA)>0);
AllInt = sum(BasicNeuronInfo.neuroclass(BasicNeuronInfo.idx_SUA)<0);
id_Int = find(BasicNeuronInfo.neuroclass(BasicNeuronInfo.idx_SUA)<0);

% Plot
f = figure('units', 'normalized', 'outerposition', [0 0 1 1]); % units normalized and stuff
PyrInt_Axes = axes('position', [0.03 0.73 0.25 0.25]);
LRatio_Axes = axes('position', [0.03 0.55 0.15 0.15]);
IsoD_Axes = axes('position', [0.21 0.55 0.15 0.15]);
SubQ_Axes = axes('position', [0.39 0.55 0.15 0.15]);

CompareFR_Axes = axes('position', [0.04 0.07 0.27 0.37]);
OverallFR_Axes = axes('position', [0.28 0.755 0.22 0.17]);
FourHz_Axes = axes('position', [0.35 0.325 0.2 0.18]);
FourHzPie_Axes = axes('position', [0.48 0.425 0.07 0.07]);
EnvelopeFourHz_Axes = axes('position', [0.35 0.25 0.2 0.03]);
PhaseNeuronsFourHz_Axes = axes('position', [0.35 0.06 0.2 0.18]);
% old
% WakeFR_Axes = axes('position', [0.3 0.3 0.22 0.17]);
% SWSFR_Axes = axes('position', [0.03 0.07 0.22 0.17]);
% REMFR_Axes = axes('position', [0.3 0.07 0.22 0.17]);
% old

Kappa_Axes = axes('position', [0.55 0.8 0.2 0.17]);
ModPye_Axes = axes('position', [0.68 0.65 0.07 0.07]);
CumZ_Axes = axes('position', [0.58 0.50 0.17 0.22]);
Phase_Axes = axes('position', [0.8 0.78 0.18 0.18]);
EnvelopeTheta_Axes = axes('position', [0.78 0.71 0.2 0.04]);
PhaseNeuronsTheta_Axes = axes('position', [0.78 0.5 0.2 0.2])

RippleMod = axes('position', [0.59 0.33 0.37 0.12]);
RippleModInd = axes('position', [0.59 0.06 0.37 0.265]);

%% Plot pie Pyr/int
axes(PyrInt_Axes);
pieid  = pie([AllPyr AllInt]);
legend ({'Principal cells', 'Interneurons'}, 'Position', [0.28 0.96 0.17 0.035]);
pieid(1).FaceColor = [0 0.484 0.87];
pieid(3).FaceColor = [0.85 0.094 0.082];
title(['MUA/SUA = ', sprintf('%2i', num_MUA), '/', sprintf('%2i', num_SUA),...
    ', N of PC = ', sprintf('%3i',AllPyr), ', N of Int = ', sprintf('%3i',AllInt)]);

%% Plot Isolation Distance of all SUA
axes(IsoD_Axes);
IsoDistance_SUA = Quality.IsoDistance(BasicNeuronInfo.idx_SUA);
fisod = nhist(IsoDistance_SUA,'minbins', 25, 'color', [0 0 0], 'noerror', 'stdtimes', 6);
title('Isolation distances of all SUA');

%% Plot LRatio of all SUA
axes(LRatio_Axes);
LRatio_SUA = Quality.LRatio(BasicNeuronInfo.idx_SUA);
flra = nhist(LRatio_SUA,'minbins', 25, 'color', [0 0 0], 'noerror', 'stdtimes', 6);
title('LRatios of all SUA');

%% Plot subjective quality of all SUA
% axes(SubQ_Axes);
% SubjectiveQuality_SUA = str2double(Quality.MyMark(BasicNeuronInfo.idx_SUA));
% fsubq = nhist(SubjectiveQuality_SUA, 'color', [0 0 0], 'noerror', 'stdtimes', 6);
% title('Subjective qualities of all SUA');

%% Plot histogram of firingrates
axes(OverallFR_Axes); %%% Pyr/Int
firingrates_SUA = BasicNeuronInfo.firingrate(BasicNeuronInfo.idx_SUA);
frid = nhist(firingrates_SUA, 'minbins', 25, 'color', [0 0 0], 'noerror', 'stdtimes', 6); 
title('Firing rates of all SUA');

% Old
% axes(CompareFR_Axes); %% Pyr/Int
% FR_tocomp = [mean(BasicNeuronInfo.FR_WakeTheta(id_Pyr)) mean(BasicNeuronInfo.FR_WakeTheta(id_Int));...
%     mean(BasicNeuronInfo.FR_WakeNoTheta(id_Pyr)) mean(BasicNeuronInfo.FR_WakeNoTheta(id_Int));...
%     mean(BasicNeuronInfo.FR_SWS(id_Pyr)) mean(BasicNeuronInfo.FR_SWS(id_Int));...
%     mean(BasicNeuronInfo.FR_REM(id_Pyr)) mean(BasicNeuronInfo.FR_REM(id_Int))];
% FR_std_up = [std(BasicNeuronInfo.FR_WakeTheta(id_Pyr)) std(BasicNeuronInfo.FR_WakeTheta(id_Int));...
%     std(BasicNeuronInfo.FR_WakeNoTheta(id_Pyr)) std(BasicNeuronInfo.FR_WakeNoTheta(id_Int));...
%     std(BasicNeuronInfo.FR_SWS(id_Pyr)) std(BasicNeuronInfo.FR_SWS(id_Int));...
%     std(BasicNeuronInfo.FR_REM(id_Pyr)) std(BasicNeuronInfo.FR_REM(id_Int))];
% FR_std_down = [0 0; 0 0; 0 0; 0 0];
% FR_std = cat(3,FR_std_down,FR_std_up);
% [ffr,ffrer] = barwitherr(FR_std,FR_tocomp);
% set(gca,'Xtick',[1:4],'XtickLabel',{'WakeTheta','WakeNoTheta','SWS','REM'});
% set(gca, 'FontSize', 12, 'FontWeight',  'bold');
% ffr(1).CData = [0 0 1];
% ffr(2).CData = [1 0 0];
% set(ffr, 'LineWidth', 3);
% set(ffrer, 'LineWidth', 3);
% ylabel('Mean firing rate (Hz)');
% title('Firing rate of all SUA in different states', 'FontSize', 14);
% legend('Pyr', 'Int', 'Location', 'NorthEast')
% % [p_fr,h_fr, her_fr] = PlotErrorBarN_DB(FR_tocomp, 'barcolors', [0 0 0], 'barwidth', 0.6, 'newfig', 0,'showpoints',0);

% Box plot
FR_tocomp = {BasicNeuronInfo.FR_WakeTheta(id_Pyr); BasicNeuronInfo.FR_WakeTheta(id_Int);...
    BasicNeuronInfo.FR_WakeNoTheta(id_Pyr); BasicNeuronInfo.FR_WakeNoTheta(id_Int);...
    BasicNeuronInfo.FR_SWS(id_Pyr); BasicNeuronInfo.FR_SWS(id_Int);...
    BasicNeuronInfo.FR_REM(id_Pyr); BasicNeuronInfo.FR_REM(id_Int)};
% FR_std_up = [std(BasicNeuronInfo.FR_WakeTheta(id_Pyr)) std(BasicNeuronInfo.FR_WakeTheta(id_Int));...
%     std(BasicNeuronInfo.FR_WakeNoTheta(id_Pyr)) std(BasicNeuronInfo.FR_WakeNoTheta(id_Int));...
%     std(BasicNeuronInfo.FR_SWS(id_Pyr)) std(BasicNeuronInfo.FR_SWS(id_Int));...
%     std(BasicNeuronInfo.FR_REM(id_Pyr)) std(BasicNeuronInfo.FR_REM(id_Int))];
% FR_std_down = [0 0; 0 0; 0 0; 0 0];
% FR_std = cat(3,FR_std_down,FR_std_up);
Cols = {[0 0 0.9], [0.9 0 0], [0 0 0.9], [0.9 0 0],[0 0 0.9], [0.9 0 0],[0 0 0.9], [0.9 0 0]};

addpath(genpath('/home/mobsrick/Dima/MatlabToolbox-master/'));

axes(CompareFR_Axes); %% Pyr/Int
MakeSpreadAndBoxPlot_SB(FR_tocomp,Cols,[1:8]);
set(gca,'LineWidth',3,'FontWeight','bold','FontSize',10,'XTick',1.5:2:7.5,...
    'XTickLabel',{'WakeTheta','WakeNoTheta','SWS','REM'})
ylabel('Mean firing rate (Hz)', 'FontWeight','bold','FontSize',14);
title('Firing rate of all SUA in different states', 'FontSize', 14);
% legend('Pyr', 'Int', 'Location', 'NorthEast')

rmpath(genpath('/home/mobsrick/Dima/MatlabToolbox-master/'));
% [p_fr,h_fr, her_fr] = PlotErrorBarN_DB(FR_tocomp, 'barcolors', [0 0 0], 'barwidth', 0.6, 'newfig', 0,'showpoints',0);

%% old
% FRWake_SUA = BasicNeuronInfo.FRWake(BasicNeuronInfo.idx_SUA);
% frwid = nhist(FRWake_SUA, 'minbins', 25, 'color', [0 0.4 0.4], 'noerror', 'stdtimes', 6); 
% title('Firing rates of all SUA during wake');for i=1:num_SUA

% 
% axes(SWSFR_Axes);
% FR_SWS_SUA = BasicNeuronInfo.FR_SWS(BasicNeuronInfo.idx_SUA);
% frsid = nhist(FR_SWS_SUA, 'minbins', 25, 'color', [0 0.4 0.4], 'noerror', 'stdtimes', 6); 
% title('Firing rates of all SUA during SWS');
% 
% axes(REMFR_Axes);
% FR_REM_SUA = BasicNeuronInfo.FR_REM(BasicNeuronInfo.idx_SUA);
% frrid = nhist(FR_REM_SUA, 'minbins', 25, 'color', [0 0.4 0.4], 'noerror', 'stdtimes', 6); 
% title('Firing rates of all SUA during REM');

%% Plot modulation theta
axes(Kappa_Axes);
fdkap = nhist(BasicNeuronInfo.kappatheta.Transf,'minbins', 25, 'color', [0 0 0], 'noerror', 'stdtimes', 6);
title('Distribution of kappa: Neurons modulated by theta');


perc_sign = length(find(log(BasicNeuronInfo.Z.Transf(BasicNeuronInfo.idx_SUA)) > log(-log(0.05))))/num_SUA*100;
perc_nonsign = 100-perc_sign;

% Plot cumulative distribution
axes(CumZ_Axes); %% Pyr/Int
[Y1,X1]=(hist(log(BasicNeuronInfo.Z.Transf(id_Pyr)),100));
plot(X1,100*(1-cumsum(Y1)/sum(Y1)),'b','linewidth',2)
xlim([0 max(log((BasicNeuronInfo.Z.Transf)))]);
hold on
[Y2,X2]=(hist(log(BasicNeuronInfo.Z.Transf(id_Int)),100));
plot(X2,100*(1-cumsum(Y2)/sum(Y2)),'r','linewidth',2)
hold on
line([log(-log(0.05)) log(-log(0.05))],ylim,'linestyle','--','color','k')
legend('Pyr', 'Int', 'Location', 'SouthWest')
xlabel('ln(Z)')
ylabel('% Neurons')
title([num2str(round(perc_sign)), '% neurons theta-modulated']);

% Plot pie graph
axes(ModPye_Axes);
pieid  = pie([perc_sign perc_nonsign]);
if perc_sign<100
    pieid(1).FaceColor = [1 0 0];
    pieid(3).FaceColor = [1 1 1];
else
    pieid(1).FaceColor = [1 0 0];
end


% Plot rose graph
axes(Phase_Axes);
allphasetheta = [];
for i=1:num_SUA
    mu = BasicNeuronInfo.mutheta.Transf(BasicNeuronInfo.idx_SUA(i));
    allphasetheta = [allphasetheta; mu];
end
rose(allphasetheta);
title('Angles of theta-modulated SUA neurons')
    
% To plot imagesc of phases - how to sort them??? - according to mu
axes(PhaseNeuronsTheta_Axes)
for i=1:num_SUA
    n(i,1:60) = zscore(histcounts(BasicNeuronInfo.phasestheta.Transf{i}, 60));
end
[a,phasesid] = sort(BasicNeuronInfo.mutheta.Transf(BasicNeuronInfo.idx_SUA)); % sorted by mu
phaseindiv = n(phasesid,:);
imagesc(linspace(0,2*pi,60), 1:num_SUA, phaseindiv)
set(gca,'XTick',[]); set(gca,'XTickLabel',[]);
xticks([0 pi/2 pi 1.5*pi 2*pi])
xticklabels({'0', '90', '180', '270', '360'})
ylabel('#Neuron')
xlabel('Phase in deg')

axes(EnvelopeTheta_Axes);
t = linspace(0,2*pi,60);
e = cos(linspace(0,2*pi,60));
plot(t,e, 'k','linewidth', 2);
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
xlim([0 2*pi])
ylim([-1.5 1.5])
title('Theta modulation - neurons')

%% Plot modulation 4Hz
if isfield(BasicNeuronInfo,'Z4Hz')
    perc_sign = length(find(log(BasicNeuronInfo.Z4Hz.Transf(BasicNeuronInfo.idx_SUA)) > log(-log(0.05))))/num_SUA*100;
    perc_nonsign = 100-perc_sign;
    
    % Plot cumulative distribution
    axes(FourHz_Axes); %% Pyr/Int
    [Y1,X1]=(hist(log(BasicNeuronInfo.Z4Hz.Transf(id_Pyr)),100));
    plot(X1,100*(1-cumsum(Y1)/sum(Y1)),'b','linewidth',2)
    xlim([0 max(log((BasicNeuronInfo.Z4Hz.Transf)))]);
    hold on
    [Y2,X2]=(hist(log(BasicNeuronInfo.Z4Hz.Transf(id_Int)),100));
    plot(X2,100*(1-cumsum(Y2)/sum(Y2)),'r','linewidth',2)
    hold on
    line([log(-log(0.05)) log(-log(0.05))],ylim,'linestyle','--','color','k')
    legend('Pyr', 'Int', 'Location', 'SouthWest')
    xlabel('ln(Z)')
    ylabel('% Neurons')
    title([num2str(round(perc_sign)), '% neurons 4Hz-modulated']);
    
    % Plot pie graph
    axes(FourHzPie_Axes);
    pieid  = pie([perc_sign perc_nonsign]);
    if perc_sign<100
        pieid(1).FaceColor = [1 0 0];
        pieid(3).FaceColor = [1 1 1];
    else
        pieid(1).FaceColor = [1 0 0];
    end
    
    % To plot imagesc of phases - how to sort them???
    axes(PhaseNeuronsFourHz_Axes);
    for i=1:num_SUA
        n(i,1:60) = zscore(histcounts(BasicNeuronInfo.phases4Hz.Transf{i}, 60));
    end
    [a,phasesid] = sort(BasicNeuronInfo.mu4Hz.Transf(BasicNeuronInfo.idx_SUA)); % sorted by mu
    phaseindiv = n(phasesid,:);
    imagesc(linspace(0,2*pi,60), 1:num_SUA, phaseindiv)
    set(gca,'XTick',[]); set(gca,'XTickLabel',[]);
    xticks([0 pi/2 pi 1.5*pi 2*pi])
    xticklabels({'0', '90', '180', '270', '360'})
    ylabel('#Neuron')
    xlabel('Phase in deg')
    
    axes(EnvelopeFourHz_Axes);
    t = linspace(0,2*pi,60);
    e = cos(linspace(0,2*pi,60));
    plot(t,e, 'k','linewidth', 2);
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    xlim([0 2*pi])
    ylim([-1.5 1.5])
end
%% Plot response to ripples
axes(RippleMod);
plot(BasicNeuronInfo.BT_neurip(BasicNeuronInfo.idx_SUA,:),...
    runmean(mean(zscore((BasicNeuronInfo.CC_neurip(BasicNeuronInfo.idx_SUA,:))')'),1),'k','linewidth',3)
title('Modulation by ripples - SUA');
set(gca,'XTickLabel', []);

axes(RippleModInd);
imagesc(BasicNeuronInfo.BT_neurip(1,:),1:num_SUA,zscore(BasicNeuronInfo.CC_neurip')');
xlabel('Time');
ylabel('Neurons');
%% Save plot
if saveplot == 1
    res = pwd;
    saveas(gcf, [res '/BasicNeuronInfo.fig']);
    saveFigure(gcf,'BasicNeuronInfo',res);
end

%% Close the figure
% close(f);

end