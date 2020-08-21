%%%
%
% TODO Incoroporate two things:
%  - Simple thresholding
%  - Stimulation simulation according to EL and time Thresholding
%  - REmove hardcoded stuff from auxiliary functions
%  - Clean everything
% 
%  - MAIN THING IS TH DELAY  - IN DEVELOPEMENT RIGHT NOW

%% Initialize the maze
mazePosMat = [22 21; 22 52; 58 52; 58 21; 45 21; 45 45; 35 45; 35 21; 22 21];
ShockZonePosMat = [58 21; 45 21; 45 35; 58 35; 58 21];
mazeDecoded = [0.4 0.42; 0.4 0.93; 1 0.93; 1 0.42; 0.8 0.42; 0.8 0.82; 0.6 0.82; 0.6 0.42; 0.4 0.42];
shockZoneDecoded = [0.8 0.42; 0.8 0.65; 1 0.65; 1 0.42; 0.8 0.42];

%% Load
load('behavResources.mat');

%% Prepare array
StimData = DecodedData(DecodedStimIDX,:);
TimeStim = ts(double(StimData(:,1)*1e4));
% RealStimTime = ts(double(Start(StimEpoch)));

OEStim = DecodedData(DecodedStimIDX,1)*1e4;
RealStimTime = double(Start(StimEpoch));

%% Find one to one correspondance of HW timestamps and stimulation instances
ITimes = nan(length(RealStimTime),1);
for i = 1:length(RealStimTime)
    tempdd = RealStimTime(i) - OEStim;
    [~,id] = min(tempdd(tempdd>0));  % Find closest preceding value
    ITimes(i) = OEStim(id);
end
ToCheck = [ITimes RealStimTime];

%% Histogram of delays

dd = ToCheck(:,2) - ToCheck(:,1);

f2 = figure('Units', 'normalized', 'OuterPosition', [0.3 0.5 0.4 0.5]);
h = histogram(dd/10,25);
h.FaceColor = 'k';
set(gca, 'FontWeight', 'bold', 'FontSize', 14)
title('Delay between location decoded and following stimulation', 'FontSize', 16)
xlabel('Delay (ms)')
xline(36, 'Color', 'r', 'Linewidth', 2)
legend(['N=' num2str(length(dd)) ' stims'], '36 ms')
% 
% saveas(f2, '/home/mobsrick/Dropbox/MOBS_workingON/Dima/Ongoing_results/OnlineDecoding/August2020/delays_2h.fig');
% saveFigure(f2, 'delays_2h', '/home/mobsrick/Dropbox/MOBS_workingON/Dima/Ongoing_results/OnlineDecoding/August2020/');

%% Histogram of error
f1 = figure('Units', 'normalized', 'OuterPosition', [0.3 0.5 0.4 0.5]);
h = histogram(DecodedData(:,4), 50);
h.FaceColor = 'k';
line([0.015 0.015], ylim, 'Color', 'r', 'LineWidth', 2)
line([0.01 0.01], ylim, 'Color', 'b', 'LineWidth', 2)
line([0.005 0.005], ylim, 'Color', 'g', 'LineWidth', 2)
set(gca, 'FontWeight', 'bold', 'FontSize', 14)
title('Network confidence distribution', 'FontSize', 16)

% saveas(f1, '/home/mobsrick/Dropbox/MOBS_workingON/Dima/Ongoing_results/OnlineDecoding/July2020/dist_EL.fig');
% saveFigure(f1, 'dist_EL', '/home/mobsrick/Dropbox/MOBS_workingON/Dima/Ongoing_results/OnlineDecoding/July2020/');

%% Static real vs inferred
f2 = figure('Units', 'normalized', 'OuterPosition', [0.1 0.3 0.8 0.6])
subplot(121)
scatter(PosMat(:,3), PosMat(:,2), 14, PosMat(:,1));
hold on
scatter(PosMat(1,3), PosMat(1,2), 72, PosMat(1,1), 'MarkerFaceColor', 'r');
xlim([20 60]);
ylim([20 53]);
plot(mazePosMat(:,1),mazePosMat(:,2),'k','LineWidth',3)
plot(ShockZonePosMat(:,1),ShockZonePosMat(:,2),'r','LineWidth',3)
hold off
title('Real position', 'FontSize', 16)
set(gca, 'XTickLabel', {}, 'YTickLabel', {});
subplot(122)
scatter(DecodedData(DecodedData(:,4)<0.015,2),DecodedData(DecodedData(:,4)<0.015,3), 14, DecodedData(DecodedData(:,4)<0.015,1));
xlim([0.37 1]);
ylim([0.39 0.95])
title('Online inference', 'FontSize', 16)
hold on
plot(mazeDecoded(:,1),mazeDecoded(:,2),'k','LineWidth',3)
plot(shockZoneDecoded(:,1),shockZoneDecoded(:,2),'r','LineWidth',3)
hold off
c = colorbar;
set(get(c, 'title'), 'string', 'Time (s)', 'FontSize', 14)
set(gca, 'XTickLabel', {}, 'YTickLabel', {});

% saveas(f2, '/home/mobsrick/Dropbox/MOBS_workingON/Dima/Ongoing_results/OnlineDecoding/July2020/RealVsInferred.fig');
% saveFigure(f2, 'RealVsInferred', '/home/mobsrick/Dropbox/MOBS_workingON/Dima/Ongoing_results/OnlineDecoding/July2020/');

%%
[stimtimes1, stim_idx1] = SimulateStim_PositionDecoder(DecodedData, 0.015, 1);
[stimtimes2, stim_idx2] = SimulateStim_PositionDecoder(DecodedData, 0.01, 1);
[stimtimes3, stim_idx3] = SimulateStim_PositionDecoder(DecodedData, 0.005, 1);


%% Location of stim vs real position
f3 = figure('Units', 'normalized', 'OuterPosition', [0 0 0.5 1])
% Real data
subplot(421)
scatter(Data(Restrict(Xtsd, TimeStim, 'align', 'closest')), Data(Restrict(Ytsd, TimeStim, 'align', 'closest')));
hold on
plot(mazePosMat(:,1),mazePosMat(:,2),'k','LineWidth',3)
plot(ShockZonePosMat(:,1),ShockZonePosMat(:,2),'r','LineWidth',3)
hold off
xlim([20 60]);
ylim([20 53]);
set(gca, 'XTickLabel', {}, 'YTickLabel', {});
title('Real position', 'FontSize', 16)
subplot(422)
scatter(StimData(:,2), StimData(:,3))
hold on
plot(mazeDecoded(:,1),mazeDecoded(:,2),'k','LineWidth',3)
plot(shockZoneDecoded(:,1),shockZoneDecoded(:,2),'r','LineWidth',3)
hold off
xlim([0.37 1]);
ylim([0.39 0.95])
set(gca, 'XTickLabel', {}, 'YTickLabel', {});
title('Online inference', 'FontSize', 16)
annotation('textbox', [0.4 0.87 0.1 0.1], 'String', 'Thresh = 0.015 (Real)', 'FontWeight', 'bold', 'FontSize', 12,...
    'EdgeColor', 'none')

% Threshold 1
TimeStim1 = ts(double(stimtimes1*1e4));
subplot(423)
scatter(Data(Restrict(Xtsd, TimeStim1, 'align', 'closest')), Data(Restrict(Ytsd, TimeStim1, 'align', 'closest')));
xlim([20 60]);
ylim([20 53]);
hold on
plot(mazePosMat(:,1),mazePosMat(:,2),'k','LineWidth',3)
plot(ShockZonePosMat(:,1),ShockZonePosMat(:,2),'r','LineWidth',3)
hold off
set(gca, 'XTickLabel', {}, 'YTickLabel', {});
subplot(424)
scatter(DecodedData(stim_idx1,2), DecodedData(stim_idx1,3))
xlim([0.37 1]);
ylim([0.39 0.95])
hold on
plot(mazeDecoded(:,1),mazeDecoded(:,2),'k','LineWidth',3)
plot(shockZoneDecoded(:,1),shockZoneDecoded(:,2),'r','LineWidth',3)
hold off
set(gca, 'XTickLabel', {}, 'YTickLabel', {});
annotation('textbox', [0.42 0.65 0.1 0.1], 'String', 'Thresh = 0.015 (Simulated)', 'FontWeight', 'bold', 'FontSize', 12,...
    'EdgeColor', 'none')

% Threshold 2
TimeStim2 = ts(double(stimtimes2*1e4));
subplot(425)
scatter(Data(Restrict(Xtsd, TimeStim2, 'align', 'closest')), Data(Restrict(Ytsd, TimeStim2, 'align', 'closest')));
xlim([20 60]);
ylim([20 53]);
hold on
plot(mazePosMat(:,1),mazePosMat(:,2),'k','LineWidth',3)
plot(ShockZonePosMat(:,1),ShockZonePosMat(:,2),'r','LineWidth',3)
hold off
set(gca, 'XTickLabel', {}, 'YTickLabel', {});
subplot(426)
scatter(DecodedData(stim_idx2,2), DecodedData(stim_idx2,3))
xlim([0.37 1]);
ylim([0.39 0.95])
hold on
plot(mazeDecoded(:,1),mazeDecoded(:,2),'k','LineWidth',3)
plot(shockZoneDecoded(:,1),shockZoneDecoded(:,2),'r','LineWidth',3)
hold off
set(gca, 'XTickLabel', {}, 'YTickLabel', {});
annotation('textbox', [0.42 0.43 0.1 0.1], 'String', 'Thresh = 0.01 (Simulated)', 'FontWeight', 'bold', 'FontSize', 12,...
    'EdgeColor', 'none')

% Threshold 3
TimeStim3 = ts(double(stimtimes3*1e4));
subplot(427)
scatter(Data(Restrict(Xtsd, TimeStim3, 'align', 'closest')), Data(Restrict(Ytsd, TimeStim3, 'align', 'closest')));
xlim([20 60]);
ylim([20 53]);
hold on
plot(mazePosMat(:,1),mazePosMat(:,2),'k','LineWidth',3)
plot(ShockZonePosMat(:,1),ShockZonePosMat(:,2),'r','LineWidth',3)
hold off
set(gca, 'XTickLabel', {}, 'YTickLabel', {});
subplot(428)
scatter(DecodedData(stim_idx3,2), DecodedData(stim_idx3,3))
xlim([0.37 1]);
ylim([0.39 0.95])
hold on
plot(mazeDecoded(:,1),mazeDecoded(:,2),'k','LineWidth',3)
plot(shockZoneDecoded(:,1),shockZoneDecoded(:,2),'r','LineWidth',3)
hold off
set(gca, 'XTickLabel', {}, 'YTickLabel', {});
annotation('textbox', [0.42 0.21 0.1 0.1], 'String', 'Thresh = 0.005 (Simulated)', 'FontWeight', 'bold', 'FontSize', 12,...
    'EdgeColor', 'none')
% 
% saveas(f3, '/home/mobsrick/Dropbox/MOBS_workingON/Dima/Ongoing_results/OnlineDecoding/July2020/RealVsInferred_thresholds.fig');
% saveFigure(f3, 'RealVsInferred_thresholds', '/home/mobsrick/Dropbox/MOBS_workingON/Dima/Ongoing_results/OnlineDecoding/July2020/');


%% Real stimulation
% Real data
subplot(221)
scatter(Data(Restrict(Xtsd, TimeStim, 'align', 'closest')), Data(Restrict(Ytsd, TimeStim, 'align', 'closest')));
hold on
plot(mazePosMat(:,1),mazePosMat(:,2),'k','LineWidth',3)
plot(ShockZonePosMat(:,1),ShockZonePosMat(:,2),'r','LineWidth',3)
hold off
xlim([20 60]);
ylim([20 53]);
set(gca, 'XTickLabel', {}, 'YTickLabel', {});
title('Real position', 'FontSize', 16)
subplot(222)
scatter(StimData(:,2), StimData(:,3))
hold on
plot(mazeDecoded(:,1),mazeDecoded(:,2),'k','LineWidth',3)
plot(shockZoneDecoded(:,1),shockZoneDecoded(:,2),'r','LineWidth',3)
hold off
xlim([0.37 1]);
ylim([0.39 0.95])
set(gca, 'XTickLabel', {}, 'YTickLabel', {});
title('Online inference', 'FontSize', 16)
annotation('textbox', [0.4 0.87 0.1 0.1], 'String', 'Thresh = 0.015 (Real)', 'FontWeight', 'bold', 'FontSize', 12,...
    'EdgeColor', 'none')


% Threshold 1
subplot(223)
scatter(Data(Restrict(Xtsd, RealStimTime, 'align', 'closest')), Data(Restrict(Ytsd, RealStimTime, 'align', 'closest')));
xlim([20 60]);
ylim([20 53]);
hold on
plot(mazePosMat(:,1),mazePosMat(:,2),'k','LineWidth',3)
plot(ShockZonePosMat(:,1),ShockZonePosMat(:,2),'r','LineWidth',3)
hold off
set(gca, 'XTickLabel', {}, 'YTickLabel', {});
subplot(224)
scatter(StimData(:,2), StimData(:,3))
xlim([0.37 1]);
ylim([0.39 0.95])
hold on
plot(mazeDecoded(:,1),mazeDecoded(:,2),'k','LineWidth',3)
plot(shockZoneDecoded(:,1),shockZoneDecoded(:,2),'r','LineWidth',3)
hold off
set(gca, 'XTickLabel', {}, 'YTickLabel', {});
annotation('textbox', [0.42 0.65 0.1 0.1], 'String', 'Thresh = 0.015 (Simulated)', 'FontWeight', 'bold', 'FontSize', 12,...
    'EdgeColor', 'none')


