function fh = FigureRealVSInferred(folder, thresh)
% 
% This funciton will plot the figure with two sunplots:
% - Real position in test set
% - Inferred position in the test set
% Inferred position are filtered acording to the evaluated loss (confidence of the network)
% 
%  INPUT
%  
%         folder          folder with inferring.mat
%         thresh          threshold on evaluated loss (data LESS than <thresh> are taken)
%         
%  OUTPUT
%  
%         fh              figure handle
% 
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 23/06/2020
% github.com/bryzgalovdm
        
%% Load data
load([folder '/inferring.mat']);

%% Figure
fh = figure('units', 'normalized', 'outerposition', [0.1 0.2 0.4 0.85]);

%% Real position
% Plot straight away
subplot(211)
scatter(pos(:,1), pos(:,2), 14, times)
% plot(pos(:,1), pos(:,2))
xl = xlim;
yl = ylim;
title('Real position')
%% Inferred position
% Remove all [0,0] positions (noise)
inferringPre = inferring(~inferring(:,1)==0,:);
timesInf = times(~inferring(:,1)==0);
% Remove all the points with evaluated loss lower than 0.01 (keep best ~20% of data)
inferringPre = inferringPre(inferringPre(:,3)<thresh,:);
timesInf = timesInf(inferringPre(:,3)<thresh);
% Smooth the data (factor=10)
X = runmean(inferringPre(:,1), 10);
Y = runmean(inferringPre(:,2), 10);
inferringSmooth_lowErr = [X Y inferringPre(:,3)];

% Plot the figure
subplot(212)
% plot(inferringSmooth_lowErr(:,1), inferringSmooth_lowErr(:,2))
scatter(inferringSmooth_lowErr(:,1), inferringSmooth_lowErr(:,2), 14, timesInf)
xlim(xl)
ylim(yl)
title('Inferred position')

%% Save the figure
saveas(fh,'RealvsInf.fig');
saveFigure(fh,'RealvsInf',pwd);
