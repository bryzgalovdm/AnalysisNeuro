function fh = CheckThresholdDetection(datadir, ChanList, ChannelThresholded, ThresholdUsed, filterfreq, savedir)

% This funtion performs basic check-ups of quality of INTAN-based spike
% thresholding. Each detectoin should be followed by TTL that indicate
% stimulation
% You will need to provide the folder with Open-Ephys-recorded events (.mat)
% and ExpeInfo with correct info about digital inputs
% 
%  INPUT
% 
%     folder          folder with Open-Ephys-recorded events (.mat)
%     ExpeInfo        ExpeInfo with correct info about digital inputs
%
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 10/09/2020 based on earlier function written by me from 2019


%%%%% ------- M994 Ch46 ThresholdUsed -------- %%%%%%%%%%

%% Allocate memory
channels = cell(length(ChanList), 1);
channelsF = cell(length(ChanList), 1);

%% Load data
for ichan = 1:length(ChanList)
    if exist([datadir '/UnfilteredChans/Chan' num2str(ChanList(ichan)) '.mat'], 'file') == 2
        load([datadir '/UnfilteredChans/Chan' num2str(ChanList(ichan)) '.mat']);
        channels{ichan} = raw;
    else
        MakeData_Detection(datadir, ChanList(ichan));
        load([datadir '/UnfilteredChans/Chan' num2str(ChanList(ichan)) '.mat']);
        channels{ichan} = raw;
    end
end
load([datadir '/behavResources.mat'],'StimEpoch');

%% Filter channels
for ichan = 1:length(ChanList)
    chanF = INTANfilter((Data(channels{ichan}))',0,filterfreq);
    channelsF{ichan} = tsd(Range(channels{ichan}), chanF');
end

%% Get Latencies
ISTemp = intervalSet(Start(StimEpoch)-300,Start(StimEpoch));
NotFound = 0;
L = nan(1,length(Start(StimEpoch)));
NotFound_idx = zeros(1,length(Start(StimEpoch)));
NegPeak = zeros(1,length(Start(StimEpoch)));

for i=1:length(Start(StimEpoch))
    PerTemp = Restrict(channelsF{ChannelThresholded},subset(ISTemp,i));
    
    Crossed{i} = thresholdIntervals(PerTemp,ThresholdUsed/0.195+1,'Direction','Below');
    if length(Start(Crossed{i})) > 1
        for j=1:length(Start(Crossed{i}))
            try
                PPP(j) = min(Data(Restrict(PerTemp,subset(Crossed{i},j))));
            catch
                datatemp = Data(PerTemp); % sometimes the very last point crosses threshold - than Restrict doesn't work
                timetemp = Range(PerTemp);
                PPP(j) = datatemp(timetemp == Start(subset(Crossed{i},j)));
                clear datatemp timetemp
            end
            idxL(j) = find(Data(PerTemp)==PPP(j));
        end
        [idxR,ii] = max(idxL);
        NegPeak(i) = PPP(ii);
        t{i} = Range(PerTemp);
        L(i) = (End(subset(ISTemp,i)) - t{i}(idxR))/10;
    elseif length(Start(Crossed{i})) == 1
        try
            NegPeak(i) = min(Data(Restrict(PerTemp,Crossed{i})));
        catch
            datatemp = Data(PerTemp); % sometimes the very last point crosses threshold - than Restrict doesn't work
            timetemp = Range(PerTemp);
            NegPeak(i) = datatemp(timetemp == Start(Crossed{i}));
            clear datatemp timetemp
        end
        idx = find(Data(PerTemp)==NegPeak(i));
        t{i} = Range(PerTemp);
        L(i) = (End(subset(ISTemp,i)) - t{i}(idx))/10;
    else
        NotFound = NotFound+1;
        NotFound_idx(i) = i;
    end
end
NotFound_idx = nonzeros(NotFound_idx);

% clear IsTemp NotFound PerTemp NegPeak t L PP
clear IsTemp PerTemp t PP

%% Get ERPs
% % RECODE
% [M59,T59] = PlotRipRaw(Chan59F,Start(StimEpoch)/1e4,[-4 4]);
% [M45,T45] = PlotRipRaw(Chan45F,Start(StimEpoch)/1e4,[-4 4]);
% [M46,T46] = PlotRipRaw(Chan46F,Start(StimEpoch)/1e4,[-4 4]);
% [M47,T47] = PlotRipRaw(Chan47F,Start(StimEpoch)/1e4,[-4 4]);
% close all
% 
% 
% % Reallign spikes
% Chan59FAligned = zeros(length(L),32);
% Chan45FAligned = zeros(length(L),32);
% Chan46FAligned = zeros(length(L),32);
% Chan47FAligned = zeros(length(L),32);
% for isig=1:4 %for each line
%     for k = 1:length(L)
%         [z,id] = min(T46(k,:));
%         if ~isempty(id)
%             Chan59FAligned(k,:) = T59(k,id-13:id+18);
%             Chan45FAligned(k,:) = T45(k,id-13:id+18);
%             Chan46FAligned(k,:) = T46(k,id-13:id+18);
%             Chan47FAligned(k,:) = T47(k,id-13:id+18);
%         end
%     end
% end
% CH_Online = {Chan59FAligned; Chan45FAligned; Chan46FAligned; Chan47FAligned};

%% Compare offline and online
% Two units 
allfiles = dir([pwd '/Waveforms/**/*.mat']);

for i=1:length(allfiles)
    %     if isempty(strfind(allfiles(i).name,'c1.mat'))
    W{i} = load([allfiles(i).folder '/' allfiles(i).name]);
    W{i} = W{i}.W;
    %     end
end
    
% Find all the events thresspassing the threshold
for i=1:length(W)
    if ~isempty(W{i})
        idx_W{i} = zeros(1,size(W{i},1));
        for k=1:size(W{i},1)
            if sum(squeeze(W{i}(k,3,:))<(ThresholdUsed/0.195))>0
                idx_W{i}(k) = k;
            end
        end
        idx_W{i} = nonzeros(idx_W{i});
    end
end

% Number of threshold exceeding spikes
NSpikesClustered = 0;
for i=1:length(idx_W)
    if ~isempty(idx_W{i})
        NSpikesClustered = NSpikesClustered + length(idx_W{i});
    end
end
MissedPerc = (NSpikesClustered - size(T46,1))/NSpikesClustered*100;
DetectedPerc = 100 - MissedPerc;


% Pool all tresspassing spikes
AllSpikesClustered = zeros(NSpikesClustered,4,32);
increment = 0;
for i=1:length(idx_W)
    if ~isempty(idx_W{i})
        for k = 1:length(idx_W{i})
            AllSpikesClustered(increment+k,:,:) = W{i}(idx_W{i}(k),:,:);
        end
        increment = increment+length(idx_W{i});
    end
end

% Percentage of different clusters
Perc_Cl = zeros(1,length(W));
for i = 1:length(W)
    Perc_Cl(i) = round(length(idx_W{i})/NSpikesClustered*100,1);
end
Perc_Cl = nonzeros(Perc_Cl);

%% Plot histogram of delays
% Check how many detections are > 2 ms
if length(find(L>2)) > 1
    warning('There is more than 1 huge outliers from expected quality');
else
    L_toplot = L;
    L_toplot(L>2) = [];
end

fh = figure('rend','painters','pos',[1 1 650 400]);
% fh = figure('units', 'normalized', 'outerposition', [0 0 0.7 0.6]);

hist(L_toplot,20);
h = findobj(gca,'Type','patch');
h.FaceColor = [0 0 0];
xlim([0 1.5])
set(gca,'LineWidth',5,'FontSize',18,'FontWeight','bold');
title('Detection-Stimulation delays')
xlabel('Time (ms)')
ylabel('Count')
box off


saveas(fh,[savedir '/DelayHist_Striatum.fig']);
saveFigure(fh,'DelayHist_Striatum',savedir);
    
%% Plot ERP

fh = figure('rend','painters','pos',[1 1 800 600]);
shadedErrorBar(M(:,1)*1000, mean(T,1),std(T,1),{'k','LineWidth',3});
set(gca,'LineWidth',5,'FontSize',18,'FontWeight','bold');
title('Filtered data triggered on the stimulation')
xlabel('Time (ms)')
ylabel('Volage (mkv)')
box off

saveas(fh,['/home/mobsrick/Dropbox/MOBS_workingON/Dima/Ongoing results/Spikes/Thibault_test/October2019/ERP.fig']);
    saveFigure(fh,'ERP',...
        '/home/mobsrick/Dropbox/MOBS_workingON/Dima/Ongoing results/Spikes/Thibault_test/October2019/');
    
%% Plot the spikes offline
fi = figure('rend','painters','pos',[1 1 650 1200]);
for i = 1:size(AllSpikesClustered,2)
    subplot(4,1,i)
    add1 = plot(squeeze(AllSpikesClustered(1:20:end,i,:))','color',[0.6 0.6 0.6]);
    hold on
    main1 = plot(nanmean(squeeze(AllSpikesClustered(:,i,:))),'color','k','linewidth',2);
    ylim([-5000 1000])
    xlim([0 32])
    hold on
    set(gca,'FontWeight','bold','FontSize',16,'LineWidth',3);
    box off
    if i == 1
       title(['Spikes detected offline'],'FontSize',18);
       set(gca, 'YTickLabel',[],'XTickLabel',[]);
    elseif i == 3
        line(xlim,[ThresholdUsed/0.195 ThresholdUsed/0.195],'Color','r','LineWidth',3);
        set(gca, 'YTickLabel',[],'XTickLabel',[]);
    elseif i == 4
        xlabel('Time in samples')
    else
        set(gca, 'YTickLabel',[],'XTickLabel',[]);
    end
end

saveas(fi,['/home/mobsrick/Dropbox/MOBS_workingON/Dima/Ongoing results/Spikes/Thibault_test/October2019/SpikesOffline.fig']);
    saveFigure(fi,'SpikesOffline',...
        '/home/mobsrick/Dropbox/MOBS_workingON/Dima/Ongoing results/Spikes/Thibault_test/October2019/');
    
    
    
%% Plot the spikes online
fi = figure('rend','painters','pos',[1 1 650 1200]);
for i = 1:length(CH_Online)
    subplot(4,1,i)
    add1 = plot(CH_Online{i}(1:20:end,:)','color',[0.6 0.6 0.6]);
    hold on
    main1 = plot(nanmean(CH_Online{i}),'color','k','linewidth',2);
    ylim([-5000 1000])
    xlim([0 32])
    hold on
    set(gca,'FontWeight','bold','FontSize',16,'LineWidth',3);
    box off
    if i == 1
       title(['Spikes detected online'],'FontSize',18);
       set(gca, 'YTickLabel',[],'XTickLabel',[]);
    elseif i == 3
        line(xlim,[ThresholdUsed/0.195 ThresholdUsed/0.195],'Color','r','LineWidth',3);
        set(gca, 'YTickLabel',[],'XTickLabel',[]);
    elseif i == 4
        xlabel('Time in samples')
    else
        set(gca, 'YTickLabel',[],'XTickLabel',[]);
    end
end

saveas(fi,['/home/mobsrick/Dropbox/MOBS_workingON/Dima/Ongoing results/Spikes/Thibault_test/October2019/SpikesOnline.fig']);
saveFigure(fi,'SpikesOnline',...
    '/home/mobsrick/Dropbox/MOBS_workingON/Dima/Ongoing results/Spikes/Thibault_test/October2019/');

%% Plot the spikes offline and online
fi = figure('rend','painters','pos',[1 1 920 1200]);
for i = 1:size(AllSpikesClustered,2)
    subplot(4,2,i*2-1)
    add1 = plot(squeeze(AllSpikesClustered(1:20:end,i,:))','color',[0.6 0.6 0.6]);
    hold on
    main1 = plot(nanmean(squeeze(AllSpikesClustered(:,i,:))),'color','k','linewidth',2);
    ylim([-5000 1000])
    xlim([0 32])
    hold on
    set(gca,'FontWeight','bold','FontSize',16,'LineWidth',3);
    box off
    if i == 1
        title(['Spikes detected offline'],'FontSize',18);
        set(gca, 'YTickLabel',[],'XTickLabel',[]);
    elseif i == 3
        line(xlim,[ThresholdUsed/0.195 ThresholdUsed/0.195],'Color','r','LineWidth',3);
        set(gca, 'YTickLabel',[],'XTickLabel',[]);
    elseif i == 4
        xlabel('Time in samples')
    else
        set(gca, 'YTickLabel',[],'XTickLabel',[]);
    end
end
for i = 1:length(CH_Online)
    subplot(4,2,2*i)
    add1 = plot(CH_Online{i}(1:20:end,:)','color',[0.6 0.6 0.6]);
    hold on
    main1 = plot(nanmean(CH_Online{i}),'color','k','linewidth',2);
    ylim([-5000 1000])
    xlim([0 32])
    hold on
    set(gca,'FontWeight','bold','FontSize',16,'LineWidth',3);
    box off
    if i == 1
       title(['Spikes detected online'],'FontSize',18);
       set(gca, 'YTickLabel',[],'XTickLabel',[]);
    elseif i == 3
        line(xlim,[ThresholdUsed/0.195 ThresholdUsed/0.195],'Color','r','LineWidth',3);
        set(gca, 'YTickLabel',[],'XTickLabel',[]);
    elseif i == 4
        xlabel('Time in samples')
    else
        set(gca, 'YTickLabel',[],'XTickLabel',[]);
    end
end

saveas(fi,['/home/mobsrick/Dropbox/MOBS_workingON/Dima/Ongoing results/Spikes/Thibault_test/October2019/SpikesOfflineOnline.fig']);
saveFigure(fi,'SpikesOfflineOnline',...
    '/home/mobsrick/Dropbox/MOBS_workingON/Dima/Ongoing results/Spikes/Thibault_test/October2019/');


%% Bar detected thresh
fi = figure('rend','painters','pos',[1 1 300 600]);
a = [MissedPerc DetectedPerc; nan(1,2)];
h = bar(a,'stacked');
set(h,{'FaceColor'},{[0.9 0.1 0.1];[0.1 0.1 0.9]});
set(gca,'XTick',[], 'XtickLabel', {},'LineWidth',5,'FontSize',18,'FontWeight','bold');
ylim([0 110])
xlim([0 2])
ylabel('Percentage')
box off
% legend('Missed spikes','Detected spikes')
% legend('boxoff')

saveas(fi,['/home/mobsrick/Dropbox/MOBS_workingON/Dima/Ongoing results/Spikes/Thibault_test/October2019/MissDet.fig']);
    saveFigure(fi,'MissDet',...
        '/home/mobsrick/Dropbox/MOBS_workingON/Dima/Ongoing results/Spikes/Thibault_test/October2019/');


%% Pie different spikes
fi = figure('rend','painters','pos',[1 1 500 500]);
pieid  = pie(Perc_Cl);
colormap jet
for i=2:2:length(pieid)
    pieid(i).FontSize = 14;
    pieid(i).FontWeight = 'bold';
end
legend ({'MUA', 'Cl11','Cl13','Cl6','Cl8'}, 'Position', [0.8 0.85 0.17 0.035],'FontSize',14);
legend('boxoff')

saveas(fi,['/home/mobsrick/Dropbox/MOBS_workingON/Dima/Ongoing results/Spikes/Thibault_test/October2019/Pie_Clusters.fig']);
    saveFigure(fi,'Pie_Clusters',...
        '/home/mobsrick/Dropbox/MOBS_workingON/Dima/Ongoing results/Spikes/Thibault_test/October2019/');