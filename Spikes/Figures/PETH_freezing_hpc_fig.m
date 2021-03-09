clear all

%% Parameters
IsSave=0;
% mice = [797 798 828 861 882 905 906 911 912 977 994 1117 1124 1161 1162 1168];
mice = [797 798 828 861 882 905 906 911 912];
% mice = [994 1117 1124 1161 1162 1168];
% mice = 994;
wi = 2; % in sec
binsize = 0.05; % in sec
nbins = 81;

%% Get data
Dir = PathForExperimentsERC_Dima('UMazePAG');
Dir = RestrictPathForExperiment(Dir,'nMice', mice);

% Allocate arrays
b = cell(length(Dir.path),1); % behavior
s = cell(length(Dir.path),1); % spikes
CondEpoch = cell(length(Dir.path),1);
FreezeEpoch = cell(length(Dir.path),1);

for imouse = 1:length(Dir.path)
    b{imouse} = load([Dir.path{imouse}{1} '/behavResources.mat'], 'SessionEpoch', 'FreezeAccEpoch', 'MovAcctsd');
    s{imouse} = load([Dir.path{imouse}{1} '/SpikeData.mat'], 'S', 'BasicNeuronInfo', 'PlaceCells');
end

%% Organize epochs and data
% Epochs
for imouse = 1:length(Dir.path)
    [~, ~, CondEpoch{imouse}] = ReturnMnemozyneEpochs(b{imouse}.SessionEpoch);
    
    FreezeEpoch{imouse} = and(b{imouse}.FreezeAccEpoch,CondEpoch{imouse});
    FreezeEpoch{imouse} = dropShortIntervals(FreezeEpoch{imouse}, 4*1E4);
    FreezeEpoch{imouse} = mergeCloseIntervals(FreezeEpoch{imouse},1*1E4);
end
% Data
for imouse = 1:length(Dir.path)
    Q = MakeQfromS(s{imouse}.S, binsize*1e4);
    time = Restrict(Q, CondEpoch{imouse});
    tempdat = zscore(full(Data(Restrict(Q, CondEpoch{imouse}))));
    for icell = 1:length(s{imouse}.S)
        Neurons{imouse}{icell} = [Range(time)/1e4 tempdat(:,icell)];
    end
end
    
%% Calculate PETH
cnt=1;
for imouse = 1:length(Dir.path)
    for icell = 1:length(Neurons{imouse})
        [r,i] = Sync(Neurons{imouse}{icell}, Start(FreezeEpoch{imouse})/1e4 ,'durations',[-wi wi]);
        T = SyncMap(r,i,'durations',[-wi wi],'nbins',nbins,'smooth',0);
        OnResp(cnt,:) = mean(T);
        
        [r,i] = Sync(Neurons{imouse}{icell}, End(FreezeEpoch{imouse})/1e4 ,'durations',[-wi wi]);
        T = SyncMap(r,i,'durations',[-wi wi],'nbins',nbins,'smooth',0);
        OffResp(cnt,:) = mean(T);
        
        cnt = cnt + 1;
    end
end

%% Calculate accelerometer trace
Accelero = nan(length(Dir.path), 201*2);
for imouse = 1:length(Dir.path)
    M_st = PlotRipRaw(b{imouse}.MovAcctsd,Start(FreezeEpoch{imouse},'s'),2000,1,0,0);
    M_end = PlotRipRaw(b{imouse}.MovAcctsd,End(FreezeEpoch{imouse},'s'),2000,1,0,0);
    
    time = [(M_st(:,1))' (M_end(:,1)+4)'];
    Accelero(imouse,:) = [(M_st(:,2))' (M_end(:,2))'];
end
    

%% Plot figure
f1 = figure('units', 'normalized', 'outerposition', [0.1 0.3 0.4 0.8]);
DatNormZ = [OnResp,OffResp];
SustVal = nanmean(DatNormZ(:,60:100),2);
UseForTresh = SustVal;
[a,ind] = sort(UseForTresh);
subplot(4,1,1:3)
imagesc(-2:0.05:6,1:size(DatNormZ,1),DatNormZ(ind,:))
caxis([-0.3 0.3])
line([0 0],ylim,'color','k')
line([4 4],ylim,'color','k')
makepretty
ylabel('Neuron #')
xlabel('Time to freeze on (s)')
% Accelero meter
subplot(4,1,4)
if size(Accelero,1) > 1
    shadedErrorBar(time,mean(Accelero, 1), std(Accelero, 1), 'k')
else
    plot(time,Accelero,'color','k')
end
makepretty
ylabel('Accelero')
xlabel('Time to freeze on (s)')

% Save figure
if IsSave
    foldertosave = ChooseFolderForFigures_DB('Spikes');
    saveas(f1,[foldertosave '/PETHonFreezing_hpc.fig']);
    saveFigure(f1, 'PETHonFreezing_hpc', foldertosave);
end