clear all

%% Parameters
IsSave=0;
mice = [797 798 828 861 882 905 906 911 912 977 994 1117 1124 1161 1162 1168];
% mice = [994];
wi = 2; % in sec
binsize = 0.05; % in sec
nbins = 81;
timeatTransition = 3;
tps=[0.05:0.05:1];

%% Get data
Dir = PathForExperimentsERC_Dima('UMazePAG');
Dir = RestrictPathForExperiment(Dir,'nMice', mice);

% Allocate arrays
b = cell(length(Dir.path),1); % behavior
s = cell(length(Dir.path),1); % spikes
CondEpoch = cell(length(Dir.path),1);
FreezeEpoch = cell(length(Dir.path),1);
Neurons = cell(length(Dir.path),1);

for imouse = 1:length(Dir.path)
    b{imouse} = load([Dir.path{imouse}{1} '/behavResources.mat'], 'SessionEpoch', 'FreezeAccEpoch', 'MovAcctsd');
    s{imouse} = load([Dir.path{imouse}{1} '/SpikeData.mat'], 'S', 'BasicNeuronInfo', 'PlaceCells');
end

%% Organize epochs and data
% Epochs
for imouse = 1:length(Dir.path)
    [~, ~, CondEpoch{imouse}] = ReturnMnemozyneEpochs(b{imouse}.SessionEpoch);
    
    
    % Preprocess freezing: remove short epochs and badly separated epochs
    FreezeEpoch{imouse} = and(b{imouse}.FreezeAccEpoch,CondEpoch{imouse});
    FreezeEpoch{imouse} = dropShortIntervals(FreezeEpoch{imouse}, 4*1E4);
    FreezeEpoch{imouse} = mergeCloseIntervals(FreezeEpoch{imouse},1*1E4);    
    for fr_ep=1:length(Start(FreezeEpoch{imouse}))
        LitEp=subset(FreezeEpoch{imouse},fr_ep);
        NotLitEp=FreezeEpoch{imouse}-LitEp;
        StopEp=intervalSet(Stop(LitEp)-timeatTransition*1e4,Stop(LitEp)+timeatTransition*1e4);
        if not(isempty(Data(Restrict(b{imouse}.MovAcctsd,and(NotLitEp,StopEp)))))
            FreezeEpochBis=FreezeEpoch{imouse}-LitEp;
        end
    end
    FreezeEpoch{imouse}=CleanUpEpoch(FreezeEpochBis);
    
end

% Data
for imouse = 1:length(Dir.path)
    Q = MakeQfromS(s{imouse}.S, binsize*1e4);
    time = Range(Restrict(Q, CondEpoch{imouse}));
    tempdat = zscore(full(Data(Restrict(Q, CondEpoch{imouse}))));
%     tempdat = full(Data(Restrict(Q, CondEpoch{imouse})));
    Neurons{imouse} = tsd(time, tempdat);
end

%% Normalize firing histograms in time
cnt=1;
for imouse = 1:length(Dir.path)
    if ~(isempty(Neurons{imouse}))
        for icell = 1:length(s{imouse}.S)
            cnt_epoch = 1;
            for iepoch=1:length(Start(FreezeEpoch{imouse}))-1
                
                ActualEpoch=subset(FreezeEpoch{imouse},iepoch);
                LittleEpoch=intervalSet(Start(ActualEpoch),Stop(ActualEpoch));
                LittleEpochPre=intervalSet(Start(ActualEpoch)-timeatTransition*1e4,Start(ActualEpoch));
                LittleEpochPost=intervalSet(Stop(ActualEpoch),Stop(ActualEpoch)+timeatTransition*1e4);
                
                TempData=Data(Restrict(Neurons{imouse},LittleEpoch));
                TempDataPre=Data(Restrict(Neurons{imouse},LittleEpochPre));
                TempDataPost=Data(Restrict(Neurons{imouse},LittleEpochPost));
                
                if ~isempty(TempData) && ~isempty(TempDataPre) && ~isempty(TempDataPost)
                    neuron = TempData(:, icell);
                    neuron = interp1([1/length(neuron):1/length(neuron):1],neuron,tps);
                    
                    neuron_pre = TempDataPre(:, icell);
                    neuron_pre = interp1([1/length(neuron_pre):1/length(neuron_pre):1],neuron_pre,[0.1:0.1:1]);
                    
                    neuron_post = TempDataPost(:,icell);
                    neuron_post = interp1([1/length(neuron_post):1/length(neuron_post):1],neuron_post,[0.1:0.1:1]);
                    
                    NormEpochs(cnt_epoch,:) = [neuron_pre neuron neuron_post];
                    cnt_epoch = cnt_epoch + 1;
                end
            end
            NormSpikes(cnt,:) = mean(NormEpochs);
            cnt = cnt + 1;
        end
    else
        NormSpikes = [];
    end
end

%% Plot correlation matrix
Corrs = corr(NormSpikes);

f1 = figure('units', 'normalized', 'outerposition', [0.1 0.3 0.35 0.6]);
imagesc(Corrs)
axis xy
colormap(fliplr(brewermap(200,'Spectral')')')
caxis([-0.5 0.5])
set(gca, 'XTick', [10 30], 'XTickLabel', {'0', '1'});
set(gca, 'YTick', [10 30], 'YTickLabel', {'0', '1'});
ylabel('Normalized freezing time')
xlabel('Normalized freezing time')
title('Neuronal correlations during freezing time')
colorbar
makepretty

% Save figure
if IsSave
    foldertosave = ChooseFolderForFigures_DB('Spikes');
    saveas(f1,[foldertosave '/CorrFreezing_hpc.fig']);
    saveFigure(f1, 'CorrFreezing_hpc', foldertosave);
end