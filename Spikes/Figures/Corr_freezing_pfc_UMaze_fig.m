clear all
%
% Correlation matrix for neuronal activity during freeezing (PAG EmbReact)
%

% Coded by Dima Bryzgalov and Samuel Laventure, MOBS team, Paris, France
% 03/2021
% github.com/bryzgalovdm
%% Parameters
IsSave=0;
mice = [490 507 508 509 510 512 514];
% mice = [994];
wi = 2; % in sec
binsize = 0.1; % in sec
nbins = 81;
timeatTransition = 3; % in sec
tps=[0.05:0.05:1];


%% Get data
% Cond
Dir_temp = PathForExperimentsEmbReact('UMazeCond');
n=1;
for imouse = 1:length(Dir_temp.ExpeInfo)
    if ismember(Dir_temp.ExpeInfo{imouse}{1}.nmouse,mice)
        Dir_Cond.path{n} = Dir_temp.path{imouse};
        Dir_Cond.ExpeInfo{n} = Dir_temp.ExpeInfo{imouse};
        n=n+1;
    end
end

% Allocate arrays
s_cond = cell(length(Dir_Cond.path), 1); % Spikes
r_cond = cell(length(Dir_Cond.path), 1); %
SessEpoch = cell(length(Dir_Cond.path), 1);
f_cond = cell(length(Dir_Cond.path), 1);
pfc_id = cell(length(Dir_Cond.path), 1);
neuroclass = cell(length(Dir_Cond.path), 1);
clu_pfc = cell(length(Dir_Cond.path), 1);
Neurons = cell(length(Dir_Cond.path),1);
FR = nan(5e4,1);
NeuronClass = nan(5e4,1);

for imouse = 1:length(Dir_Cond.path)
    s_cond{imouse} = ConcatenateDataFromFolders_SB(Dir_Cond.path{imouse},'Spikes');
    f_cond{imouse} = ConcatenateDataFromFolders_SB(Dir_Cond.path{imouse},'Epoch', 'Epochname', 'freezeepoch');
    SessEpoch{imouse} = ConcatenateDataFromFolders_SB(Dir_Cond.path{imouse},'Epoch','epochname','sessiontype');
    StimEpoch{imouse} = ConcatenateDataFromFolders_SB(Dir_Cond.path{imouse},'Epoch','epochname','stimepoch');
    accelero{imouse} = ConcatenateDataFromFolders_SB(Dir_Cond.path{imouse},'accelero');
    pfc_id{imouse} = load([Dir_Cond.path{imouse}{1} 'SpikesToAnalyse/PFCx_Neurons.mat']);
    neuroclass{imouse} = load([Dir_Cond.path{imouse}{1} 'WFIdentit.mat']);
end

%% Organize epochs and data
% Epochs
for imouse = 1:length(Dir_Cond.path)
    % Preprocess freezing: remove short epochs and badly separated epochs
    FreezeEpoch{imouse} = f_cond{imouse};
    FreezeEpoch{imouse} = dropShortIntervals(FreezeEpoch{imouse}, 4*1E4);
    FreezeEpoch{imouse} = mergeCloseIntervals(FreezeEpoch{imouse},1*1E4); 
    % Remove epochs co-terminated with stims
    FreezeEpoch{imouse} = GetStandaloneFreezeEpochs(FreezeEpoch{imouse}, StimEpoch{imouse}, 2.5);
    
    for fr_ep=1:length(Start(FreezeEpoch{imouse}))
        LitEp=subset(FreezeEpoch{imouse},fr_ep);
        NotLitEp=FreezeEpoch{imouse}-LitEp;
        StopEp=intervalSet(Stop(LitEp)-timeatTransition*1e4,Stop(LitEp)+timeatTransition*1e4);
        if not(isempty(Data(Restrict(accelero{imouse},and(NotLitEp,StopEp)))))
            FreezeEpochBis=FreezeEpoch{imouse}-LitEp;
        end
    end
    FreezeEpoch{imouse}=CleanUpEpoch(FreezeEpochBis);
    
end

% Data
for imouse = 1:length(Dir_Cond.path)
    clu_pfc{imouse} = s_cond{imouse}(pfc_id{imouse}.number);
    NClass{imouse} = neuroclass{imouse}.UnitID(pfc_id{imouse}.number,1);
    
    Neurons{imouse} = cell(length(clu_pfc),1);
    Q = MakeQfromS(clu_pfc{imouse}, binsize*1e4);
    time = Range(Q);
    tempdat = zscore(full(Data(Q)));
    Neurons{imouse} = tsd(time, tempdat);
end

%% Normalize firing histograms in time
cnt=1;
for imouse = 1:length(Dir_Cond.path)
    if ~(isempty(Neurons{imouse}))
        for icell = 1:length(clu_pfc{imouse})
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
            NormSpikes(cnt,:) = nanmean(NormEpochs);
            cnt = cnt + 1;
        end
    else
        NormSpikes = [];
    end
end

%% Plot correlation matrix
Corrs = corr(NormSpikes);

f1 = figure('units', 'normalized', 'outerposition', [0.1 0.3 0.35 0.6]);
imagesc(Corrs')
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
    saveas(f1,[foldertosave '/Freezing/CorrFreezing_pfc_UMaze.fig']);
    saveFigure(f1, 'CorrFreezing_pfc_UMaze', [foldertosave '/Freezing']);
end