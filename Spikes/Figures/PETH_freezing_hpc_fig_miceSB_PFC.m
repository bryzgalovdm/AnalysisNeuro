function f1 = PETH_freezing_hpc_fig_miceSB_PFC(varargin)
% 
% This function plots PSTH for the [-2 2] s period of freezing onset and
% [-2 2] s period of freezing offset
% 
% INPUT
%   
%   IsSaveFig (optional)            Boolean: save all figures or not (defalut = false)
%   IsSaveFile (optional)           Boolean: save file with data or data (defalut = false)
%   IsCalcAdditional (optional)     Boolean: record additional parameters for clusters or not (defalut = false)
%                                   Parameters are: spatialInfo, PlaceCells
%                                   or not, firing rate, tentative neuron
%                                   class
%   Redo (optional)                 Boolean: if true recalcalculate
%                                   everything from scratch (default=true)
%   Mice (optional)                 Number of mice in the analysis (default - all PAG)
% 
% OUTPUT
% 
%   f1                              figure with PETH
% 
% EXAMPLE
% 
%   f1 = PETH_freezing_hpc_fig('IsSaveFile', 1, 'IsCalcAdditional', 1, 'ReDo', 1);
% 
% 
% By Dima Bryzgalov, MOBS team, Paris,
% ~/03/2021
% github.com/bryzgalovdm

%% Parameters
IsSaveFig=0;
IsSaveFile = 0;
Redo = 1;
IsCalcAdditional = 0;

mice = [490 507 508 509 510 512 514];

wi = 2; % in sec
binsize = 0.5; % in sec
nbins = 81;
speed_thresh = 5; % in cm/s

%% Optional parameters handling
for i=1:2:length(varargin)    
    switch(lower(varargin{i})) 
        case 'issavefig'
            IsSaveFig = varargin{i+1};
            if IsSaveFig ~= 0 && IsSaveFig ~= 1 
                error('Incorrect value for property ''IsSaveFig'' (type ''help PETH_freezing_hpc_fig'' for details).');
            end
        case 'issavefile'
            IsSaveFile = varargin{i+1};
            if IsSaveFile ~= 0 && IsSaveFile ~= 1 
                error('Incorrect value for property ''IsSaveFile'' (type ''help PETH_freezing_hpc_fig'' for details).');
            end
        case 'iscalcadditional'
            IsCalcAdditional = varargin{i+1};
            if IsCalcAdditional ~= 0 && IsCalcAdditional ~= 1 
                error('Incorrect value for property ''IsCalcAdditional'' (type ''help PETH_freezing_hpc_fig'' for details).');
            end
        case 'redo'
            Redo = varargin{i+1};
            if Redo ~= 0 && Redo ~= 1 
                error('Incorrect value for property ''Redo'' (type ''help PETH_freezing_hpc_fig'' for details).');
            end
         case 'mice'
            mice = varargin{i+1};
            if ~isa(mice, 'numeric')
                error('Incorrect value for property ''mice'' (type ''help PETH_freezing_hpc_fig'' for details).');
            end
    end
end

if Redo
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
    accelero = cell(length(Dir_Cond.path), 1);
    neuroclass = cell(length(Dir_Cond.path), 1);
    clu_pfc = cell(length(Dir_Cond.path), 1);
    Neurons = cell(length(Dir_Cond.path),1);
    FR = nan(5e4,1);
    NeuronClass = nan(5e4,1);

    for imouse = 1:length(Dir_Cond.path)
        s_cond{imouse} = ConcatenateDataFromFolders_SB(Dir_Cond.path{imouse},'Spikes');
        r_cond{imouse} = ConcatenateDataFromFolders_SB(Dir_Cond.path{imouse},'Ripples');
        f_cond{imouse} = ConcatenateDataFromFolders_SB(Dir_Cond.path{imouse},'Epoch', 'Epochname', 'freezeepoch');
        SessEpoch{imouse} = ConcatenateDataFromFolders_SB(Dir_Cond.path{imouse},'Epoch','epochname','sessiontype');
        StimEpoch{imouse} = ConcatenateDataFromFolders_SB(Dir_Cond.path{imouse},'Epoch','epochname','stimepoch');
        accelero{imouse} = ConcatenateDataFromFolders_SB(Dir_Cond.path{imouse},'accelero');
        pfc_id{imouse} = load([Dir_Cond.path{imouse}{1} 'SpikesToAnalyse/PFCx_Neurons.mat']);
        neuroclass{imouse} = load([Dir_Cond.path{imouse}{1} 'WFIdentit.mat']);
    end
    
    
    
    
    %% Organize epochs and data
    % Epoch
    for imouse = 1:length(Dir_Cond.path)
        CondEpoch{imouse} = or(SessEpoch{imouse}.UMazeCond{1}, SessEpoch{imouse}.UMazeCond{2});
        CondEpoch{imouse} = or(CondEpoch{imouse}, SessEpoch{imouse}.UMazeCond{3});
        CondEpoch{imouse} = or(CondEpoch{imouse}, SessEpoch{imouse}.UMazeCond{4});
        CondEpoch{imouse} = or(CondEpoch{imouse}, SessEpoch{imouse}.UMazeCond{5});
        
        % Freezeepoch
        FreezeEpoch{imouse} = f_cond{imouse};
        FreezeEpoch{imouse} = dropShortIntervals(FreezeEpoch{imouse}, 4*1E4);
        FreezeEpoch{imouse} = mergeCloseIntervals(FreezeEpoch{imouse},1*1E4);
        % Remove epochs co-terminated with stims
        FreezeEpoch{imouse} = GetStandaloneFreezeEpochs(FreezeEpoch{imouse}, StimEpoch{imouse}, 2.5);
        
    end
    
    % Data
    for imouse = 1:length(Dir_Cond.path)
        clu_pfc{imouse} = s_cond{imouse}(pfc_id{imouse}.number);
        NClass{imouse} = neuroclass{imouse}.UnitID(pfc_id{imouse}.number,1);
        
        Neurons{imouse} = cell(length(clu_pfc),1);
        Q = MakeQfromS(clu_pfc{imouse}, binsize*1e4);
        time = Range(Q, 's');
        tempdat = zscore(full(Data(Q)));
        for icell = 1:length(clu_pfc{imouse})
            Neurons{imouse}{icell} = [time tempdat(:,icell)]; % in sec
        end
    end
    
    %% Calculate PETH
    cnt=1;
    for imouse = 1:length(Dir_Cond.path)
        for icell = 1:length(Neurons{imouse})
            [r,i] = Sync(Neurons{imouse}{icell}, Start(FreezeEpoch{imouse}, 's') ,'durations',[-wi wi]);
            T = SyncMap(r,i,'durations',[-wi wi],'nbins',nbins,'smooth',0);
            OnResp(cnt,:) = mean(T);
            
            [r,i] = Sync(Neurons{imouse}{icell}, End(FreezeEpoch{imouse}, 's') ,'durations',[-wi wi]);
            T = SyncMap(r,i,'durations',[-wi wi],'nbins',nbins,'smooth',0);
            OffResp(cnt,:) = mean(T);
            
            % Calculate auxiliary info
            if IsCalcAdditional
                % Firing rate
                FR(cnt) = GetFiringRate(clu_pfc{imouse}(icell), CondEpoch{imouse});
                % Neuron class
                NeuronClass(cnt) = NClass{imouse}(icell);
            end
            
            cnt = cnt + 1;
            
        end
        
    end
    
    if IsCalcAdditional
        FR(isnan(FR)) = [];
        NeuronClass(isnan(NeuronClass)) = [];
    end
    % Create data
    DatNormZ = [OnResp,OffResp];
    
    %% Calculate accelerometer trace
    Accelero = nan(length(Dir_Cond.path), 201*2);
    for imouse = 1:length(Dir_Cond.path)
        M_st = PlotRipRaw(accelero{imouse},Start(f_cond{imouse}, 's'),2000,1,0,0);
        M_end = PlotRipRaw(accelero{imouse},End(f_cond{imouse}, 's'),2000,1,0,0);
        
        time = [(M_st(:,1))' (M_end(:,1)+4)'];
        Accelero(imouse,:) = [(M_st(:,2))' (M_end(:,2))'];
    end
    
else
    folderdata = ChooseFolderForFigures_DB('Data');
    load([folderdata '/peth_freezing_pfc.mat'], 'DatNormZ', 'time', 'ind_sort', 'PlaceCells', 'SpInfo', 'FR', 'NeuronClass', 'Accelero');
end

%% Plot figure
f1 = figure('units', 'normalized', 'outerposition', [0.1 0.3 0.4 0.8]);
SustVal = nanmean(DatNormZ(:,60:100),2);
% SustVal = nanmean(DatNormZ(:,121:141),2);
% SustVal = nanmean(DatNormZ(:,40:120),2)./abs(nanmean(DatNormZ(:,121:141),2));
UseForTresh = SustVal;
[~,ind_sort] = sort(UseForTresh);
subplot(4,1,1:3)
imagesc(-2:0.05:6,1:size(DatNormZ,1),DatNormZ(ind_sort,:))
caxis([-0.45 0.45])
set(gca, 'XTick', [-2 -1 0 1.5 2.5 4 5 6], 'XTickLabels',...
    {'-2', '-1', '0', '1.5', '-1.5', '0', '1', '2'});
line([0 0],ylim,'color',[.9856, .7372, .2537])
line([4 4],ylim,'color','m')
line(xlim, [50 50], 'Color','w', 'LineStyle', '--', 'LineWidth', 3);
line(xlim, [280 280], 'Color', 'w', 'LineStyle', '--', 'LineWidth', 3)
makepretty
ylabel('Neuron #')
xlabel('Time to freeze on (s)                        Time to freeze off (s)')

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
if IsSaveFig
    foldertosave = ChooseFolderForFigures_DB('Spikes');
    saveas(f1,[foldertosave '/Freezing/PETHonFreezing_pfc.fig']);
    saveFigure(f1, 'PETHonFreezing_pfc', [foldertosave '/Freezing/']);
end


%% Find out how to divide in groups
A = mean(DatNormZ(ind_sort,:));
f2 = figure('units', 'normalized', 'outerposition', [0 0 0.5 0.6]);
plot(mean(DatNormZ(ind_sort,60:100),2))
ylim([-0.6 0.6])
l1 = line(xlim, [mean(A) mean(A)], 'Color', 'b', 'LineWidth', 1);
l2 = line(xlim, [-0.1 -0.1], 'Color', 'r', 'LineWidth', 2);
line([50 50], ylim, 'Color', 'r', 'LineWidth', 2);
l3 = line(xlim, [0.08 0.08], 'Color', 'm', 'LineWidth', 2);
line([280 280], ylim, 'Color', 'm', 'LineWidth', 2);
xlabel('Neuron number')
ylabel('Firing at freezing zscored')
legend([l1 l2 l3], 'Mean', 'OFF-Freezing', 'ON-Freezing', 'Location', 'NorthWest');
makepretty

%% Groups mean
group_data = DatNormZ(ind_sort,:);
group_data = {group_data(1:50,:), group_data(51:279,:), group_data(280:end,:)};
titles = {'OFF on freezing', 'Not affected', 'ON on freezing'};
f3 = figure('units', 'normalized', 'outerposition', [0 0 0.4 1]);
for igroup = 1:length(group_data)
    subplot(3,1,igroup)
    plot(-2.05:0.05:6, mean(group_data{igroup}), 'Color', 'k');
    xlim([-2 6])
    ylim([-0.25 0.2]) 
    line([0 0],ylim,'color',[.9856, .7372, .2537])
    line([4 4],ylim,'color','m')
    ylabel('Firing zscored')
    xlabel('Time to freeze on (s)')
    title(titles{igroup});
    makepretty
end

if IsSaveFig
    foldertosave = ChooseFolderForFigures_DB('Spikes');
    saveas(f3,[foldertosave '/Freezing/SpikesOnFreezing_pfc_groups.fig']);
    saveFigure(f3, 'SpikesOnFreezing_pfc_groups', [foldertosave '/Freezing/']);
end

%% Try first 3 PCs 
[~, sc, ~,~,explained] = pca(DatNormZ');
f4 = figure('units', 'normalized', 'outerposition', [0 0 0.4 1]);
for ipc = 1:3
    subplot(3,1,ipc)
    plot(-2.05:0.05:6, sc(:, ipc), 'Color', 'k');
    xlim([-2 6])
    ylim([-2 4.5])
    line([0 0],ylim,'color',[.9856, .7372, .2537])
    line([4 4],ylim,'color','m')
    ylabel(['PC#' num2str(ipc)])
    xlabel('Time to freeze on (s)')
    title(['PC#' num2str(ipc) ', ' num2str(round(explained(ipc),1)) '% of data explained']);
    makepretty
end

if IsSaveFig
    foldertosave = ChooseFolderForFigures_DB('Spikes');
    saveas(f4,[foldertosave '/Freezing/PCAOnFreezing_pfc_groups.fig']);
    saveFigure(f4, 'PCAOnFreezing_pfc_groups', [foldertosave '/Freezing/']);
end

%% Plot
if (Redo && IsCalcAdditional) || (~Redo && exist('PlaceCells', 'var'))
    FR_sorted = FR(ind_sort);
    NeuronClass_sorted = NeuronClass(ind_sort);
    
    % Figures specs
    cols = {[0.2 0.2 0.8], [.6 .6 .6], [.8 .2 .2]};
    legs = {'OFF', '~~', 'ON'};
    
    f5 = figure('units', 'normalized', 'outerposition', [0 0 0.6 0.5]);
    subplot(121)
    MakeBoxPlot_DB({FR_sorted(1:50), FR_sorted(51:279), FR_sorted(280:end)}, cols, 1:3, legs, 0);
    ylabel('Firing rate (Hz)')
    title('Firing rate')
    makepretty
    subplot(122)
    NeuronClass_OFF = NeuronClass_sorted(1:50);
    NeuronClass_N = NeuronClass_sorted(51:279);
    NeuronClass_ON = NeuronClass_sorted(280:end);
    Int_perc = [sum(NeuronClass_OFF<0)/length(NeuronClass_OFF)*100 sum(NeuronClass_N<0)/length(NeuronClass_N)*100 ...
        sum(NeuronClass_ON<0)/length(NeuronClass_ON)*100];
    [~,h] = PlotErrorBarN_DB(Int_perc, 'barcolors', [1 0 0], 'barwidth', 0.6, 'newfig', 0, 'showpoints', 0);
    h.FaceColor = 'flat';
    h.CData(1,:) = [0 0 1];
    h.CData(2,:) = [.8 .8 .8];
    set(gca, 'XTick', 1:3, 'XTickLabels', {'OFF', '~~', 'ON'});
    ylabel('% of interneurons in subpop')
    title('Neuron class')
    makepretty
    
    if IsSaveFig
        foldertosave = ChooseFolderForFigures_DB('Spikes');
        saveas(f5,[foldertosave '/Freezing/Freezing_pfc_group_props.fig']);
        saveFigure(f5, 'Freezing_pfc_group_props', [foldertosave '/Freezing/']);
    end
end
%% SaveFile
if IsSaveFile
    folderdata = ChooseFolderForFigures_DB('Data');
    save([folderdata '/peth_freezing_pfc.mat'], 'DatNormZ', 'time', 'ind_sort', 'PlaceCells', 'SpInfo', 'FR', 'NeuronClass', 'Accelero');
end

end