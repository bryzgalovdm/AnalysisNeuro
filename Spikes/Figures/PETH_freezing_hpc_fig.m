function f1 = PETH_freezing_hpc_fig(varargin)
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
%                                   everything from scratch (default  =false)
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
Redo = 0;
IsCalcAdditional = 0;

mice = [797 798 828 861 882 905 906 911 912 977 994 1117 1124 1161 1162 1168 1182 1186 1199];

wi = 2; % in sec
binsize = 0.05; % in sec
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
    Dir = PathForExperimentsERC_Dima('UMazePAG');
    Dir = RestrictPathForExperiment(Dir,'nMice', mice);
    
    % Allocate arrays
    b = cell(length(Dir.path),1); % behavior
    s = cell(length(Dir.path),1); % spikes
    CondEpoch = cell(length(Dir.path),1);
    UMazeEpoch = cell(length(Dir.path),1);
    FreezeEpoch = cell(length(Dir.path),1);
    Neurons = cell(length(Dir.path),1);
    PlaceCells = nan(5e4,1);
    PlaceCells_SZ = nan(5e4,1);
    SpInfo = zeros(5e4,1);
    FR = nan(5e4,1);
    NeuronClass = nan(5e4,1);
    kappa_theta = nan(5e4,1);
    kappa4Hz = nan(5e4,1);
    pval_theta = nan(5e4,1);
    pval4Hz = nan(5e4,1);    
    
    for imouse = 1:length(Dir.path)
        b{imouse} = load([Dir.path{imouse}{1} '/behavResources.mat'], 'SessionEpoch', 'FreezeAccEpoch', 'MovAcctsd', ...
            'AlignedXtsd', 'AlignedYtsd', 'Vtsd', 'TTLInfo');
        s{imouse} = load([Dir.path{imouse}{1} '/SpikeData.mat'], 'S', 'BasicNeuronInfo', 'PlaceCells', 'TT', 'RippleGroups');
    end
    
    %% Organize epochs and data
    % Epochs
    for imouse = 1:length(Dir.path)
        [~, ~, ~, CondEpoch{imouse}] = ReturnMnemozyneEpochs(b{imouse}.SessionEpoch);
        [~, ~, UMazeEpoch{imouse}] = ReturnMnemozyneEpochs(b{imouse}.SessionEpoch,...
            'Speed', b{i}.Vtsd, 'SpeedThresh', speed_thresh);
        
        FreezeEpoch{imouse} = and(b{imouse}.FreezeAccEpoch,CondEpoch{imouse});
        FreezeEpoch{imouse} = dropShortIntervals(FreezeEpoch{imouse}, 4*1E4);
        FreezeEpoch{imouse} = mergeCloseIntervals(FreezeEpoch{imouse},1*1E4);
        % Remove epochs co-terminated with stims
        stims = ts(Start(b{imouse}.TTLInfo.StimEpoch));
        stims = Restrict(stims, CondEpoch{imouse});
        stims = intervalSet(Range(stims), Range(stims)+100);
        FreezeEpoch{imouse} = GetStandaloneFreezeEpochs(FreezeEpoch{imouse}, stims, 2.5);
        
    end
    
    
    
    % Data
    for imouse = 1:length(Dir.path)
        % Get cells from pyramidal layer
        id_ripples = PyramLayerSorting(s{imouse}.RippleGroups, s{imouse}.TT);
        selected_cells{imouse} = s{imouse}.S(id_ripples);
        
        Neurons{imouse} = cell(length(selected_cells{imouse}),1);
        Q = MakeQfromS(selected_cells{imouse}, binsize*1e4);
        time = Restrict(Q, CondEpoch{imouse});
        tempdat = zscore(full(Data(Restrict(Q, CondEpoch{imouse}))));
        for icell = 1:length(selected_cells{imouse})
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
            
            % Calculate auxiliary info
            if IsCalcAdditional
                % Spatial Info
                if ~isempty(s{imouse}.PlaceCells.stats{icell})
                    if ~isempty(s{imouse}.PlaceCells.stats{icell}.spatialInfo)
                        SpInfo(cnt) = s{imouse}.PlaceCells.stats{icell}.spatialInfo;
                    else
                        SpInfo(cnt) = nan;
                    end
                else
                   SpInfo(cnt) = nan; 
                end
                % Place Cell or not
                if sum(s{imouse}.PlaceCells.idx == icell) > 0
                    PlaceCells(cnt) = true;
                else
                    PlaceCells(cnt) = false;
                end
                if sum(s{imouse}.PlaceCells.SZ == icell) > 0
                    PlaceCells_SZ(cnt) = true;
                else
                    PlaceCells_SZ(cnt) = false;
                end
                % Firing rate
                FR(cnt) = s{imouse}.BasicNeuronInfo.firingrate(icell);
                % Neuron class
                NeuronClass(cnt) = s{imouse}.BasicNeuronInfo.neuroclass(icell);
                % Theta kappa
                kappa_theta(cnt) = s{imouse}.BasicNeuronInfo.kappatheta.Transf(icell);
                % 4Hz kappa
                if isfield(s{imouse}.BasicNeuronInfo, 'kappa4Hz')
                    kappa4Hz(cnt) = s{imouse}.BasicNeuronInfo.kappa4Hz.Transf(icell);
                end
                % Theta pval
                pval_theta(cnt) = s{imouse}.BasicNeuronInfo.pvaltheta.Transf(icell);
                % 4Hz pval
                if isfield(s{imouse}.BasicNeuronInfo, 'pval4Hz')
                    pval4Hz(cnt) = s{imouse}.BasicNeuronInfo.pval4Hz.Transf(icell);
                end
            end
            
            cnt = cnt + 1;
            
        end
        
    end
    
    if IsCalcAdditional
        SpInfo = nonzeros(SpInfo);
        PlaceCells(isnan(PlaceCells)) = [];
        PlaceCells_SZ(isnan(PlaceCells_SZ)) = [];
        FR(isnan(FR)) = [];
        ids_all = length(FR);
        NeuronClass(isnan(NeuronClass)) = [];
        kappa_theta(ids_all+1:end) = [];
        kappa4Hz(ids_all+1:end) = [];
        pval_theta(isnan(pval_theta)) = [];
        pval4Hz(ids_all+1:end) = [];
    end
    % Create data
    DatNormZ = [OnResp,OffResp];
    
    %% Calculate accelerometer trace
    Accelero = nan(length(Dir.path), 201*2);
    for imouse = 1:length(Dir.path)
        M_st = PlotRipRaw(b{imouse}.MovAcctsd,Start(FreezeEpoch{imouse},'s'),2000,1,0,0);
        M_end = PlotRipRaw(b{imouse}.MovAcctsd,End(FreezeEpoch{imouse},'s'),2000,1,0,0);
        
        time = [(M_st(:,1))' (M_end(:,1)+4)'];
        Accelero(imouse,:) = [(M_st(:,2))' (M_end(:,2))'];
    end
    
else
    folderdata = ChooseFolderForFigures_DB('Data');
    load([folderdata '/peth_freezing.mat'], 'DatNormZ', 'time', 'ind_sort', 'PlaceCells', 'SpInfo', 'FR', 'NeuronClass', 'Accelero');
end

%% Plot figure
f1 = figure('units', 'normalized', 'outerposition', [0.1 0.3 0.35 0.85]);
SustVal = nanmean(DatNormZ(:,60:100),2);
% SustVal = nanmean(DatNormZ(:,121:141),2);
% SustVal = nanmean(DatNormZ(:,40:120),2)./abs(nanmean(DatNormZ(:,121:141),2));
UseForTresh = SustVal;
[~,ind_sort] = sort(UseForTresh);
subplot(3,1,1:2)
imagesc(-2:0.05:6,1:size(DatNormZ,1),DatNormZ(ind_sort,:))
caxis([-0.45 0.45])
set(gca, 'XTick', [-2 -1 0 1.5 2.5 4 5 6], 'XTickLabels',...
    {'-2', '-1', '0', '1.5', '-1.5', '0', '1', '2'});
line([0 0],ylim,'color',[.9856, .7372, .2537])
line([4 4],ylim,'color','m')
line(xlim, [90 90], 'Color','w', 'LineStyle', '--', 'LineWidth', 3);
line(xlim, [630 630], 'Color', 'w', 'LineStyle', '--', 'LineWidth', 3)
makepretty
ylabel('Neuron #')
% xlabel('Time to freeze on (s)                        Time to freeze off (s)')

% Accelero meter
subplot(3,1,3)
if size(Accelero,1) > 1
    shadedErrorBar(time,mean(Accelero, 1), std(Accelero, 1), 'k')
else
    plot(time,Accelero,'color','k')
end
set(gca, 'XTick', [-2 -1 0 1.5 2.5 4 5 6], 'XTickLabels',...
    {'-2', '-1', '0', '1.5', '-1.5', '0', '1', '2'});
makepretty
ylabel('Accelero')
% xlabel('Time to freeze on (s)')

% Save figure
if IsSaveFig
    foldertosave = ChooseFolderForFigures_DB('Spikes');
    saveas(f1,[foldertosave '/Freezing/PETHonFreezing_hpc.fig']);
    saveFigure(f1, 'PETHonFreezing_hpc', [foldertosave '/Freezing/']);
end


%% Find out how to divide in groups
A = mean(DatNormZ(ind_sort,:));
f2 = figure('units', 'normalized', 'outerposition', [0 0 0.5 0.6]);
plot(mean(DatNormZ(ind_sort,60:100),2))
ylim([-0.6 0.6])
l1 = line(xlim, [mean(A) mean(A)], 'Color', 'b', 'LineWidth', 1);
l2 = line(xlim, [-0.1 -0.1], 'Color', 'r', 'LineWidth', 2);
line([100 100], ylim, 'Color', 'r', 'LineWidth', 2);
l3 = line(xlim, [0.08 0.08], 'Color', 'm', 'LineWidth', 2);
line([600 600], ylim, 'Color', 'm', 'LineWidth', 2);
xlabel('Neuron number')
ylabel('Firing at freezing zscored')
legend([l1 l2 l3], 'Mean', 'OFF-Freezing', 'ON-Freezing', 'Location', 'NorthWest');
makepretty

%% Groups mean
group_data = DatNormZ(ind_sort,:);
group_data = {group_data(1:70,:), group_data(711:729,:), group_data(730:end,:)};
titles = {'OFF on freezing', 'Not affected', 'ON on freezing'};
f3 = figure('units', 'normalized', 'outerposition', [0 0 0.4 1]);
for igroup = 1:length(group_data)
    subplot(3,1,igroup)
    plot(-2.05:0.05:6, mean(group_data{igroup}), 'Color', 'k');
    xlim([-2 6])
    ylim([-0.25 0.3]) 
    line([0 0],ylim,'color',[.9856, .7372, .2537])
    line([4 4],ylim,'color','m')
    set(gca, 'XTick', [-2 -1 0 1.5 2.5 4 5 6], 'XTickLabels',...
    {'-2', '-1', '0', '1.5', '-1.5', '0', '1', '2'});
    ylabel('Firing zscored')
%     xlabel('Time to freeze on (s)')
    title(titles{igroup});
    makepretty
end

if IsSaveFig
    foldertosave = ChooseFolderForFigures_DB('Spikes');
    saveas(f3,[foldertosave '/Freezing/SpikesOnFreezing_hpc_groups.fig']);
    saveFigure(f3, 'SpikesOnFreezing_hpc_groups', [foldertosave '/Freezing/']);
end

%% Try first 3 PCs 
[~, sc, ~,~,explained] = pca(DatNormZ');
f4 = figure('units', 'normalized', 'outerposition', [0 0 0.4 1]);
for ipc = 1:3
    subplot(3,1,ipc)
    plot(-2.05:0.05:6, sc(:, ipc), 'Color', 'k');
    xlim([-2 6])
    ylim([-2.5 4.5])
    line([0 0],ylim,'color',[.9856, .7372, .2537])
    line([4 4],ylim,'color','m')
    set(gca, 'XTick', [-2 -1 0 1.5 2.5 4 5 6], 'XTickLabels',...
    {'-2', '-1', '0', '1.5', '-1.5', '0', '1', '2'});
    ylabel(['PC#' num2str(ipc)])
%     xlabel('Time to freeze on (s)')
    title(['PC#' num2str(ipc) ', ' num2str(round(explained(ipc),1)) '% of data explained']);
    makepretty
end

if IsSaveFig
    foldertosave = ChooseFolderForFigures_DB('Spikes');
    saveas(f4,[foldertosave '/Freezing/PCAOnFreezing_hpc_groups.fig']);
    saveFigure(f4, 'PCAOnFreezing_hpc_groups', [foldertosave '/Freezing/']);
end

%% Plot
if (Redo && IsCalcAdditional) || (~Redo && exist('PlaceCells', 'var'))
    SpInfo_sorted = SpInfo(ind_sort);
    FR_sorted = FR(ind_sort);
    PlaceCells_sorted = PlaceCells(ind_sort);
    PlaceCellsSZ_sorted = PlaceCells_SZ(ind_sort);
    NeuronClass_sorted = NeuronClass(ind_sort);
    kappa_theta_sorted = kappa_theta(ind_sort);
    kappa4Hz_sorted = kappa4Hz(ind_sort);
    pval_theta_sorted = pval_theta(ind_sort);
    pval4Hz_sorted = pval4Hz(ind_sort);
    
    % Figures specs
    cols = {[0.2 0.2 0.8], [.6 .6 .6], [.8 .2 .2]};
    legs = {'OFF', '~~', 'ON'};
    
    f5 = figure('units', 'normalized', 'outerposition', [0 0 1 0.5]);
    subplot(161)
    MakeBoxPlot_DB({SpInfo_sorted(1:70), SpInfo_sorted(71:729), SpInfo_sorted(730:end)}, cols, 1:3, legs, 0);
    ylabel('Spatial info')
    title('Spatial information')
    makepretty
    subplot(162)
    PlaceCells_OFF = PlaceCells_sorted(1:70);
    PlaceCells_N = PlaceCells_sorted(71:729);
    PlaceCells_ON = PlaceCells_sorted(730:end);
    PC_perc = [sum(PlaceCells_OFF==1)/length(PlaceCells_OFF)*100 sum(PlaceCells_N==1)/length(PlaceCells_N)*100 ...
        sum(PlaceCells_ON==1)/length(PlaceCells_ON)*100];
    [~,h] = PlotErrorBarN_DB(PC_perc, 'barcolors', [1 0 0], 'barwidth', 0.6, 'newfig', 0, 'showpoints', 0);
    h.FaceColor = 'flat';
    h.CData(1,:) = [0 0 1];
    h.CData(2,:) = [.8 .8 .8];
    set(gca, 'XTick', 1:3, 'XTickLabels', {'OFF', '~~', 'ON'});
    ylabel('% of place cells in subpop')
    title('Place cells')
    makepretty
    subplot(163)
    PlaceCells_OFF = PlaceCells_sorted(1:70);
    PlaceCells_N = PlaceCells_sorted(71:729);
    PlaceCells_ON = PlaceCells_sorted(730:end);
    PC_perc = [sum(PlaceCells_OFF==1)/sum(PlaceCells)*100 sum(PlaceCells_N==1)/sum(PlaceCells)*100 ...
        sum(PlaceCells_ON==1)/sum(PlaceCells)*100];
    [~,h] = PlotErrorBarN_DB(PC_perc, 'barcolors', [1 0 0], 'barwidth', 0.6, 'newfig', 0, 'showpoints', 0);
    h.FaceColor = 'flat';
    h.CData(1,:) = [0 0 1];
    h.CData(2,:) = [.8 .8 .8];
    set(gca, 'XTick', 1:3, 'XTickLabels', {'OFF', '~~', 'ON'});
    ylabel('% of place cells in subpop of all place cells')
    title('Place cells')
    makepretty
    subplot(164)
    PlaceCellsSZ_OFF = PlaceCellsSZ_sorted(1:70);
    PlaceCellsSZ_N = PlaceCellsSZ_sorted(71:729);
    PlaceCellsSZ_ON = PlaceCellsSZ_sorted(730:end);
    PC_percSZ = [sum(PlaceCellsSZ_OFF==1)/sum(PlaceCells_SZ)*100 sum(PlaceCellsSZ_N==1)/sum(PlaceCells_SZ)*100 ...
        sum(PlaceCellsSZ_ON==1)/sum(PlaceCells_SZ)*100];
    [~,h] = PlotErrorBarN_DB(PC_percSZ, 'barcolors', [1 0 0], 'barwidth', 0.6, 'newfig', 0, 'showpoints', 0);
    h.FaceColor = 'flat';
    h.CData(1,:) = [0 0 1];
    h.CData(2,:) = [.8 .8 .8];
    set(gca, 'XTick', 1:3, 'XTickLabels', {'OFF', '~~', 'ON'});
    ylabel('% of SZ place cells in subpop of all place cells')
    title('SZ Place cells')
    makepretty
    subplot(165)
    MakeBoxPlot_DB({FR_sorted(1:70), FR_sorted(71:729), FR_sorted(730:end)}, cols, 1:3, legs, 0);
    ylabel('Firing rate (Hz)')
    title('Firing rate')
    makepretty
    subplot(166)
    NeuronClass_OFF = NeuronClass_sorted(1:70);
    NeuronClass_N = NeuronClass_sorted(71:729);
    NeuronClass_ON = NeuronClass_sorted(730:end);
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
        saveas(f5,[foldertosave '/Freezing/Freezing_hpc_group_props.fig']);
        saveFigure(f5, 'Freezing_hpc_group_props', [foldertosave '/Freezing/']);
    end

    f6 = figure('units', 'normalized', 'outerposition', [0 0 0.6 0.5]);
    subplot(121)
    MakeBoxPlot_DB({FR_sorted(1:70), FR_sorted(71:729), FR_sorted(730:end)}, cols, 1:3, legs, 0);
    ylabel('Firing rate (Hz)')
    title('Firing rate')
    makepretty
    subplot(122)
    NeuronClass_OFF = NeuronClass_sorted(1:70);
    NeuronClass_N = NeuronClass_sorted(71:729);
    NeuronClass_ON = NeuronClass_sorted(730:end);
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
        saveas(f6,[foldertosave '/Freezing/Freezing_hpc_group_props_short.fig']);
        saveFigure(f6, 'Freezing_hpc_group_props_short', [foldertosave '/Freezing/']);
    end
    
    f7 = figure('units', 'normalized', 'outerposition', [0 0 1 0.5]);
    subplot(141)
    MakeBoxPlot_DB({kappa_theta_sorted(1:70), kappa_theta_sorted(71:729), kappa_theta_sorted(730:end)}, cols, 1:3, legs, 0);
    ylabel('Kappa theta')
    title('Kappa theta')
    makepretty
    subplot(142)
    MakeBoxPlot_DB({kappa4Hz_sorted(1:70), kappa4Hz_sorted(71:729), kappa4Hz_sorted(730:end)}, cols, 1:3, legs, 0);
    ylabel('Kappa 4Hz')
    title('Kappa 4Hz')
    makepretty
    subplot(143)
    pval_theta_OFF = pval_theta_sorted(1:70);
    pval_theta_N = pval_theta_sorted(71:729);
    pval_theta_ON = pval_theta_sorted(730:end);
    Int_perc = [sum(pval_theta_OFF<0.05)/length(pval_theta_OFF)*100 sum(pval_theta_N<0.05)/length(pval_theta_N)*100 ...
        sum(pval_theta_ON<0.05)/length(pval_theta_ON)*100];
    [~,h] = PlotErrorBarN_DB(Int_perc, 'barcolors', [1 0 0], 'barwidth', 0.6, 'newfig', 0, 'showpoints', 0);
    h.FaceColor = 'flat';
    h.CData(1,:) = [0 0 1];
    h.CData(2,:) = [.8 .8 .8];
    set(gca, 'XTick', 1:3, 'XTickLabels', {'OFF', '~~', 'ON'});
    ylabel('% of cell sign. modulated by theta')
    title('Theta modulation')
    makepretty
    subplot(144)
    pval4Hz_OFF = pval4Hz_sorted(1:70);
    pval4Hz_N = pval4Hz_sorted(71:729);
    pval4Hz_ON = pval4Hz_sorted(730:end);
    Int_perc = [sum(pval4Hz_OFF<0.05)/length(pval4Hz_OFF)*100 sum(pval4Hz_N<0.05)/length(pval4Hz_N)*100 ...
        sum(pval4Hz_ON<0.05)/length(pval4Hz_ON)*100];
    [~,h] = PlotErrorBarN_DB(Int_perc, 'barcolors', [1 0 0], 'barwidth', 0.6, 'newfig', 0, 'showpoints', 0);
    h.FaceColor = 'flat';
    h.CData(1,:) = [0 0 1];
    h.CData(2,:) = [.8 .8 .8];
    set(gca, 'XTick', 1:3, 'XTickLabels', {'OFF', '~~', 'ON'});
    ylabel('% of cell sign. modulated by 4Hz')
    title('4Hz modulation')
    makepretty
    if IsSaveFig
        foldertosave = ChooseFolderForFigures_DB('Spikes');
        saveas(f7,[foldertosave '/Freezing/Freezing_hpc_group_spectralmod.fig']);
        saveFigure(f7, 'Freezing_hpc_group_spectralmod', [foldertosave '/Freezing/']);
    end

end
%% SaveFile
if IsSaveFile
    folderdata = ChooseFolderForFigures_DB('Data');
    save([folderdata '/peth_freezing.mat'], 'DatNormZ', 'time', 'ind_sort', 'PlaceCells', 'SpInfo', 'FR', 'NeuronClass', 'Accelero');
end

end