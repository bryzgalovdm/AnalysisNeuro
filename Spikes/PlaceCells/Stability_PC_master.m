function Stability_PC_master(mice, expe, varargin)


%% Parameters
sav=0;
% Sleep time to restrict
SleepTimeToRestrict = 2*60*60*1e4; % 2 hours
% Speed thresh
speed_thresh = 6;
% Spatial information threshold
si_thresh = 1.2;

sizemap = 50; % Default 50

% Names of epochs
NamesEpochs = {'BaselineExplo', '1stHalfBaseline', '2ndHalfBaseline', 'Conditioning', 'PostTests', 'Extinction'};

% Maze borders
mazeMap = [6 7; 6 59; 59 59; 59 7; 39 7; 39 42; 24 42; 24 7; 6 7];
ShockZoneMap = [6 7; 6 30; 24 30; 24 7; 6 7];
mazeSpikes = [0 0; 0 1; 1 1; 1 0; 0.65 0; 0.65 0.75; 0.35 0.75; 0.35 0; 0 0];
shockZoneSpikes = [0 0; 0 0.35; 0.35 0.35; 0.35 0; 0 0]; 

%% Optional Arguments
for i=1:2:length(varargin)    
    switch(lower(varargin{i})) 
        case 'issave'
            sav = varargin{i+1};
            if sav ~= 1 && sav ~= 0
                error('Incorrect value for property ''IsSave'' (type ''help PreTestCharacteristics'' for details).');
            end
    end
end

%% Optional Arguments
for i=1:2:length(varargin)    
    switch(lower(varargin{i})) 
        case 'issave'
            sav = varargin{i+1};
            if sav ~= 1 && sav ~= 0
                error('Incorrect value for property ''IsSave'' (type ''help PreTestCharacteristics'' for details).');
            end
        case 'showpoints'
            SP = varargin{i+1};
            if SP ~= 1 && SP ~= 0
                error('Incorrect value for property ''ShowPoints'' (type ''help PreTestCharacteristics'' for details).');
            end
    end
end

%% Manage experiment
if strcmp(expe, 'PAG')
    fetchpaths = 'UMazePAG';
    FigName = 'RipplesPAG';
elseif strcmp(expe, 'MFB')
    fetchpaths = 'StimMFBWake';
    FigName = 'RipplesMFB';
end

%% Allocate data arrays
% Get paths of each individual mousesum(Stability{1}>0.5)/length(Stability{1})
Dir = PathForExperimentsERC(fetchpaths);
Dir = RestrictPathForExperiment(Dir,'nMice',mice);
numsessions = CountNumSesionsERC(Dir);

% Data
Session = cell(length(numsessions), 1);
% Data
Neurons = cell(length(numsessions), 1);
%% Load the data

cnt=1;
for imouse=1:length(Dir.path)
    for isession = 1:length(Dir.path{imouse})
        Session{cnt} = load([Dir.path{imouse}{isession} 'behavResources.mat'], 'SessionEpoch', 'Vtsd', 'AlignedXtsd',...
            'AlignedYtsd');
        Neurons{cnt} = load([Dir.path{imouse}{isession} 'SpikeData.mat'],'S','PlaceCells', 'TT', 'BasicNeuronInfo');
        MouseNumber{cnt} = Dir.name{imouse};
        cnt=cnt+1;
    end
end

%% Allocate memory
Epochs = cell(numsessions, 1);
PlaceCells = cell(numsessions, 1);

%% Create epochs
for isession = 1:numsessions
    % Speed-filtered epochs
    [~, ~, UMazeEpoch, CondEpoch, ~, PostTestEpoch] =...
        ReturnMnemozyneEpochs(Session{isession}.SessionEpoch, 'Speed', Session{isession}.Vtsd,...
         'SpeedThresh', speed_thresh);
    % Speedsum(Stability{1}>0.5)/length(Stability{1})
    smoothspeed  = tsd(Range(Session{isession}.Vtsd),movmedian(Data(Session{isession}.Vtsd),5));
    LocomotionEpoch = thresholdIntervals(smoothspeed,speed_thresh,'Direction','Above');
    % Half of baseline exploration
    [~, ~, TempBaseEpoch, ~, ~, ~] = ReturnMnemozyneEpochs(Session{isession}.SessionEpoch);
    st = Start(TempBaseEpoch);
    en = End(TempBaseEpoch);
    Base1Half = intervalSet(st,((en-st)/2)+st);
    Base2Half = intervalSet(((en-st)/2)+st,en);
    MovingBase1Half = and(LocomotionEpoch,Base1Half);
    MovingBase2Half = and(LocomotionEpoch,Base2Half);
    if isfield(Session{isession}.SessionEpoch, 'Extinct')
        MovingExtinction = and(LocomotionEpoch, Session{isession}.SessionEpoch.Extinct);
    else
        MovingExtinction = [];
    end
    
    Epochs{isession} = {UMazeEpoch, MovingBase1Half, MovingBase2Half, CondEpoch, PostTestEpoch, MovingExtinction};
end    

%% Extract place cells
cnt=1;
for isession = 1:numsessions
    if ~isempty(Neurons{isession}.PlaceCells.idx)
        for ineuron = 1:length(Neurons{isession}.PlaceCells.idx)
            for iepoch = 1:length(Epochs{isession})
                if ~isempty(Epochs{isession}{iepoch})
                    [map{cnt}{iepoch}, ~, stats{cnt}{iepoch}, px{cnt}{iepoch}, py{cnt}{iepoch},...
                        FR{cnt}{iepoch}, ~, ~, lEpoch{cnt}{iepoch}] = PlaceField_DB(Neurons{isession}.S{Neurons{isession}.PlaceCells.idx(ineuron)},...
                        Session{isession}.AlignedXtsd, Session{isession}.AlignedYtsd, 'Epoch', Epochs{isession}{iepoch},...
                        'plotresults', 0, 'plotpoisson', 0, 'SizeMap', sizemap);
                else
                    map{cnt}{iepoch} = [];
                    stats{cnt}{iepoch} = [];
                    px{cnt}{iepoch} = [];
                    py{cnt}{iepoch} = [];
                    FR{cnt}{iepoch} = [];
                    lEpoch{cnt}{iepoch} = [];
                end
            end
            if sum(Neurons{isession}.PlaceCells.SZ == Neurons{isession}.PlaceCells.idx(ineuron)) > 0
                isSZ(cnt) = true;
            else
                isSZ(cnt) = false;
            end
            traj_id(cnt) = isession;
            ClNum{cnt} = Neurons{isession}.TT{Neurons{isession}.PlaceCells.idx(ineuron)};
            FR_base(cnt) = Neurons{isession}.BasicNeuronInfo.firingrate(Neurons{isession}.PlaceCells.idx(ineuron));
            cnt=cnt+1;
        end
    end
end

% leangth cond
for icell = 1:length(map)
    LL{icell} = lEpoch{icell}{4};
end
LCond = unique(cell2mat(LL));



%% Calculate the stability metrics
Explostab_epochs = cell(numsessions, 1);
for isession = 1:length(Epochs)
    stab_epochs{isession} = Epochs(1:4);
end
Stability = cell(2,1); % Within hab-stability and Hab-Cond stability


for icell = 1:length(map)
    for istab = 1:length(Stability)
        if ~isempty(map{icell}{2}) && ~isempty(map{icell}{3})
            cf = corrcoef(map{icell}{2}.rate, map{icell}{3}.rate); % Within hab-stability
            Stability{1}(icell) = cf(2,1);
        end
        if ~isempty(map{icell}{1}) && ~isempty(map{icell}{4})
            cf = corrcoef(map{icell}{1}.rate, map{icell}{4}.rate); % Hab-Cond-Stability
            Stability{2}(icell) = cf(2,1);
        end
    end
end
% Find really nice place cells
for icell=1:length(stats)
    if stats{icell}{1}.spatialInfo>si_thresh
        isGood(icell) = true;
    else
        isGood(icell) = false;
    end
end
% ge tsi in a vector
for icell=1:length(stats)
    if ~isempty(stats{icell}{1}.spatialInfo)
        si(icell) = stats{icell}{1}.spatialInfo;
    else
        si(icell)=nan;
    end
end
[rankedECStability, idranked] = sort(Stability{2}, 'descend');
startrank = find(rankedECStability<0.5,1, 'first');


% Plot it
cols = {[.2 .2 .2], [.9 0 0]};
xticks = {'Explo', 'Explo-Cond'};
fs = figure ('units', 'normalized','outerposition', [0 0 1 1]);
subplot(131)
b = MakeBoxPlot_DB(Stability, cols, [1 2], xticks, 0);
ylabel('Rate map correlation');
title(['All place cells,N=' num2str(length(map))])
for iplot = 1:length(b)
    b{iplot}.boxAlpha = .45;
end
p = DoWilcoxonOnArray(Stability, {[1 2]});
if p < 0.05
    sigstar_DB({[1 2]},p,0,'LineWigth',16,'StarSize',24);
end
makepretty
subplot(132)
Stability_noSZ = Stability;
Stability_noSZ{1}(isSZ) = [];
Stability_noSZ{2}(isSZ) = [];
b = MakeBoxPlot_DB(Stability_noSZ, cols, [1 2], xticks, 0);
ylabel('Rate map correlation');
title(['Place cells outside of SZ, N=' num2str(sum(~isSZ))])
for iplot = 1:length(b)
    b{iplot}.boxAlpha = .45;
end
p = DoWilcoxonOnArray(Stability_noSZ, {[1 2]});
if p < 0.05
    sigstar_DB({[1 2]},p,0,'LineWigth',16,'StarSize',24);
end
makepretty
subplot(133)
Stability_noSZ_highSI = Stability;
Stability_noSZ_highSI{1}(or(isSZ, ~isGood)) = [];
Stability_noSZ_highSI{2}(or(isSZ, ~isGood)) = [];
b = MakeBoxPlot_DB(Stability_noSZ_highSI, cols, [1 2], xticks, 0);
ylabel('Rate map correlation');
title(['Very good place cells outside of SZ, N=' num2str(sum(~or(isSZ, ~isGood)))])
for iplot = 1:length(b)
    b{iplot}.boxAlpha = .45;
end
p = DoWilcoxonOnArray(Stability_noSZ_highSI, {[1 2]});
if p < 0.05
    sigstar_DB({[1 2]},p,0,'LineWigth',16,'StarSize',24);
end
makepretty

if sav
    foldertochoose = ChooseFolderForFigures_DB('PlaceFieldsStability');
    saveas(fs,[foldertochoose filesep 'Stability_average.fig']);
    saveFigure(fs, 'Stability_average', foldertochoose);
end

fs1 = figure ('units', 'normalized','outerposition', [0 0 .35 .65]);
h1 = histogram(Stability{1}, 30, 'FaceColor', [0.2 .2 .2], 'FaceAlpha', .5, 'EdgeAlpha',.5);
hold on
h2 = histogram(Stability{2}, 30, 'FaceColor', [0.8 0 0], 'FaceAlpha', .3, 'EdgeAlpha',.3);
xlabel('Stability coefficient')
ylabel('Count')
legend([h1 h2], {'Within exploration stability', 'Exploration-Learning stability'}, 'Location', 'NorthWest')
makepretty

if sav
    foldertochoose = ChooseFolderForFigures_DB('PlaceFieldsStability');
    saveas(fs1,[foldertochoose filesep 'Stability_hist.fig']);
    saveFigure(fs1, 'Stability_hist', foldertochoose);
end

fs2 = figure ('units', 'normalized','outerposition', [0 0 1 .65]);
subplot(121)
scatter(si, Stability{1}, 100, 'k', 'filled');
[rho,p] = corr(si', Stability{1}', 'Type', 'Spearman', 'rows', 'complete');
l1 = lsline;
l1.Color = 'k';
t1 = text(.75, .4, ['rho = ' num2str(round(rho, 3))], 'sc', 'FontSize', 15);
t2 = text(.75, .3, ['p = ' num2str(round(p, 3))], 'sc', 'FontSize', 15);
xlabel('Spatial info')
ylabel('Explo stability')
title('Explo stability')
makepretty
subplot(122)
idnan = isnan(Stability{2});
sin = si(~idnan);
Stabilityn = Stability{2}(~idnan);
scatter(sin, Stabilityn, 100, 'k', 'filled');
[rho,p] = corr(sin', Stabilityn', 'Type', 'Spearman', 'rows', 'complete');
l1 = lsline;
l1.Color = 'k';
t1 = text(.75, .4, ['rho = ' num2str(round(rho, 3))], 'sc', 'FontSize', 15);
t2 = text(.75, .3, ['p = ' num2str(round(p, 3))], 'sc', 'FontSize', 15);
xlabel('Spatial info')
ylabel('Explo-Cond stability')
title('Explo-Cond stability')
makepretty

if sav
    foldertochoose = ChooseFolderForFigures_DB('PlaceFieldsStability');
    saveas(fs2,[foldertochoose filesep 'corr_si_stab.fig']);
    saveFigure(fs2, 'corr_si_stab', foldertochoose);
end

%% Remove all cells with less than a minute
time_thresh = 60;
for icell = 1:length(LL)
    if isempty(LL{icell}) || LL{icell} < time_thresh
        id_g(icell) = false;
    else
        id_g(icell) = true;
    end
end

fs = figure ('units', 'normalized','outerposition', [0 0 .4 1]);
Stability_Lepoch = Stability;
Stability_Lepoch{1}(id_g) = [];
Stability_Lepoch{2}(id_g) = [];
b = MakeBoxPlot_DB(Stability_Lepoch, cols, [1 2], xticks, 0);
ylabel('Rate map correlation');
title(['All place cells,N=' num2str(length(Stability_Lepoch{1}))])
for iplot = 1:length(b)
    b{iplot}.boxAlpha = .45;
end
p = DoWilcoxonOnArray(Stability_Lepoch, {[1 2]});
if p < 0.05
    sigstar_DB({[1 2]},p,0,'LineWigth',16,'StarSize',24);
end
makepretty


%% Corellation FR
% get data
for icell=1:length(FR)
    if ~isempty(FR{icell}{1})
        FR_hab(icell) = FR{icell}{1};
    else
        FR_hab(icell) = 0;
    end
    if ~isempty(stats{icell}{4}) && ~isempty(stats{icell}{4}.peak)
        FRpeak_hab(icell) = max(stats{icell}{1}.peak);
    else
        FRpeak_hab(icell) = nan;
    end

    if ~isempty(FR{icell}{4})
        FR_cond(icell) = FR{icell}{4};
    else
        FR_cond(icell) = 0;
    end
    if ~isempty(stats{icell}{4}) && ~isempty(stats{icell}{4}.peak)
        FRpeak_cond(icell) = max(stats{icell}{4}.peak);
    else
        FRpeak_cond(icell) = nan;
    end
end
fc = figure ('units', 'normalized','outerposition', [0 0 .9 .6]);
subplot(121)
scatter(FR_hab,FR_cond, 100, 'k', 'filled');
[rho,p] = corr(FR_hab', FR_cond', 'Type', 'Spearman', 'rows', 'complete');
l1 = lsline;
l1.Color = 'k';
t1 = text(.75, .8, ['rho = ' num2str(round(rho, 3))], 'sc', 'FontSize', 15);
t2 = text(.75, .7, ['p = ' num2str(round(p, 3))], 'sc', 'FontSize', 15);
xlabel('FR during free exploaration')
ylabel('FR during conditioning')
title('Average Firing rates')
makepretty
subplot(122)
scatter(FRpeak_hab,FRpeak_cond, 100, 'k', 'filled');
[rho,p] = corr(FRpeak_hab', FRpeak_cond', 'Type', 'Spearman', 'rows', 'complete');
l1 = lsline;
l1.Color = 'k';
t1 = text(.75, .8, ['rho = ' num2str(round(rho, 3))], 'sc', 'FontSize', 15);
t2 = text(.75, .7, ['p = ' num2str(round(p, 3))], 'sc', 'FontSize', 15);
xlabel('Peak FR during free exploration')
ylabel('Peak FR during conditioning')
title('Peak Firing rates')
makepretty

foldertochoose = ChooseFolderForFigures_DB('PlaceFieldsStability');
saveas(fc,[foldertochoose filesep 'FRcorrelations.fig']);
saveFigure(fc, 'FRcorrelations', foldertochoose);


fc1 = figure ('units', 'normalized','outerposition', [0 0 1 1]);
subplot(121)
b = MakeBoxPlot_DB({FR_hab, FR_cond}, cols, [1 2], xticks, 0);
ylabel('Average firing rate');
title('Average firing rate');
for iplot = 1:length(b)
    b{iplot}.boxAlpha = .45;
end
p = DoWilcoxonOnArray({FR_hab, FR_cond}, {[1 2]});
if p < 0.05
    sigstar_DB({[1 2]},p,0,'LineWigth',16,'StarSize',24);
end
makepretty
subplot(122)
b = MakeBoxPlot_DB({FRpeak_hab, FRpeak_cond}, cols, [1 2], xticks, 0);
ylabel('Peak firing rate');
title('Peak firing rate');
for iplot = 1:length(b)
    b{iplot}.boxAlpha = .45;
end
p = DoWilcoxonOnArray({FRpeak_hab, FRpeak_cond}, {[1 2]});
if p < 0.05
    sigstar_DB({[1 2]},p,0,'LineWigth',16,'StarSize',24);
end
makepretty

foldertochoose = ChooseFolderForFigures_DB('PlaceFieldsStability');
saveas(fc1,[foldertochoose filesep 'FRcorrelations_mean.fig']);
saveFigure(fc1, 'FRcorrelations_mean', foldertochoose);


fc2 = figure ('units', 'normalized','outerposition', [0 0 .9 .6]);
subplot(121)
scatter((FR_hab-FR_cond), Stability{2}, 100, 'k', 'filled');
[rho,p] = corr((FR_hab-FR_cond)',  Stability{2}', 'Type', 'Spearman', 'rows', 'complete');
l1 = lsline;
l1.Color = 'k';
t1 = text(.75, .8, ['rho = ' num2str(round(rho, 3))], 'sc', 'FontSize', 15);
t2 = text(.75, .7, ['p = ' num2str(round(p, 3))], 'sc', 'FontSize', 15);
xlabel('Difference mean FR')
ylabel('Explo-Cond stability')
% title('Average Firing rates')
makepretty
subplot(122)
scatter((FRpeak_hab-FRpeak_cond),Stability{2}, 100, 'k', 'filled');
[rho,p] = corr((FRpeak_hab-FRpeak_cond)', Stability{2}', 'Type', 'Spearman', 'rows', 'complete');
l1 = lsline;
l1.Color = 'k';
t1 = text(.75, .8, ['rho = ' num2str(round(rho, 3))], 'sc', 'FontSize', 15);
t2 = text(.75, .7, ['p = ' num2str(round(p, 3))], 'sc', 'FontSize', 15);
xlabel('Difference peak FR')
ylabel('Explo-Cond stability')
% title('Peak Firing rates')
makepretty

fc3 = figure ('units', 'normalized','outerposition', [0 0 1 .6]);
subplot(131)
scatter(FR_base, Stability{2}, 100, 'k', 'filled');
[rho,p] = corr(FR_base',  Stability{2}', 'Type', 'Spearman', 'rows', 'complete');
l1 = lsline;
l1.Color = 'k';
t1 = text(.75, .8, ['rho = ' num2str(round(rho, 3))], 'sc', 'FontSize', 15);
t2 = text(.75, .7, ['p = ' num2str(round(p, 3))], 'sc', 'FontSize', 15);
xlabel('Full day (base) mean FR')
ylabel('Explo-Cond stability')
makepretty
subplot(132)
scatter((FR_base-FR_hab), Stability{2}, 100, 'k', 'filled');
[rho,p] = corr((FR_base-FR_hab)',  Stability{2}', 'Type', 'Spearman', 'rows', 'complete');
l1 = lsline;
l1.Color = 'k';
t1 = text(.2, .3, ['rho = ' num2str(round(rho, 3))], 'sc', 'FontSize', 15);
t2 = text(.2, .2, ['p = ' num2str(round(p, 3))], 'sc', 'FontSize', 15);
xlabel('Difference Base-Explo')
ylabel('Explo-Cond stability')
makepretty
subplot(133)
scatter((FR_base-FR_cond), Stability{2}, 100, 'k', 'filled');
[rho,p] = corr((FR_base-FR_cond)',  Stability{2}', 'Type', 'Spearman', 'rows', 'complete');
l1 = lsline;
l1.Color = 'k';
t1 = text(.2, .3, ['rho = ' num2str(round(rho, 3))], 'sc', 'FontSize', 15);
t2 = text(.2, .2, ['p = ' num2str(round(p, 3))], 'sc', 'FontSize', 15);
xlabel('Difference Base-Cond')
ylabel('Explo-Cond stability')
makepretty

%% Calculate spatial correlations
cells = [1 6];
for icell=1:length(cells)
    [co_intra{icell},lags] = xcorr(map{cells(icell)}{1}.rate(:),map{cells(icell)}{4}.rate(:));
    co_inter{icell} = xcorr(map{cells(icell)}{2}.rate(:),map{cells(icell)}{3}.rate(:));
end
fcc = figure ('units', 'normalized','outerposition', [0 0 .8 .6]);
for icell = 1:length(cells)
    subplot(1,2,icell)
    p1 = plot(lags, co_inter{icell});
    hold on
    p2 = plot(lags, co_intra{icell});
    xlabel('Lags')
    ylabel('Cross-correlation')
    title(['Explotab=' num2str(round(Stability{1}(cells(icell)),2))...
        ',Explo-Cond Stab=' num2str(round(Stability{2}(cells(icell)),2))])
    legend([p1 p2],  {'Within explo', 'Explo-Cond'})
    makepretty
end

% Calculate zero lag for all
for icell = 1:length(map)
    if ~isempty(map{icell}{4})
        CO_intra(icell) = xcorr(map{icell}{2}.rate(:),map{icell}{3}.rate(:), 0, 'coeff');
        CO_inter(icell) = xcorr(map{icell}{1}.rate(:),map{icell}{4}.rate(:), 0, 'coeff');
    else
        CO_intra(icell) = nan;
        CO_inter(icell) = nan;
    end
end

fcc1 = figure ('units', 'normalized','outerposition', [0 0 .45 .5]);
scatter(CO_intra,CO_inter, 100, 'k', 'filled');
[rho,p] = corr(CO_intra', CO_inter', 'Type', 'Spearman', 'rows', 'complete');
l1 = lsline;
l1.Color = 'k';
t1 = text(.2, .8, ['rho = ' num2str(round(rho, 3))], 'sc', 'FontSize', 15);
t2 = text(.2, .7, ['p = ' num2str(round(p, 3))], 'sc', 'FontSize', 15);
xlabel('Zero-lag cross-corr within explo')
ylabel('Zero-lag cross-corr between explo and cond')
title('Zero-lag cross correlation')
makepretty

foldertochoose = ChooseFolderForFigures_DB('PlaceFieldsStability');
saveas(fcc1,[foldertochoose filesep 'CC_spatial.fig']);
saveFigure(fcc1, 'CC_spatial', foldertochoose);


fcc2 = figure ('units', 'normalized','outerposition', [0 0 .4 .8]);
b = MakeBoxPlot_DB({CO_intra, CO_inter}, cols, [1 2], xticks, 0);
ylabel('Zero-lag cross-corr');
title('Zero-lag cross-correlation');
for iplot = 1:length(b)
    b{iplot}.boxAlpha = .45;
end
p = DoWilcoxonOnArray({CO_intra, CO_inter}, {[1 2]});
if p < 0.05
    sigstar_DB({[1 2]},p,0,'LineWigth',16,'StarSize',24);
end
makepretty

foldertochoose = ChooseFolderForFigures_DB('PlaceFieldsStability');
saveas(fcc2,[foldertochoose filesep 'CC_spatial_mean.fig']);
saveFigure(fcc2, 'CC_spatial_mean', foldertochoose);

%% Remove bad ExploStability
id_badstab = find(Stability{1}<0.5);
cols = {[.2 .2 .2], [.9 0 0]};
xticks = {'Explo', 'Explo-Cond'};
fs = figure ('units', 'normalized','outerposition', [0 0 .4 1]);
Stability_wo = Stability;
Stability_wo{1}(id_badstab) = [];
Stability_wo{2}(id_badstab) = [];
b = MakeBoxPlot_DB(Stability_wo, cols, [1 2], xticks, 0);
ylabel('Rate map correlation');
title(['All place cells,N=' num2str(length(map))])
for iplot = 1:length(b)
    b{iplot}.boxAlpha = .45;
end
p = DoWilcoxonOnArray(Stability_wo, {[1 2]});
if p < 0.05
    sigstar_DB({[1 2]},p,0,'LineWigth',16,'StarSize',24);
end
        

%% Plot good examples

goodcells = idranked([6:8 11]);
badcells = idranked([75 76 89 96]);


% Good place cells
fe1 = figure ('units', 'normalized','outerposition', [0 0 1 1]);
cnt=1;
for icell = goodcells
    subplot(4,5,(cnt-1)*5+1)
    imagesc(map{icell}{1}.rate);
    hold on
    plot(mazeMap(:,1),mazeMap(:,2),'w','LineWidth',1.5)
    plot(ShockZoneMap(:,1),ShockZoneMap(:,2),'r','LineWidth',1.5)
    set(gca,'XTickLabel',{},'YTickLabel',{});
    axis xy
    colormap jet
    
    subplot(4,5,(cnt-1)*5+2)
    imagesc(map{icell}{4}.rate);
    hold on
    plot(mazeMap(:,1),mazeMap(:,2),'w','LineWidth',1.5)
    plot(ShockZoneMap(:,1),ShockZoneMap(:,2),'r','LineWidth',1.5)
    set(gca,'XTickLabel',{},'YTickLabel',{});
    axis xy
    colormap jet


    subplot(4,5,(cnt-1)*5+4)
    imagesc(map{icell}{2}.rate);
    hold on
    plot(mazeMap(:,1),mazeMap(:,2),'w','LineWidth',1.5)
    plot(ShockZoneMap(:,1),ShockZoneMap(:,2),'r','LineWidth',1.5)
    set(gca,'XTickLabel',{},'YTickLabel',{});
    axis xy
    colormap jet
    
    subplot(4,5,(cnt-1)*5+5)
    imagesc(map{icell}{3}.rate);
    hold on
    plot(mazeMap(:,1),mazeMap(:,2),'w','LineWidth',1.5)
    plot(ShockZoneMap(:,1),ShockZoneMap(:,2),'r','LineWidth',1.5)
    set(gca,'XTickLabel',{},'YTickLabel',{});
    axis xy
    colormap jet
    cnt=cnt+1;
end

foldertochoose = ChooseFolderForFigures_DB('PlaceFieldsStability');
saveas(fe1,[foldertochoose filesep 'disser' filesep 'Good_rate.fig']);
saveFigure(fe1, 'Good_rate', [foldertochoose filesep 'disser']);

% spikes traj
fe2 = figure ('units', 'normalized','outerposition', [0 0 1 1]);
cnt=1;
for icell = goodcells
%     subplot(4,5,(cnt-1)*5+1)
    subplot(121)
    plot(Data(Restrict(Session{traj_id(icell)}.AlignedXtsd, Epochs{traj_id(icell)}{1})),...
        Data(Restrict(Session{traj_id(icell)}.AlignedYtsd, Epochs{traj_id(icell)}{1})),...
        'Color',[0.8 0.8 0.8], 'LineWidth', 1.5);
    hold on
    if ~isempty(px{icell}{1})
        plot(px{icell}{1},py{icell}{1},'r.', 'MarkerSize', 80);
    end
    plot(mazeSpikes(:,1),mazeSpikes(:,2),'k','LineWidth',1.5)
    plot(shockZoneSpikes(:,1),shockZoneSpikes(:,2),'r','LineWidth', 1.5)
    set(gca,'XTickLabel',{},'YTickLabel',{});
    xlim([-.1 1.1]);
    ylim([-.1 1.1]);
    
%     subplot(4,5,(cnt-1)*5+2)
    subplot(122)
    plot(Data(Restrict(Session{traj_id(icell)}.AlignedXtsd, Epochs{traj_id(icell)}{4})),...
        Data(Restrict(Session{traj_id(icell)}.AlignedYtsd, Epochs{traj_id(icell)}{4})),...
        'Color',[0.8 0.8 0.8], 'LineWidth', 1.5);
    hold on
    if ~isempty(px{icell}{4})
        plot(px{icell}{4},py{icell}{4},'r.', 'MarkerSize', 80);
    end
    plot(mazeSpikes(:,1),mazeSpikes(:,2),'k','LineWidth',1.5)
    plot(shockZoneSpikes(:,1),shockZoneSpikes(:,2),'r','LineWidth', 1.5)
    set(gca,'XTickLabel',{},'YTickLabel',{});
    xlim([-.1 1.1]);
    ylim([-.1 1.1]);
    

    subplot(4,5,(cnt-1)*5+4)
    plot(Data(Restrict(Session{traj_id(icell)}.AlignedXtsd, Epochs{traj_id(icell)}{2})),...
        Data(Restrict(Session{traj_id(icell)}.AlignedYtsd, Epochs{traj_id(icell)}{2})),...
        'Color',[0.8 0.8 0.8]);
    hold on
    if ~isempty(px{icell}{2})
        plot(px{icell}{2},py{icell}{2},'r.');
    end
    plot(mazeSpikes(:,1),mazeSpikes(:,2),'k','LineWidth',1.5)
    plot(shockZoneSpikes(:,1),shockZoneSpikes(:,2),'r','LineWidth', 1.5)
    set(gca,'XTickLabel',{},'YTickLabel',{});
    xlim([-.1 1.1]);
    ylim([-.1 1.1]);
    
    
    subplot(4,5,(cnt-1)*5+5)
    plot(Data(Restrict(Session{traj_id(icell)}.AlignedXtsd, Epochs{traj_id(icell)}{3})),...
        Data(Restrict(Session{traj_id(icell)}.AlignedYtsd, Epochs{traj_id(icell)}{3})),...
        'Color',[0.8 0.8 0.8]);
    hold on
    if ~isempty(px{icell}{3})
        plot(px{icell}{3},py{icell}{3},'r.');
    end
    plot(mazeSpikes(:,1),mazeSpikes(:,2),'k','LineWidth',1.5)
    plot(shockZoneSpikes(:,1),shockZoneSpikes(:,2),'r','LineWidth', 1.5)
    set(gca,'XTickLabel',{},'YTickLabel',{});
    xlim([-.1 1.1]);
    ylim([-.1 1.1]);

    cnt=cnt+1;
end

foldertochoose = ChooseFolderForFigures_DB('PlaceFieldsStability');
saveas(fe2,[foldertochoose filesep 'disser' filesep 'Good_spikes.fig']);
saveFigure(fe2, 'Good_spikes', [foldertochoose filesep 'disser']);


% Bad place cells
fe3 = figure ('units', 'normalized','outerposition', [0 0 1 1]);
cnt=1;
for icell = badcells
    subplot(4,5,(cnt-1)*5+1)
    imagesc(map{icell}{1}.rate);
    hold on
    plot(mazeMap(:,1),mazeMap(:,2),'w','LineWidth',1.5)
    plot(ShockZoneMap(:,1),ShockZoneMap(:,2),'r','LineWidth',1.5)
    set(gca,'XTickLabel',{},'YTickLabel',{});
    axis xy
    colormap jet
    
    subplot(4,5,(cnt-1)*5+2)
    imagesc(map{icell}{4}.rate);
    hold on
    plot(mazeMap(:,1),mazeMap(:,2),'w','LineWidth',1.5)
    plot(ShockZoneMap(:,1),ShockZoneMap(:,2),'r','LineWidth',1.5)
    set(gca,'XTickLabel',{},'YTickLabel',{});
    axis xy
    colormap jet


    subplot(4,5,(cnt-1)*5+4)
    imagesc(map{icell}{2}.rate);
    hold on
    plot(mazeMap(:,1),mazeMap(:,2),'w','LineWidth',1.5)
    plot(ShockZoneMap(:,1),ShockZoneMap(:,2),'r','LineWidth',1.5)
    set(gca,'XTickLabel',{},'YTickLabel',{});
    axis xy
    colormap jet
    
    subplot(4,5,(cnt-1)*5+5)
    imagesc(map{icell}{3}.rate);
    hold on
    plot(mazeMap(:,1),mazeMap(:,2),'w','LineWidth',1.5)
    plot(ShockZoneMap(:,1),ShockZoneMap(:,2),'r','LineWidth',1.5)
    set(gca,'XTickLabel',{},'YTickLabel',{});
    axis xy
    colormap jet
    cnt=cnt+1;
end

foldertochoose = ChooseFolderForFigures_DB('PlaceFieldsStability');
saveas(fe3,[foldertochoose filesep 'disser' filesep 'Bad_rate.fig']);
saveFigure(fe3, 'Bad_rate', [foldertochoose filesep 'disser']);

% spikes traj
fe4 = figure ('units', 'normalized','outerposition', [0 0 1 1]);
cnt=1;
for icell = badcells
    subplot(4,5,(cnt-1)*5+1)
    plot(Data(Restrict(Session{traj_id(icell)}.AlignedXtsd, Epochs{traj_id(icell)}{1})),...
        Data(Restrict(Session{traj_id(icell)}.AlignedYtsd, Epochs{traj_id(icell)}{1})),...
        'Color',[0.8 0.8 0.8]);
    hold on
    if ~isempty(px{icell}{1})
        plot(px{icell}{1},py{icell}{1},'r.');
    end
    plot(mazeSpikes(:,1),mazeSpikes(:,2),'k','LineWidth',1.5)
    plot(shockZoneSpikes(:,1),shockZoneSpikes(:,2),'r','LineWidth', 1.5)
    set(gca,'XTickLabel',{},'YTickLabel',{});
    xlim([-.1 1.1]);
    ylim([-.1 1.1]);
    
    subplot(4,5,(cnt-1)*5+2)
    plot(Data(Restrict(Session{traj_id(icell)}.AlignedXtsd, Epochs{traj_id(icell)}{4})),...
        Data(Restrict(Session{traj_id(icell)}.AlignedYtsd, Epochs{traj_id(icell)}{4})),...
        'Color',[0.8 0.8 0.8]);
    hold on
    if ~isempty(px{icell}{4})
        plot(px{icell}{4},py{icell}{4},'r.');
    end
    plot(mazeSpikes(:,1),mazeSpikes(:,2),'k','LineWidth',1.5)
    plot(shockZoneSpikes(:,1),shockZoneSpikes(:,2),'r','LineWidth', 1.5)
    set(gca,'XTickLabel',{},'YTickLabel',{});
    xlim([-.1 1.1]);
    ylim([-.1 1.1]);
    

    subplot(4,5,(cnt-1)*5+4)
    plot(Data(Restrict(Session{traj_id(icell)}.AlignedXtsd, Epochs{traj_id(icell)}{2})),...
        Data(Restrict(Session{traj_id(icell)}.AlignedYtsd, Epochs{traj_id(icell)}{2})),...
        'Color',[0.8 0.8 0.8]);
    hold on
    if ~isempty(px{icell}{2})
        plot(px{icell}{2},py{icell}{2},'r.');
    end
    plot(mazeSpikes(:,1),mazeSpikes(:,2),'k','LineWidth',1.5)
    plot(shockZoneSpikes(:,1),shockZoneSpikes(:,2),'r','LineWidth', 1.5)
    set(gca,'XTickLabel',{},'YTickLabel',{});
    xlim([-.1 1.1]);
    ylim([-.1 1.1]);
    
    
    subplot(4,5,(cnt-1)*5+5)
    plot(Data(Restrict(Session{traj_id(icell)}.AlignedXtsd, Epochs{traj_id(icell)}{3})),...
        Data(Restrict(Session{traj_id(icell)}.AlignedYtsd, Epochs{traj_id(icell)}{3})),...
        'Color',[0.8 0.8 0.8]);
    hold on
    if ~isempty(px{icell}{3})
        plot(px{icell}{3},py{icell}{3},'r.');
    end
    plot(mazeSpikes(:,1),mazeSpikes(:,2),'k','LineWidth',1.5)
    plot(shockZoneSpikes(:,1),shockZoneSpikes(:,2),'r','LineWidth', 1.5)
    set(gca,'XTickLabel',{},'YTickLabel',{});
    xlim([-.1 1.1]);
    ylim([-.1 1.1]);

    cnt=cnt+1;
end

foldertochoose = ChooseFolderForFigures_DB('PlaceFieldsStability');
saveas(fe4,[foldertochoose filesep 'disser' filesep 'Bad_spikes.fig']);
saveFigure(fe4, 'Bad_spikes', [foldertochoose filesep 'disser']);


% %% Plot all place cells (8 per figure)
% numplots = floor(length(traj_id)/8);
% ids_full = arrayfun(@(x) 8*x+1:8*x+8, 0:numplots-1, 'UniformOutput', false);
% residual_plot = mod(length(traj_id), 8);
% ids_incompl = length(traj_id)-residual_plot+1:length(traj_id);
% 
% 
% % Rate maps
% for iplot = 1:numplots
%     newline = 0;
%     plots=0;
%     % Rate map figure
%     f{iplot} = figure ('units', 'normalized','outerposition', [0 0 .8 1]);
%     for shift = 1:length(NamesEpochs)
%         annotation('textbox',[.132+.134*(shift-1) .945 .1 .03],'String',NamesEpochs{shift},'FitBoxToText','off', 'FontSize', 14,...
%         'FontWeight', 'bold', 'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%     end
%     
%     celln = 1;
%     for icell=ids_full{iplot}
%         for iepoch = 1:length(NamesEpochs)
%             subplot(8, 6, iepoch+newline*6);
%             if ~isempty(map{icell}{iepoch})
%                 imagesc(map{icell}{iepoch}.rate);
%                 hold on
%                 plot(mazeMap(:,1),mazeMap(:,2),'w','LineWidth',1.5)
%                 plot(ShockZoneMap(:,1),ShockZoneMap(:,2),'r','LineWidth',1.5)
%                 set(gca,'XTickLabel',{},'YTickLabel',{});
%                 axis xy
%                 colormap jet
%                 title(['SpInfo=' num2str(round(stats{icell}{iepoch}.spatialInfo, 2)) ',FR=' num2str(round(FR{icell}{iepoch},2)) ' Hz'])
%             end
%             plots=plots+1;
%             if mod(plots,6) == 0
%                 newline=newline+1;
%             end
%         end
%         annotation('textbox',[.01 .85-.106*(celln-1) .1 .07],'String',...
%             compose([MouseNumber{traj_id(icell)} '\n' 'SG' num2str(ClNum{icell}(1)) ' Cl' num2str(ClNum{icell}(2))]),...
%             'FitBoxToText','off','FontSize', 14,'FontWeight', 'bold',  'LineStyle', 'none',...
%             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%         celln=celln+1;
%     end
% end
% 
% newline=0;
% plots=0;
% fr = figure ('units', 'normalized','outerposition', [0 0 .8 1]);
% for shift = 1:length(NamesEpochs)
%     annotation('textbox',[.132+.134*(shift-1) .945 .1 .03],'String',NamesEpochs{shift},'FitBoxToText','off', 'FontSize', 14,...
%         'FontWeight', 'bold', 'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
% end
% 
% celln = 1;
% for icell=ids_incompl
%     for iepoch = 1:length(NamesEpochs)
%         subplot(residual_plot, 6, iepoch+newline*6);
%         if ~isempty(map{icell}{iepoch})
%             imagesc(map{icell}{iepoch}.rate);
%             hold on
%             plot(mazeMap(:,1),mazeMap(:,2),'w','LineWidth',1.5)
%             plot(ShockZoneMap(:,1),ShockZoneMap(:,2),'r','LineWidth',1.5)
%             set(gca,'XTickLabel',{},'YTickLabel',{});
%             axis xy
%             colormap jet
%             title(['SpInfo=' num2str(round(stats{icell}{iepoch}.spatialInfo, 2)) ',FR=' num2str(round(FR{icell}{iepoch},2)) ' Hz'])
%         end
%         plots=plots+1;
%         if mod(plots,6) == 0
%             newline=newline+1;
%         end
%     end
%     annotation('textbox',[.01 .85-.106*(celln-1) .1 .07],'String',...
%         compose([MouseNumber{traj_id(icell)} '\n' 'SG' num2str(ClNum{icell}(1)) ' Cl' num2str(ClNum{icell}(2))]),...
%         'FitBoxToText','off','FontSize', 14,'FontWeight', 'bold',  'LineStyle', 'none',...
%         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%     celln=celln+1;
% end
% 
% 
% % Spikes on trajectories
% for iplot = 1:numplots
%     newline = 0;
%     plots=0;
%     % Rate map figure
%     f2{iplot} = figure ('units', 'normalized','outerposition', [0 0 .8 1]);
%     for shift = 1:length(NamesEpochs)
%         annotation('textbox',[.132+.134*(shift-1) .945 .1 .03],'String',NamesEpochs{shift},'FitBoxToText','off', 'FontSize', 14,...
%         'FontWeight', 'bold', 'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%     end
%     
%     celln = 1;
%     for icell=ids_full{iplot}
%         for iepoch = 1:length(NamesEpochs)
%             subplot(8, 6, iepoch+newline*6);
%             plot(Data(Restrict(Session{traj_id(icell)}.AlignedXtsd, Epochs{traj_id(icell)}{iepoch})),...
%                 Data(Restrict(Session{traj_id(icell)}.AlignedYtsd, Epochs{traj_id(icell)}{iepoch})),...
%                 'Color',[0.8 0.8 0.8]);
%             hold on
%             if ~isempty(px{icell}{iepoch})
%                 plot(px{icell}{iepoch},py{icell}{iepoch},'r.');
%             end
%             plot(mazeSpikes(:,1),mazeSpikes(:,2),'k','LineWidth',1.5)
%             plot(shockZoneSpikes(:,1),shockZoneSpikes(:,2),'r','LineWidth', 1.5)
%             
%             
%             xlim([-.1 1.1]);
%             ylim([-.1 1.1]);
%             if ~isempty(map{icell}{iepoch})
%                 title(['SpInfo=' num2str(round(stats{icell}{iepoch}.spatialInfo, 2)) ',FR=' num2str(round(FR{icell}{iepoch},2)) ' Hz'])
%             end
%             plots=plots+1;
%             
%             if mod(plots,6) == 0
%                 newline=newline+1;
%             end
%         end
%         annotation('textbox',[.01 .85-.106*(celln-1) .1 .07],'String',...
%             compose([MouseNumber{traj_id(icell)} '\n' 'SG' num2str(ClNum{icell}(1)) ' Cl' num2str(ClNum{icell}(2))]),...
%             'FitBoxToText','off','FontSize', 14,'FontWeight', 'bold',  'LineStyle', 'none',...
%             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%         celln=celln+1;
%     end
% end
% 
% 
% newline = 0;
% plots=0;
% fr2 = figure ('units', 'normalized','outerposition', [0 0 .8 1]);
% for shift = 1:length(NamesEpochs)
%         annotation('textbox',[.132+.134*(shift-1) .945 .1 .03],'String',NamesEpochs{shift},'FitBoxToText','off', 'FontSize', 14,...
%         'FontWeight', 'bold', 'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
% end
% 
% celln = 1;
% for icell=ids_incompl
%     for iepoch = 1:length(NamesEpochs)
%         subplot(residual_plot, 6, iepoch+newline*6);
%         plot(Data(Restrict(Session{traj_id(icell)}.AlignedXtsd, Epochs{traj_id(icell)}{iepoch})),...
%             Data(Restrict(Session{traj_id(icell)}.AlignedYtsd, Epochs{traj_id(icell)}{iepoch})),...
%             'Color',[0.8 0.8 0.8]);
%         hold on
%         if ~isempty(px{icell}{iepoch})
%             plot(px{icell}{iepoch},py{icell}{iepoch},'r.');
%         end
%         plot(mazeSpikes(:,1),mazeSpikes(:,2),'k','LineWidth', 1.5)
%         plot(shockZoneSpikes(:,1),shockZoneSpikes(:,2),'r','LineWidth', 1.5)
%         
%         xlim([-.1 1.1]);
%         ylim([-.1 1.1]);
%         if ~isempty(map{icell}{iepoch})
%             title(['SpInfo=' num2str(round(stats{icell}{iepoch}.spatialInfo, 2)) ',FR=' num2str(round(FR{icell}{iepoch},2)) ' Hz'])
%         end
%         plots=plots+1;
%         
%         if mod(plots,6) == 0
%             newline=newline+1;
%         end
%     end
%     annotation('textbox',[.01 .85-.106*(celln-1) .1 .07],'String',...
%         compose([MouseNumber{traj_id(icell)} '\n' 'SG' num2str(ClNum{icell}(1)) ' Cl' num2str(ClNum{icell}(2))]),...
%         'FitBoxToText','off','FontSize', 14,'FontWeight', 'bold',  'LineStyle', 'none',...
%         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%     celln=celln+1;
% end

    
% %% Bad stability neurons
% idsprepare = idranked(startrank:end);
% numplots = floor(length(idsprepare)/8);
% for iplot = 1:floor(length(idsprepare)/8)
%     ids_full{iplot} = idsprepare((iplot-1)*8+1:(iplot-1)*8+8);
% end
% residual_plot = mod(length(idsprepare), 8);
% ids_incompl = idsprepare(floor(length(idsprepare)/8)*8+1:length(idsprepare));
% 
% foldertochoose = ChooseFolderForFigures_DB('PlaceFieldsStability');
% foldertosave = [foldertochoose filesep 'LowStabEx'];
% 
% % Rate maps
% for iplot = 1:numplots
%     newline = 0;
%     plots=0;
%     % Rate map figure
%     f{iplot} = figure ('units', 'normalized','outerposition', [0 0 .8 1]);
%     for shift = 1:length(NamesEpochs)
%         annotation('textbox',[.132+.134*(shift-1) .945 .1 .03],'String',NamesEpochs{shift},'FitBoxToText','off', 'FontSize', 14,...
%         'FontWeight', 'bold', 'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%     end
%     
%     celln = 1;
%     for icell=ids_full{iplot}
%         for iepoch = 1:length(NamesEpochs)
%             subplot(8, 6, iepoch+newline*6);
%             if ~isempty(map{icell}{iepoch})
%                 imagesc(map{icell}{iepoch}.rate);
%                 hold on
%                 plot(mazeMap(:,1),mazeMap(:,2),'w','LineWidth',1.5)
%                 plot(ShockZoneMap(:,1),ShockZoneMap(:,2),'r','LineWidth',1.5)
%                 set(gca,'XTickLabel',{},'YTickLabel',{});
%                 axis xy
%                 colormap jet
%                 title(['SpInfo=' num2str(round(stats{icell}{iepoch}.spatialInfo, 2)) ',FR=' num2str(round(FR{icell}{iepoch},2)) ' Hz'])
%             end
%             plots=plots+1;
%             if mod(plots,6) == 0
%                 newline=newline+1;
%             end
%         end
%         annotation('textbox',[.01 .85-.106*(celln-1) .1 .12],'String',...
%             compose([MouseNumber{traj_id(icell)} '\n' 'SG' num2str(ClNum{icell}(1)) ' Cl' num2str(ClNum{icell}(2))...
%              '\n Stab=' num2str(round(Stability{2}(icell),2))]),...
%             'FitBoxToText','off','FontSize', 14,'FontWeight', 'bold',  'LineStyle', 'none',...
%             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%         celln=celln+1;
%     end
% 
%     saveas(f{iplot},[foldertosave filesep 'lowstabex_ratemap' num2str(iplot) '.fig']);
%     saveFigure(f{iplot}, ['lowstabex_ratemap' num2str(iplot)], foldertosave);
% end
% 
% newline=0;
% plots=0;
% fr = figure ('units', 'normalized','outerposition', [0 0 .8 1]);
% for shift = 1:length(NamesEpochs)
%     annotation('textbox',[.132+.134*(shift-1) .945 .1 .03],'String',NamesEpochs{shift},'FitBoxToText','off', 'FontSize', 14,...
%         'FontWeight', 'bold', 'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
% end
% 
% celln = 1;
% for icell=ids_incompl
%     for iepoch = 1:length(NamesEpochs)
%         subplot(residual_plot, 6, iepoch+newline*6);
%         if ~isempty(map{icell}{iepoch})
%             imagesc(map{icell}{iepoch}.rate);
%             hold on
%             plot(mazeMap(:,1),mazeMap(:,2),'w','LineWidth',1.5)
%             plot(ShockZoneMap(:,1),ShockZoneMap(:,2),'r','LineWidth',1.5)
%             set(gca,'XTickLabel',{},'YTickLabel',{});
%             axis xy
%             colormap jet
%             title(['SpInfo=' num2str(round(stats{icell}{iepoch}.spatialInfo, 2)) ',FR=' num2str(round(FR{icell}{iepoch},2)) ' Hz'])
%         end
%         plots=plots+1;
%         if mod(plots,6) == 0
%             newline=newline+1;
%         end
%     end
%     annotation('textbox',[.01 .85-.106*(celln-1) .1 .12],'String',...
%         compose([MouseNumber{traj_id(icell)} '\n' 'SG' num2str(ClNum{icell}(1)) ' Cl' num2str(ClNum{icell}(2))...
%              '\n Stab=' num2str(round(Stability{2}(icell),2))]),...
%         'FitBoxToText','off','FontSize', 14,'FontWeight', 'bold',  'LineStyle', 'none',...
%         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%     celln=celln+1;
% end
% saveas(fr,[foldertosave filesep 'lowstabex_ratemap_last' '.fig']);
% saveFigure(fr, 'lowstabex_ratemap_last', foldertosave);
% 
% 
% % Spikes on trajectories
% for iplot = 1:numplots
%     newline = 0;
%     plots=0;
%     % Rate map figure
%     f2{iplot} = figure ('units', 'normalized','outerposition', [0 0 .8 1]);
%     for shift = 1:length(NamesEpochs)
%         annotation('textbox',[.132+.134*(shift-1) .945 .1 .03],'String',NamesEpochs{shift},'FitBoxToText','off', 'FontSize', 14,...
%         'FontWeight', 'bold', 'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%     end
%     
%     celln = 1;
%     for icell=ids_full{iplot}
%         for iepoch = 1:length(NamesEpochs)
%             subplot(8, 6, iepoch+newline*6);
%             plot(Data(Restrict(Session{traj_id(icell)}.AlignedXtsd, Epochs{traj_id(icell)}{iepoch})),...
%                 Data(Restrict(Session{traj_id(icell)}.AlignedYtsd, Epochs{traj_id(icell)}{iepoch})),...
%                 'Color',[0.8 0.8 0.8]);
%             hold on
%             if ~isempty(px{icell}{iepoch})
%                 plot(px{icell}{iepoch},py{icell}{iepoch},'r.');
%             end
%             plot(mazeSpikes(:,1),mazeSpikes(:,2),'k','LineWidth',1.5)
%             plot(shockZoneSpikes(:,1),shockZoneSpikes(:,2),'r','LineWidth', 1.5)
%             
%             
%             xlim([-.1 1.1]);
%             ylim([-.1 1.1]);
%             if ~isempty(map{icell}{iepoch})
%                 title(['SpInfo=' num2str(round(stats{icell}{iepoch}.spatialInfo, 2)) ',FR=' num2str(round(FR{icell}{iepoch},2)) ' Hz'])
%             end
%             plots=plots+1;
%             
%             if mod(plots,6) == 0
%                 newline=newline+1;
%             end
%         end
%         annotation('textbox',[.01 .85-.106*(celln-1) .1 .07],'String',...
%             compose([MouseNumber{traj_id(icell)} '\n' 'SG' num2str(ClNum{icell}(1)) ' Cl' num2str(ClNum{icell}(2))...
%              '\n Stab=' num2str(round(Stability{2}(icell),2))]),...
%             'FitBoxToText','off','FontSize', 14,'FontWeight', 'bold',  'LineStyle', 'none',...
%             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%         celln=celln+1;
%     end
%     saveas(f2{iplot},[foldertosave filesep 'lowstabex_spikes' num2str(iplot) '.fig']);
%     saveFigure(f2{iplot}, ['lowstabex_spikes' num2str(iplot)], foldertosave);
% end
% 
% 
% newline = 0;
% plots=0;
% fr2 = figure ('units', 'normalized','outerposition', [0 0 .8 1]);
% for shift = 1:length(NamesEpochs)
%         annotation('textbox',[.132+.134*(shift-1) .945 .1 .03],'String',NamesEpochs{shift},'FitBoxToText','off', 'FontSize', 14,...
%         'FontWeight', 'bold', 'LineStyle', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
% end
% 
% celln = 1;
% for icell=ids_incompl
%     for iepoch = 1:length(NamesEpochs)
%         subplot(residual_plot, 6, iepoch+newline*6);
%         plot(Data(Restrict(Session{traj_id(icell)}.AlignedXtsd, Epochs{traj_id(icell)}{iepoch})),...
%             Data(Restrict(Session{traj_id(icell)}.AlignedYtsd, Epochs{traj_id(icell)}{iepoch})),...
%             'Color',[0.8 0.8 0.8]);
%         hold on
%         if ~isempty(px{icell}{iepoch})
%             plot(px{icell}{iepoch},py{icell}{iepoch},'r.');
%         end
%         plot(mazeSpikes(:,1),mazeSpikes(:,2),'k','LineWidth', 1.5)
%         plot(shockZoneSpikes(:,1),shockZoneSpikes(:,2),'r','LineWidth', 1.5)
%         
%         xlim([-.1 1.1]);
%         ylim([-.1 1.1]);
%         if ~isempty(map{icell}{iepoch})
%             title(['SpInfo=' num2str(round(stats{icell}{iepoch}.spatialInfo, 2)) ',FR=' num2str(round(FR{icell}{iepoch},2)) ' Hz'])
%         end
%         plots=plots+1;
%         
%         if mod(plots,6) == 0
%             newline=newline+1;
%         end
%     end
%     annotation('textbox',[.01 .85-.106*(celln-1) .1 .07],'String',...
%         compose([MouseNumber{traj_id(icell)} '\n' 'SG' num2str(ClNum{icell}(1)) ' Cl' num2str(ClNum{icell}(2))...
%              '\n Stab=' num2str(round(Stability{2}(icell),2))]),...
%         'FitBoxToText','off','FontSize', 14,'FontWeight', 'bold',  'LineStyle', 'none',...
%         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%     celln=celln+1;
% 
% end
% saveas(fr2,[foldertosave filesep 'lowstabex_spikes_last' '.fig']);
% saveFigure(fr2, 'lowstabex_spikes_last', foldertosave);
% 


end