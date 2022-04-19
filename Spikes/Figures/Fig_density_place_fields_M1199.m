%%% Density of tuning curves for mouse 1999 reversal 16 april 2021
% 1D


%% Load the data
nmouse = 1199;

% Paths retrieved
Dir = PathForExperimentsERC('Reversal');
Dir = RestrictPathForExperiment(Dir,'nMice',nmouse);

cd(Dir.path{1}{1})
b = load('behavResources.mat', 'SessionEpoch', 'LinearDist');
s = load('SpikeData.mat');

%% Find tuning curve of every cell
% Positions histogram
datLinear = Data(Restrict(b.LinearDist, or(b.SessionEpoch.Hab, b.SessionEpoch.Hab2)));
timeLinear = Range(Restrict(b.LinearDist, or(b.SessionEpoch.Hab, b.SessionEpoch.Hab2)));
tsLinear = ts(timeLinear);
positions = linspace(min(datLinear), max(datLinear), 50);
histPos = hist(datLinear,positions);

% Find position for every spikes
spikespos = cell(length(s.S),1);
for icell = 1:length(s.S)
    tsSpike = ts(Data(Restrict(s.S{icell}, or(b.SessionEpoch.Hab, b.SessionEpoch.Hab2))));
    [~,id] = Restrict(tsLinear, tsSpike);
    spikespos{icell} = [Range(tsSpike) datLinear(id)];
end

% Spike curves across positions
for icell = 1:length(s.S)
    spkPerPos{icell} = hist(spikespos{icell}(:,2),positions);
end

%% Tuning curves
for icell = 1:length(s.S)
    PosTuning{icell} = spkPerPos{icell}./histPos * 15;
end

% Pool tuning curves from place cells
PooledCurves = PosTuning{s.PlaceCells.idx(1)};
for ipc = 2:length(s.PlaceCells.idx)
    PooledCurves = PooledCurves + PosTuning{s.PlaceCells.idx(ipc)};
end

% Pool all tuning curves
PooledCurvesAll = PosTuning{1};
for icell = 2:length(s.S)
    PooledCurvesAll = PooledCurvesAll + PosTuning{icell};
end


%% Plot all tuning curves
% Only place cells
f(1) = figure('units', 'normalized', 'outerposition', [0.1 0.2 0.45 0.65]);
bar(positions,PooledCurves, 'FaceColor', [.1, .41, .515])
% plot(positions,PooledCurves, 'Color', [.1, .41, .515] , 'LineWidth', 2)
xlabel('Linear position')
ylabel('Pooled tuning curves from place cells')
makepretty_DB

% All neurons
f(2) = figure('units', 'normalized', 'outerposition', [0.1 0.2 0.45 0.65]);
bar(positions,PooledCurvesAll, 'FaceColor', [.1, .41, .515])
xlabel('Linear position')
ylabel('Pooled tuning curves from all neurons')
makepretty_DB

% Save
foldertosave = ChooseFolderForFigures_DB('PlaceFieldFinal');
names = {'1D_tuning_PCs', '1D_tuning_All'};
for iplot = 1:length(f)
    saveas(f(iplot),[foldertosave filesep names{iplot} '.fig']);
    saveFigure(f(iplot), names{iplot}, foldertosave);
end
