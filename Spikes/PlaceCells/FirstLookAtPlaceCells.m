function SpInfo = FirstLookAtPlaceCells(S, Xtsd, Ytsd, SessionIS, SessionNames)
%
% This function will plot you place fields of all neurons in particular session
% It will create folder PlaceCells in your current folder with 2 figures of each session:
% - Rate maps
% - Spikes over trajectories
%
% INPUT
%
%     S                spike times in tsd array
%     Xtsd             X position in tsd
%     Ytsd             Y position is tsd
%     SessionIS        cell with intervalSets to calculate PF in
%     SessionNames     cell of strings with names of sessions from SessionIS (for plotting)
%
%  OUTPUT
%
%     SpInfo           matrix of spatial information (clusters * sessions)
%
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 18/06/2020
% github.com/bryzgalovdm

%% PreAllocate
map = cell(length(S), length(SessionIS));
stats = cell(length(S), length(SessionIS));
px = cell(length(S), length(SessionIS));
py = cell(length(S), length(SessionIS));
FR = cell(length(S), length(SessionIS));

SpInfo = nan(length(S), length(SessionIS));

%% Do the job
for isess = 1:length(SessionIS)
    for isp = 1:length(S)
        [map{isp, isess}, ~, stats{isp, isess}, px{isp, isess}, py{isp, isess}, FR{isp, isess}] = PlaceField_DB(S{isp},...
            Xtsd, Ytsd, 'Epoch', SessionIS{isess}, 'PlotResults',0, 'PlotPoisson',0);
    end
end

%% Get SpInfo
for isess = 1:length(SessionIS)
    for isp = 1:length(S)
        if ~isempty(stats{isp, isess}) && ~isempty(stats{isp, isess}.spatialInfo)
            SpInfo(isp, isess) = stats{isp, isess}.spatialInfo;
        end
    end
end

%% Plot and save
mkdir([pwd '/PlaceCells']);

% Number of rows and columns
nrows = round(sqrt(length(S)));
ncols = ceil(sqrt(length(S)));

% Loop
for isess = 1:length(SessionIS)
    
    % Rate map figure
    f1 = figure ('units', 'normalized','outerposition', [0 0 1 1]);
    for isp = 1:length(S)
        subplot(nrows, ncols, isp);
        if ~isempty(map{isp, isess})
            imagesc(map{isp, isess}.rate);
            axis xy
            colormap jet
            title(['Cl' num2str(isp) ',FR=' num2str(round(FR{isp, isess},2)) ' Hz'])
        end
    end
    mtit(f1, SessionNames{isess}, 'FontSize', 16, 'xoff', 0, 'yoff', 0.01, 'zoff', 0.04);
    % Save
    saveas(f1, [pwd '/PlaceCells/RateMap_' SessionNames{isess} '.fig']);
    saveFigure(f1, ['RateMap_' SessionNames{isess}], [pwd '/PlaceCells/']);
    
    % Spike position figure
    f2 = figure ('units', 'normalized','outerposition', [0 0 1 1]);
    for isp = 1:length(S)
        subplot(nrows, ncols, isp);
        plot(Data(Restrict(Xtsd, SessionIS{isess})), Data(Restrict(Ytsd, SessionIS{isess})),...
            'Color',[0.8 0.8 0.8]);
        hold on
        if ~isempty(px{isp, isess})
            plot(px{isp, isess},py{isp, isess},'r.');
        end
        
        [xl, yl] = DefineGoodFigLimits_2D(Data(Restrict(Xtsd, SessionIS{isess})),...
            Data(Restrict(Ytsd, SessionIS{isess})));
        xlim(xl);
        ylim(yl);
        title(['Cl' num2str(isp) ',FR=' num2str(round(FR{isp, isess},2)) ' Hz'])
    end
    mtit(f2, SessionNames{isess}, 'FontSize', 16, 'xoff', 0, 'yoff', 0.01, 'zoff', 0.04);
    % Save
    saveas(f2, [pwd '/PlaceCells/SpikesTraj_' SessionNames{isess} '.fig']);
    saveFigure(f2, ['SpikesTraj_' SessionNames{isess}], [pwd '/PlaceCells/']);
    
end

end