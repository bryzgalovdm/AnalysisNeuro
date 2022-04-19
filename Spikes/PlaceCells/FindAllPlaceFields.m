function FindAllPlaceFields(directory, experiment, varargin)

%%% FindAllPlaceCells
%
% This script was specifically written to define place cells in the dataset
% of Dima Bryzgalov, MOBS team, France
%
% If you want to use it, please take time to read this note and take time to
% test the parameters on your dataset
%
% The script works with the S tsdArray from SpikeData.mat and with AlignedXtsd,
% AlignedYtsd and Vtsd from behavResources.mat.
%
% It will create a structure <PlaceCells> in SpikeData.mat that contains two fields:
% <idx>, with in indices of PCs, and <SZ>, with indices of PCs that overlap with
% the shock zone (see below)
%
% If you don't want to overwrite your existing PlaceCells strucure, set <overwrite>
% to false, and the script will back up your old data into SpikeData_old.mat
%
% It takes each spike train S{j} and restrict it to the time of habituation
% and all pre-tests in which the mouse was moving faster than <speed_thresh>
%
% If you do not want to restrict your spike train to big epoch with interruptions, you can
% set <EpochLong> to false, and S{j} will be restricted only to 15-min habituation session
%
% Only spike trains of length longer than <lThresh> and with firing rate higher than <FRthresh>
% are taken into analysis
%
% Place cells are defined as neurons that:
%  - Have FR higher than <FRthresh> in the movement epoch longer than <lThresh>
%  - Have spatial information bigger than <SIthresh>
%  - Have place fields bigger than <Sthresh_pf>
%     * If algorithm defines two place fields, the script checks whether it least one of them
%             passes <Sthresh_pf>
%           if yes, it defines this neuron as a place cell
%           if no, it checks the euclidian distance between two place fields.
%                     If it is smaller than <Dthresh_pf>, the script checks whether surface sum
%                         of both place fields > <Sthresh_pf>
%
% Maze boundaries are defined manually in <mazeMap>
%
% Place cells that overlap with shock zone are defined if number of overlapping elements
% between place fields and manually defined <ShockZoneMap> is bigger than <OverlapFactor>
%
% ____________________________________________________
%
% Yes, this script is VERY particular case - it's done on purpose
%
% By Dima Bryzgalov, MOBS team, Paris,
% 24/04/2020
% github.com/bryzgalovdm

%% Optional  Parameters

% Do you want to save the PC indices?
sav = true;

% Overwrite
overwrite = false;

% Do you want to save a figure?
savfig = false;

% Whether we do on uninterrupted 15 min (false) or on full Hab + PreTests (true)
EpochLong = true;

% SizeMap
sizemap = 50; % ---- corresponds to around ~0.8 cm per pixel - maze 39*45 cm

% Smoothing of maps
smoothing = 2;

% Firing rate threshold for PC definition (0.25 - is low threshold)
FRthresh = 0.3;

% Epoch length threshold for PC definition
lThresh = 300; % in sec

% Spatial info threshold for PC definition (0.8 - is low threshold)
SIthresh = 0.9;

% Threshold on the size of place field
Sthresh_pf = 2/100; % 2% of the whole surface

% Threshold on the euclidian distance between two fields
Dthresh_pf = 15; % Note that it's euclidian distance - it will correspond to mental distance only in case of linear mazes

% Threshold on speed (im cm/s, epoch with lower are not considered)
speed_thresh = 3;

% How many pixels in the shockzone for a PC to become SZ-overlapping
OverlapFactor = 5;

% Coordinates of the maze and the ShockZone - valid only if sizemap is 50
mazeMap = [6 7; 6 59; 59 59; 59 7; 39 7; 39 42; 24 42; 24 7; 6 7];
ShockZoneMap = [6 7; 6 30; 24 30; 24 7; 6 7];
% mazeMap2 = [24 15; 24 77; 85 77; 85 15; 63 15;  63 58; 46 58; 46 15; 24 15]; - sizemap - ?
% ShockZoneMap2 = [24 15; 24 48; 46 48; 46 15; 24 15];

%% Parameters parsing
for i=1:2:length(varargin)    
    switch(lower(varargin{i})) 
        case 'saveindices'
            sav = varargin{i+1};
            if sav ~= 1 && sav ~= 0
                error('Incorrect value for property ''SaveIndices'' (type ''help FindAllPlaceFields'' for details).');
            end
        case 'overwrite'
            overwrite = varargin{i+1};
            if overwrite ~= 1 && overwrite ~= 0
                error('Incorrect value for property ''OverWrite'' (type ''help FindAllPlaceFields'' for details).');
            end
        case 'epochlong'
            EpochLong = varargin{i+1};
            if EpochLong ~= 1 && EpochLong ~= 0
                error('Incorrect value for property ''EpochLong'' (type ''help FindAllPlaceFields'' for details).');
            end
        case 'sizemap'
            sizemap = varargin{i+1};
            if ~isa(sizemap,'numeric')
                error('Incorrect value for property ''SizeMap'' (type ''help FindAllPlaceFields'' for details).');
            end
        case 'smoothing'
            smoothing = varargin{i+1};
            if ~isa(smoothing,'numeric')
                error('Incorrect value for property ''smoothing'' (type ''help FindAllPlaceFields'' for details).');
            end
        case 'frthresh'
            FRthresh = varargin{i+1};
            if ~isa(FRthresh,'numeric')
                error('Incorrect value for property ''FRthresh'' (type ''help FindAllPlaceFields'' for details).');
            end
        case 'lTtresh'
            lThresh = varargin{i+1};
            if ~isa(lThresh,'numeric')
                error('Incorrect value for property ''lThresh'' (type ''help FindAllPlaceFields'' for details).');
            end
        case 'sithresh'
            SIthresh = varargin{i+1};
            if ~isa(SIthresh,'numeric')
                error('Incorrect value for property ''SIthresh'' (type ''help FindAllPlaceFields'' for details).');
            end
        case 'sthresh_pf'
            Sthresh_pf = varargin{i+1};
            if ~isa(Sthresh_pf,'numeric')
                error('Incorrect value for property ''Sthresh_pf'' (type ''help FindAllPlaceFields'' for details).');
            end
        case 'dthresh_pf'
            Dthresh_pf = varargin{i+1};
            if ~isa(Dthresh_pf,'numeric')
                error('Incorrect value for property ''Dthresh_pf'' (type ''help FindAllPlaceFields'' for details).');
            end
        case 'speed_thresh'
            speed_thresh = varargin{i+1};
            if ~isa(speed_thresh,'numeric')
                error('Incorrect value for property ''speed_thresh'' (type ''help FindAllPlaceFields'' for details).');
            end
        case 'overlapfactor'
            OverlapFactor = varargin{i+1};
            if ~isa(OverlapFactor,'numeric')
                error('Incorrect value for property ''OverlapFactor'' (type ''help FindAllPlaceFields'' for details).');
            end
    end
end

%% Load Data
% Overwrite
if ~overwrite
    if isunix
        system(['cp ' directory filesep 'SpikeData.mat ' directory filesep 'SpikeData_old.mat']);
    end
    if ispc
        system(['copy ' directory filesep 'SpikeData.mat ' directory filesep 'SpikeData_old.mat']);
    end
end

spikes = load([directory filesep 'SpikeData.mat']);
behav = load([directory filesep 'behavResources.mat'], 'SessionEpoch','AlignedXtsd','AlignedYtsd','Vtsd');

numpretests = sum(contains(fields(behav.SessionEpoch), 'TestPre'));
[HabEpoch, ~, UMazeEpoch] = ReturnMnemozyneEpochs(behav.SessionEpoch, 'NumberTests', numpretests);
if EpochLong
    UMazeEpoch = UMazeEpoch;
else
    UMazeEpoch = HabEpoch;
end

LocomotionEpoch = thresholdIntervals(tsd(Range(behav.Vtsd),movmedian(Data(behav.Vtsd),5)),...
    speed_thresh,'Direction','Above');

UMazeMovingEpoch = and(UMazeEpoch, LocomotionEpoch);


%% Calculate number of cells
totalcells = length(spikes.S);

%% PreAllocation
map = cell(1, totalcells);
stats_temp = cell(1, totalcells);
stats = cell(1, totalcells);
FR = cell(1, totalcells);
mapout = cell(1, totalcells);

%% Calculate rate maps and find place cells

cntAll = 0;
cnt1 = 0;
cnt2 = 0;
lEpoch = sum(End(UMazeMovingEpoch, 's')-Start(UMazeMovingEpoch, 's'));
PlaceCells_temp.idx = [];

for j=1:length(spikes.S)

    cntAll=cntAll+1;
    FRinEpoch = length(Data(Restrict(spikes.S{j},UMazeMovingEpoch)))/lEpoch;

    try
        [map{cntAll}, ~, stats_temp{cntAll}, ~, ~, FR{cntAll}]=PlaceField_DB(spikes.S{j},behav.AlignedXtsd, behav.AlignedYtsd,...
            'Epoch', UMazeMovingEpoch, 'SizeMap', sizemap, 'Smoothing', smoothing, 'PlotResults', 0, 'PlotPoisson', 0);
    catch
        stats_temp{cntAll}=[];
    end
    stats_to_save{j} = stats_temp{cntAll};


    if lEpoch > lThresh && FRinEpoch > FRthresh

        if ~isempty(stats_temp{cntAll})
            % Calculate size of the place field for thresholding
            if ~iscell(stats_temp{cntAll}.field)
                size_pf = nnz(stats_temp{cntAll}.field)/(sizemap*sizemap);
            else
                for ss = 1:2
                    size_pf(ss) = nnz(stats_temp{cntAll}.field{ss})/(sizemap*sizemap);
                end
                if size_pf(1) > Sthresh_pf || size_pf(2) > Sthresh_pf
                    size_pf = sum(size_pf);
                else
                    coord = [stats_temp{cntAll}.x(1) stats_temp{cntAll}.y(1); stats_temp{cntAll}.x(2) stats_temp{cntAll}.y(2)];
                    dist = pdist(coord, 'euclidean');
                    if dist < Dthresh_pf
                        size_pf = sum(size_pf);
                    else
                        size_pf = 0;
                    end
                end

            end
        end
        % Write down the place cell structure
        if ~isempty(stats_temp{cntAll}.spatialInfo)
            if stats_temp{cntAll}.spatialInfo > SIthresh && size_pf > Sthresh_pf
                cnt1 = cnt1+1;
                cnt2=cnt2+1;
                stats{cnt1} = stats_temp{cntAll};
                mapout{cnt1} = map{cntAll};
                idx{cnt1}=[1 j];
                PlaceCells_temp.idx(cnt2) = j;     
            end
        end
    end

end
clear map stats_temp FR size_pf dist

% Remove empty cells
stats = stats(~cellfun('isempty',stats));
mapout = mapout(~cellfun('isempty',mapout));
idx = idx(~cellfun('isempty',idx));

% Display how many you found
perc_PC = length(idx)/totalcells*100;

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
fprintf('%0.2f%% of all neurons are place cells\n', perc_PC);
fprintf('%2i total\n', length(idx));
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');

%% Find SZ-overlapping spikes

FakeSZ = zeros (62,62);
FakeSZ (7:25,8:30) = 1;

overlapCells = zeros(1, length(stats));

% Find ShockZone overlapping place fields
cnt3=0;
for i=1:length(stats)
    if iscell(stats{i}.field)
        for k=1:2
            OverlappedFields = FakeSZ & stats{i}.field{k};
            numOverlap = nnz(OverlappedFields);
            if numOverlap > OverlapFactor
                cnt3=cnt3+1;
                overlapCells(cnt3) = i;
                break
            end
        end
    else
        OverlappedFields = FakeSZ & stats{i}.field;
        numOverlap = nnz(OverlappedFields);
        if numOverlap > OverlapFactor
            cnt3=cnt3+1;
            overlapCells(cnt3) = i;
        end
    end
end
overlapCells = nonzeros(overlapCells);

% Allocate
PlaceCells_temp.SZ = [];

if ~strcmp(experiment, 'Novel') 
    % Record
    for i=1:length(overlapCells)
        PlaceCells_temp.SZ(end+1)= idx{overlapCells(i)}(2);
    end
end


%% Save the place cell information in the apppropriate folders
if sav
    PlaceCells = PlaceCells_temp;
    
    PlaceCells.stats = stats_to_save;
    PlaceCells.options.EpochLong = EpochLong;
    PlaceCells.options.sizemap = sizemap;
    PlaceCells.options.smoothing = smoothing;
    PlaceCells.options.FRthresh = FRthresh;
    PlaceCells.options.lThresh = lThresh;
    PlaceCells.options.SIthresh = SIthresh;
    PlaceCells.options.Sthresh_pf = Sthresh_pf;
    PlaceCells.options.Dthresh_pf = Dthresh_pf;
    PlaceCells.options.speed_thresh = speed_thresh;
    PlaceCells.options.OverlapFactor = OverlapFactor;
    
    save([directory filesep 'SpikeData.mat'],'PlaceCells','-append');
end

%% Figure

%Prepare an array to plot
result=zeros(62,62);

for i=1:length(stats)
    if iscell(stats{i}.field)
        for k=1:2
            result = result+stats{i}.field{k};
        end
    else
        result = result+stats{i}.field;
    end
end

fh = figure('units', 'normalized', 'outerposition', [0 1 1 1]);

imagesc(result);
axis xy
colormap jet
hold on
plot(mazeMap(:,1),mazeMap(:,2),'w','LineWidth',4)
plot(ShockZoneMap(:,1),ShockZoneMap(:,2),'r','LineWidth',4)
set(gca,'XTickLabel',{},'YTickLabel',{});

title([num2str(length(stats)) ' PCs, ', num2str(perc_PC) '% of all units found, ' ...
    num2str(length(overlapCells)) ' PCs overlapping with SZ'], 'FontWeight','bold','FontSize',18);

%% Figure

%Prepare an array to plot
result=zeros(62,62);

for i=1:length(stats)
    if iscell(stats{i}.field)
        for k=1:2
            result = result+stats{i}.field{k};
        end
    else
        result = result+stats{i}.field;
    end
end

fh = figure('units', 'normalized', 'outerposition', [0 1 1 1]);

imagesc(result);
axis xy
colormap jet
hold on
plot(mazeMap(:,1),mazeMap(:,2),'w','LineWidth',4)
plot(ShockZoneMap(:,1),ShockZoneMap(:,2),'r','LineWidth',4)
set(gca,'XTickLabel',{},'YTickLabel',{});

title([num2str(length(stats)) ' PCs, ', num2str(perc_PC) '% of all units found, ' ...
    num2str(length(overlapCells)) ' PCs overlapping with SZ'], 'FontWeight','bold','FontSize',18);

if savfig
    saveas(fh,[dropbox pathfig 'AllPlaceFields_Cur.fig']);
    saveFigure(fh,'AllPlaceFields_Cur',[dropbox pathfig]);
end

%% Optional figures - uncomment them if to help you debug
% Fields separately
fi = figure('units', 'normalized', 'outerposition', [0 1 1 1]);
for i=1:length(idx)
    subplot(9,10,i)
    if iscell(stats{i}.field)
        imagesc(stats{i}.field{1}+stats{i}.field{2})
        axis xy
        hold on
        plot(mazeMap(:,1),mazeMap(:,2),'w','LineWidth',3)
        plot(ShockZoneMap(:,1),ShockZoneMap(:,2),'r','LineWidth',3)
        title([' Cl' num2str(idx{i}(2))])
    else
        imagesc(stats{i}.field)
        axis xy
        hold on
        plot(mazeMap(:,1),mazeMap(:,2),'w','LineWidth',3)
        plot(ShockZoneMap(:,1),ShockZoneMap(:,2),'r','LineWidth',3)
        title([' Cl' num2str(idx{i}(2))])
    end
    hold off
end
%
% Plot all the rate maps
fa = figure('units', 'normalized', 'outerposition', [0 1 1 1]);
for i=1:length(stats)
    subplot(9,10,i)
    imagesc(mapout{i}.rate)
    axis xy
    hold on
    plot(mazeMap(:,1),mazeMap(:,2),'w','LineWidth',3)
    plot(ShockZoneMap(:,1),ShockZoneMap(:,2),'r','LineWidth',3)
    title([' Cl' num2str(idx{i}(2))])
    hold off
end

%% Optional figure - comparison of old and new method

if savfig
    saveas(fi,[dropbox pathfig 'AllPlaceFields_Cur.fig']);
    saveFigure(fi,'AllPlaceFields_Cur',[dropbox pathfig]);
end

end