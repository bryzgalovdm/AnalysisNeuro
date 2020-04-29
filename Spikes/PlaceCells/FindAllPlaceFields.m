%%% FindAllPlaceCells
% 
% This script was specifically written to define place cells in the dataset
% of Dima Bryzgalov, MOBS team, France
% 
% If you want to use it, please take time to read this note and take time to
% test the parameters on your dataset
% 
% The script works with the S tsdArray from SpikeData.mat and with CleanAlignedXtsd,
% CleanAlignedYtsd and CleanVtsd from behavResources.mat.
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

%% Parameters

% Do you want to save the PC indices?
sav = false;

% Overwrite
overwrite = true;

% Do you want to save a figure?
savfig = false;

% Paths and names to save
pathfig = '/MOBS_workingON/Dima/Ongoing_results/PlaceField_Final/'; % without dropbox path

% Mice in the analysis
nmouse = [797 798 828 861 882 905 906 911 912 977 994];
% nmouse = [906 912]; % Had PreMazes
% nmouse = [905 911]; % Did not have PreMazes

% Paths retrieved
% Dir = PathForExperimentsERC_Dima('UMazePAG'); 
Dir = PathForExperimentsERC_DimaMAC('UMazePAG');
Dir = RestrictPathForExperiment(Dir,'nMice',nmouse);

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

%% PreAllocation
spikes = cell(1, length(Dir.path));
behav  = cell(1, length(Dir.path));
UMazeEpoch = cell(1, length(Dir.path));
LocomotionEpoch = cell(1, length(Dir.path));
UMazeMovingEpoch = cell(1, length(Dir.path));

%% Load Data
for i=1:length(Dir.path)
    
    % Overwrite
    if ~overwrite
        if isunix
            system(['cp ' Dir.path{i}{1} 'SpikeData.mat ' Dir.path{i}{1} 'SpikeData_old.mat']);
        end
        if ispc
            system(['copy ' Dir.path{i}{1} 'SpikeData.mat ' Dir.path{i}{1} 'SpikeData_old.mat']);
        end
    end
    
    spikes{i} = load([Dir.path{i}{1} '/SpikeData.mat']);
    behav{i} = load([Dir.path{i}{1} '/behavResources.mat'], 'SessionEpoch','CleanAlignedXtsd','CleanAlignedYtsd','CleanVtsd');
    
    if EpochLong
        UMazeEpoch{i} = or(behav{i}.SessionEpoch.Hab,behav{i}.SessionEpoch.TestPre1);
        UMazeEpoch{i} = or(UMazeEpoch{i},behav{i}.SessionEpoch.TestPre2);
        UMazeEpoch{i} = or(UMazeEpoch{i},behav{i}.SessionEpoch.TestPre3);
        UMazeEpoch{i} = or(UMazeEpoch{i},behav{i}.SessionEpoch.TestPre4);
    else
        UMazeEpoch{i} = behav{i}.SessionEpoch.Hab;
    end
    
    LocomotionEpoch{i} = thresholdIntervals(tsd(Range(behav{i}.CleanVtsd),movmedian(Data(behav{i}.CleanVtsd),5)),...
        speed_thresh,'Direction','Above');
    
    UMazeMovingEpoch{i} = and(UMazeEpoch{i}, LocomotionEpoch{i});
        
end

%% Calculate number of cells
totalcells = 0;
for i=1:length(Dir.path)
        totalcells = totalcells + length(spikes{i}.S);
end

%% PreAllocation
map = cell(1, totalcells);
stats_temp = cell(1, totalcells);
stats = cell(1, totalcells);
FR = cell(1, totalcells);
mapout = cell(1, totalcells);
idx = cell(1, length(Dir.path));
PlaceCells_temp = cell(1, length(Dir.path));

%% Calculate rate maps and find place cells

cntAll = 0;
cnt1 = 0;

for i=1:length(Dir.path)
    
    cnt2 = 0;
    lEpoch = sum(End(UMazeMovingEpoch{i}, 's')-Start(UMazeMovingEpoch{i}, 's'));
    PlaceCells_temp{i}.idx = [];
    
    for j=1:length(spikes{i}.S)
        
        cntAll=cntAll+1;
        FRinEpoch = length(Data(Restrict(spikes{i}.S{j},UMazeMovingEpoch{i})))/lEpoch;
        
        if lEpoch > lThresh && FRinEpoch > FRthresh
            try
                [map{cntAll}, ~, stats_temp{cntAll}, ~, ~, FR{cntAll}]=PlaceField_DB(spikes{i}.S{j},behav{i}.CleanAlignedXtsd, behav{i}.CleanAlignedYtsd,...
                    'Epoch', UMazeMovingEpoch{i}, 'SizeMap', sizemap, 'Smoothing', smoothing, 'PlotResults', 0, 'PlotPoisson', 0);
            catch
                stats_temp{cntAll}=[];
            end
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
            if stats_temp{cntAll}.spatialInfo > SIthresh && size_pf > Sthresh_pf
                cnt1 = cnt1+1;
                cnt2=cnt2+1;
                stats{cnt1} = stats_temp{cntAll};
                mapout{cnt1} = map{cntAll};
                idx{cnt1}=[i j];
                PlaceCells_temp{i}.idx(cnt2) = j;
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
for i=1:length(Dir.path)
    PlaceCells_temp{i}.SZ = [];
end

% Record
for i=1:length(overlapCells)
    PlaceCells_temp{idx{overlapCells(i)}(1)}.SZ(end+1)= idx{overlapCells(i)}(2);
end


%% Save the place cell information in the apppropriate folders
if sav
    for i=1:length(Dir.path)
        PlaceCells = PlaceCells_temp{i}
        
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
        
        save([Dir.path{i}{1} 'SpikeData.mat'],'PlaceCells','-append');
    end
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

title([num2str(length(Dir.path)) ' mice, ' num2str(length(stats)) ' PCs, ', num2str(perc_PC) '% of all units found, ' ...
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
        title([Dir.name{idx{i}(1)} ' Cl' num2str(idx{i}(2))])
    else
        imagesc(stats{i}.field)
        axis xy
        hold on
        plot(mazeMap(:,1),mazeMap(:,2),'w','LineWidth',3)
        plot(ShockZoneMap(:,1),ShockZoneMap(:,2),'r','LineWidth',3)
        title([Dir.name{idx{i}(1)} ' Cl' num2str(idx{i}(2))])
    end
    hold off
end
% 
% % Plot all the rate maps
% fa = figure('units', 'normalized', 'outerposition', [0 1 1 1]);
% for i=1:length(stats)
%     subplot(9,10,i)
%     imagesc(mapout{i}.rate)
%     axis xy
%     hold on
%     plot(mazeMap(:,1),mazeMap(:,2),'w','LineWidth',3)
%     plot(ShockZoneMap(:,1),ShockZoneMap(:,2),'r','LineWidth',3)
%     title([Dir.name{idx{i}(1)} ' Cl' num2str(idx{i}(2))])
%     hold off
% end

%% Optional figure - comparison of old and new method

