%%% FindAllPlaceCells

% TDO: threshold the size of the place field, clean figures code, help


%% Parameters

% Overwrite
overwrite = false;

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

% Spatial info threshold for PC definition
SIthresh = 0.8;

% Firing rate threshold for PC definition
FRthresh = 0.25;

% Epoch length threshold for PC definition
lThresh = 300; % in sec

% SizeMap
sizemap = 50; % ---- corresponds to around ~0.8 cm per pixel - maze 39*45 cm

% Smoothing
smoothing = 2;

% Threshold on speed (im cm/s, epoch with lower are not considered)
speed_thresh = 3;

% How many pixels in the shockzone to become overlapping
OverlapFactor = 5;

% Coordinates of the maze and the ShockZone - valid only if sizemap is 50
mazeMap = [6 7; 6 59; 59 59; 59 7; 39 7; 39 42; 24 42; 24 7; 6 7];
ShockZoneMap = [6 7; 6 30; 24 30; 24 7; 6 7];
% mazeMap2 = [24 15; 24 77; 85 77; 85 15; 63 15;  63 58; 46 58; 46 15; 24 15]; - sizemap - ?
% ShockZoneMap2 = [24 15; 24 48; 46 48; 46 15; 24 15];

% Do you want to save a figure?
sav = true;

% Paths and names to save
pathfig = '/MOBS_workingON/Dima/Ongoing_results/PlaceField_Final/'; % without dropbox path


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
        system(['cp ' Dir.path{i}{1} 'SpikeData.mat ' Dir.path{i}{1} 'SpikeData_old.mat']);
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
                if stats_temp{cntAll}.spatialInfo > SIthresh
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
end

clear map stats_temp FR

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

for i=1:length(Dir.path)
    PlaceCells = PlaceCells_temp{i};
    save([Dir.path{i}{1} 'SpikeData.mat'],'PlaceCells','-append');
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

if sav
    saveas(fh,[dropbox pathfig 'AllPlaceFields_Cur.fig']);
    saveFigure(fh,'AllPlaceFields_Cur',[dropbox pathfig]);
end

%% Optional figures
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

% Plot all the place cells
fa = figure('units', 'normalized', 'outerposition', [0 1 1 1]);
for i=1:length(stats)
    subplot(9,10,i)
    imagesc(mapout{i}.rate)
    axis xy
    hold on
    plot(mazeMap(:,1),mazeMap(:,2),'w','LineWidth',3)
    plot(ShockZoneMap(:,1),ShockZoneMap(:,2),'r','LineWidth',3)
    hold off
end


