function id_overlap = FindIDofPC_InZone(S, id_PC, Xtsd, Ytsd, MovingEpoch, Area, varargin)
% 
% Function find indices of place cells in array <PlaceCells> that overlap
% with an area in <Area>. (It could be used to find place cells with place
% fields in particular locations.)
% 
%  INPUT
% 
%         S                         Tsd array with spike trains
%         id_PC                     Indices of place cells in a vector
%         Xtsd                      Tsd with X coordinate
%         Ytsd                      Tsd with Y coordinate
%         MovingEpoch               intervalSet with moving animal in the maze
%                                   to calculate place field
%         Area                      2*2 matrix that delineates area to look for
%                                   place fields into. Must be in units if rate
%                                   maps. First column is borders on X, second
%                                   column is borders on Y
%         SizeMap (optional)        Size of the map to calculate place fields at
%         Smoothing (optional)      Smoothing factor for PlaceField
%         OverlapFactor (optional)  Number of pixels overlapping between <Area> and
%                                   place field that would allow to call place cell
%                                   overlapping with the <Area> :)
% 
%         
%  OUTPUT
%  
%         id_overlap            indices of place cells that overlap with
%                               <Area>  
% 
% 
%  EXAMPLE
% 
%           id_overlap = FindIDofPC_InZone(S, PlaceCells.idx, CleanAlignedXtsd, CleanAlignedYtsd, UMazeMovingEpoch, [7 8; 25 30]);
%           id_overlap = FindIDofPC_InZone(S, PlaceCells.idx,...
%           CleanAlignedXtsd, CleanAlignedYtsd, UMazeMovingEpoch, [7 8; 25 30], 'OverlapFactor', 10);
% 
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 14/12/2020
% github.com/bryzgalovdm

%% Defaults
sizemap = 50; % ---- corresponds to around ~0.8 cm per pixel - maze 39*45 cm
smoothing = 2;
OverlapFactor = 10; % How many pixels in the shockzone for a PC to become SZ-overlapping

%% Parse optional parameters
for i = 1:2:length(varargin)
    if ~ischar(varargin{i})
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help FindIDofPC_InZone">FindIDofPC_InZone</a>'' for details).']);
    end
    switch(lower(varargin{i}))
		case 'sizemap'
			sizemap = varargin{i+1};
            if isdscalar(sizemap,'<0')
				error('Incorrect value for property ''SizeMap'' (type ''help <a href="matlab:help FindIDofPC_InZone">FindIDofPC_InZone</a>'' for details).');
            end
		case 'smoothing'
			smoothing = varargin{i+1};
            if isdscalar(smoothing,'<0')
				error('Incorrect value for property ''Smoothing'' (type ''help <a href="matlab:help FindIDofPC_InZone">FindIDofPC_InZone</a>'' for details).');
            end
        case 'overlapfactor'
			OverlapFactor = varargin{i+1};
            if isdscalar(OverlapFactor,'<0')
				error('Incorrect value for property ''OverlapFactor'' (type ''help <a href="matlab:help FindIDofPC_InZone">FindIDofPC_InZone</a>'' for details).');
            end
    end
end

%% Allocate memory
map = cell(length(id_PC), 1);
stats = cell(length(id_PC), 1);

%% Get place field
for icell = 1:length(id_PC)
    [map{icell}, ~, stats{icell}] = PlaceField_DB(S{id_PC(icell)}, Xtsd, Ytsd,...
     'Epoch', MovingEpoch, 'SizeMap', sizemap, 'Smoothing', smoothing, 'PlotResults', 0, 'PlotPoisson', 0);
end

%% Check whether it overlaps with the area
MaskedMap = zeros(size(map{find(~cellfun('isempty', map), 1)}.rate));
MaskedMap(Area(1,1):Area(2,1),Area(1,2):Area(2,2)) = 1;

overlapCells = zeros(1, length(stats));

% Find ShockZone overlapping place fields
for icell=1:length(stats)
    if ~isempty(stats{icell})
        IsOverlap = FindOverlapArea_PC(stats{icell}, MaskedMap, OverlapFactor);
        if IsOverlap
            overlapCells(icell) = icell;
        end
    end
end
overlapCells = nonzeros(overlapCells);

%% Output figure

id_overlap = overlapCells;


end