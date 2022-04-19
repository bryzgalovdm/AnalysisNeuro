function [numPCs, numperM] = CountPlaceCells(Dir, varargin)

% CountPlaceCells
% 
% Function that quickly counts all recorded place cells
% gets them from SpikeData.PlaceCells field
% 
% See also FindAllPlaceFields

%% Optional arguments
for i=1:2:length(varargin)
    
    switch(lower(varargin{i}))
        
        case 'verbose'
            verbose = varargin{i+1};
            if ~isa(verbose,'logical')
                error('Incorrect value for property ''Verbose'' (type ''help PlaceCellCrossValidation'' for details).');
            end
    end
end

try
    verbose;
catch
    verbose = true;
end

%% Count
numPCs = 0;
numperM = zeros(1, length(Dir.path));
for j=1:length(Dir.path)
    for k=1:length(Dir.path{i})
        
        % Get the data
        cd(Dir.path{j}{k});
        load('SpikeData.mat','PlaceCells');
        if isfield(PlaceCells,'idx')
            numPCs = numPCs + length(PlaceCells.idx);
            numperM(j) = length(PlaceCells.idx);
            clear PlaceCells
        end
    end
end

if verbose
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
    fprintf('%2i place cells in the analysis\n', numPCs);
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
end

end