function numPCs = CountPlaceCells(Dir)

% CountPlaceCells
% 
% Script that quickly counts all recorded place cells
% gets them from SpikeData.PlaceCells field

%% Count
numPCs = 0;
for j=1:length(Dir.path)
    
    % Get the data
    cd(Dir.path{j}{1});
    load('SpikeData.mat','PlaceCells');
    if isfield(PlaceCells,'idx')
        numPCs = numPCs + length(PlaceCells.idx);
        clear PlaceCells
    end
end

end