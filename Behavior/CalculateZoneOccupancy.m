function occup = CalculateZoneOccupancy(behavResources_slice)
%
% Calculates occupancy in each zone of the Umaze
%
% INPUT
%
%     behavResources_slice      one session structure of concatenated behavioral data
% 
%  OUTPUT
%
%     occup                     occupancies of each zone as ratios
%
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 04/11/2020
% github.com/bryzgalovdm

occup = nan(length(behavResources_slice.CleanZoneIndices), 1);

for izone=1:length(behavResources_slice.Zone)
    occup(izone) = size(behavResources_slice.CleanZoneIndices{izone},1)./...
        size(Data(behavResources_slice.CleanXtsd),1);
end

end