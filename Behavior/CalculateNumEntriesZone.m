function NumEntries = CalculateNumEntriesZone(behavResources_slice)
%
% Calculates number of entries in each zone of the Umaze
%
% INPUT
%
%     behavResources_slice      one session structure of concatenated behavioral data
% 
%  OUTPUT
%
%     NumEntries                number of entries in each zone
%
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 04/11/2020
% github.com/bryzgalovdm

NumEntries = zeros(length(behavResources_slice.CleanZoneIndices), 1);

for izone = 1:length(behavResources_slice.CleanZoneIndices)
    if ~isempty(behavResources_slice.CleanZoneIndices{izone})
        NumEntries(izone) = length(find(diff(behavResources_slice.CleanZoneIndices{izone})>1))+1;
    end
end


end