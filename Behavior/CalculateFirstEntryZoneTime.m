function EntryTime = CalculateFirstEntryZoneTime(behavResources_slice, maxtime)
%
% Calculates latency to entry each zone of the Umaze
%
% INPUT
%
%     behavResources_slice      one session structure of concatenated behavioral data
% 
%  OUTPUT
%
%     EntryTime                     entry time in each zone
%
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 04/11/2020
% github.com/bryzgalovdm

EntryTime = nan(length(behavResources_slice.CleanZoneIndices), 1);

for izone = 1:length(behavResources_slice.CleanZoneIndices)
    if isempty(behavResources_slice.CleanZoneIndices{izone})
        EntryTime(izone) = maxtime;
    else
        temp_id = behavResources_slice.CleanZoneIndices{izone}(1);
         EntryTime(izone) = behavResources_slice.CleanPosMat(temp_id(1),1) - ...
            behavResources_slice.CleanPosMat(1,1);
    end
end


end