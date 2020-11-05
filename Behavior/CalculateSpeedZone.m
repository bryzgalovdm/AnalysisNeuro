function Speed = CalculateSpeedZone(behavResources_slice)
%
% Calculates speed of a mouse in each zone of the Umaze
%
% INPUT
%
%     behavResources_slice      one session structure of concatenated behavioral data
% 
%  OUTPUT
%
%     Speed                speed in each zone
%
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 04/11/2020
% github.com/bryzgalovdm

Speed = zeros(length(behavResources_slice.CleanZoneIndices), 1);

Vtemp=Data(behavResources_slice.CleanVtsd);
for izone = 1:length(behavResources_slice.CleanZoneIndices)
    if ~isempty(behavResources_slice.CleanZoneIndices{izone})
        VZone=Vtemp(behavResources_slice.CleanZoneIndices{izone}(1:end-1),1);
        Speed(izone)=nanmean(VZone,1);
    end
end


end