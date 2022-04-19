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

% load('E:\ERC_data\M905R\Retracking\TestPost\TestPost4\behavResources_Offline');
% NumEntries = zeros(length(ZoneEpoch), 1);
% for izone = 1:length(ZoneEpoch)
%     if ~isempty(ZoneEpoch{izone})
%         preparedepoch = mergeCloseIntervals(ZoneEpoch{izone}, 2e4); % 2s
%         NumEntries(izone) = length(Start(preparedepoch));
%     end
% end

NumEntries = zeros(length(behavResources_slice.ZoneEpoch), 1);

for izone = 1:length(behavResources_slice.ZoneEpoch)
    if ~isempty(behavResources_slice.ZoneEpoch{izone})
        preparedepoch = mergeCloseIntervals(behavResources_slice.ZoneEpoch{izone}, 2e4); % 2s
        NumEntries(izone) = length(Start(preparedepoch));
    end
end

end