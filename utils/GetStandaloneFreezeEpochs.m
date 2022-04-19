function FreezeEpoch_alone = GetStandaloneFreezeEpochs(FreezeEpoch, StimEpoch, timefactor)

% Get arrays
EndFreeze = End(FreezeEpoch, 's');
stims = Start(StimEpoch, 's');
% Calculate difference betwen end of freezing and closest stim
diffs = nan(length(EndFreeze),1);
for iepoch = 1:length(diffs)
    diffs(iepoch) = min(abs(EndFreeze(iepoch) - stims));
end
% Get freezing epochs that are not co-terminated with stim
id_epoch = find(diffs > timefactor);
% And create new freezing epochs only containing them
st = Start(FreezeEpoch);
en = End(FreezeEpoch);
FreezeEpoch_alone = intervalSet(st(id_epoch), en(id_epoch));

end