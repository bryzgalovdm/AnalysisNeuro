function SwapShockSafeZones(dirs)
%
% This function swaps Shock and Safe zone (in case you swap wrongly during the experiment)
%
% INPUT
%
%     dirs      directories with behavResources.mat with wrong zones
% 
%  OUTPUT
%
%     It saves swapped variables in beahvResources.mat
%
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 16/10/2020
% github.com/bryzgalovdm

for idir = 1:length(dirs)
    DoShockSafeZonesSwap(dirs{idir});
end

end



function DoShockSafeZonesSwap(dir)

load([dir 'behavResources.mat']);

Zone_temp = Zone;
ZoneEpoch_temp = ZoneEpoch;
ZoneIndices_temp = ZoneIndices;

%%%% Zone
Zone{1} = Zone_temp{2};
Zone{2} = Zone_temp{1};
Zone{4} = Zone_temp{5};
Zone{5} = Zone_temp{4};
Zone{8} = Zone_temp{9};
Zone{9} = Zone_temp{8};

%%%% ZoneEpoch
ZoneEpoch{1} = ZoneEpoch_temp{2};
ZoneEpoch{2} = ZoneEpoch_temp{1};
ZoneEpoch{4} = ZoneEpoch_temp{5};
ZoneEpoch{5} = ZoneEpoch_temp{4};
ZoneEpoch{8} = ZoneEpoch_temp{9};
ZoneEpoch{9} = ZoneEpoch_temp{8};

%%%% ZoneIndices
ZoneIndices{1} = ZoneIndices_temp{2};
ZoneIndices{2} = ZoneIndices_temp{1};
ZoneIndices{4} = ZoneIndices_temp{5};
ZoneIndices{5} = ZoneIndices_temp{4};
ZoneIndices{8} = ZoneIndices_temp{9};
ZoneIndices{9} = ZoneIndices_temp{8};

if exist('CleanZoneEpoch', 'var')
    CleanZoneEpoch_temp = CleanZoneEpoch;
    CleanZoneIndices_temp = CleanZoneIndices;
    
    
    %%%% CleanZoneEpoch
    CleanZoneEpoch{1} = CleanZoneEpoch_temp{2};
    CleanZoneEpoch{2} = CleanZoneEpoch_temp{1};
    CleanZoneEpoch{4} = CleanZoneEpoch_temp{5};
    CleanZoneEpoch{5} = CleanZoneEpoch_temp{4};
    CleanZoneEpoch{8} = CleanZoneEpoch_temp{9};
    CleanZoneEpoch{9} = CleanZoneEpoch_temp{8};
    
    %%%% CleanZoneIndices
    CleanZoneIndices{1} = CleanZoneIndices_temp{2};
    CleanZoneIndices{2} = CleanZoneIndices_temp{1};
    CleanZoneIndices{4} = CleanZoneIndices_temp{5};
    CleanZoneIndices{5} = CleanZoneIndices_temp{4};
    CleanZoneIndices{8} = CleanZoneIndices_temp{9};
    CleanZoneIndices{9} = CleanZoneIndices_temp{8};
end

save([dir 'behavResources.mat'], 'Zone', 'ZoneEpoch', 'ZoneIndices', '-append');

if exist('CleanZoneEpoch', 'var')
    save([dir 'behavResources.mat'], 'CleanZoneEpoch', 'CleanZoneIndices', '-append');
end


end