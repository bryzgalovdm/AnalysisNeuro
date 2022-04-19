function [freezing_shocks, escape_shocks] = SeparateAfterShocks(AlignedXtsd, AlignedYtsd, StimEpoch, timefactor)
%
% This functions separates PAG shocks on the basis of the behavior they
% induce. It checks the position of the mouse <timefactor> sec after
% the shock. If an animal is outside if the shock zone, it is assumed that
% it escaped, otherwise - it did not.
% 
%
%  INPUT
%       
%       AlignedXtsd      Aligned X coordinate (tsd)
%       AlignedYtsd      Aligned Y coordinate (tsd)
%       StimEpoch        interval set with stimulation times
%       timefactor       time after which to compare position with the
%       position of shock in sec (recommended - 3 sec)
% 
% 
%  OUTPUT
%
%       freezing_shocks     serial numbers of shocks after which animal
%       stayed in the zone
%       escape_shocks       serial numbers of shocks after which animal
%       left the zone
% 
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 06/09/2021
% github.com/bryzgalovdm


%% Parameters
bordersShockX = [0 .4];
bordersShockY = [0 .5];

%% Process inputs
% Stimulations
stims = ts(Start(StimEpoch));
afterstims = ts(Start(StimEpoch)+timefactor*1e4);

%% Find locations
% Find locations of the shock
Xshock = unique(Data(Restrict(AlignedXtsd, stims, 'align', 'closest')));
Yshock = unique(Data(Restrict(AlignedYtsd, stims, 'align', 'closest')));
                
% Find location of the mouse <timefactor> sec after the shock
Xaftershock = unique(Data(Restrict(AlignedXtsd, afterstims, 'align', 'closest')));
Yaftershock = unique(Data(Restrict(AlignedYtsd, afterstims, 'align', 'closest')));

%% Separate shocks
cnt1=1;
cnt2=1;
for ishock = 1:length(Xshock)
    if inpolygon(Xshock(ishock),Yshock(ishock),bordersShockX,bordersShockY)
        if inpolygon(Xaftershock(ishock),Yaftershock(ishock),bordersShockX,bordersShockY)
            freezing_shocks(cnt1) = ishock;
            cnt1=cnt2+1;
        else
            escape_shocks(cnt2) = ishock;
            cnt2=cnt2+1;
        end
    end
end
        

end