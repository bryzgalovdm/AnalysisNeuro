function firingrates = GetFiringRate(S, Epoch)
%
% Calculates firing rates of neurons in a given Epoch
%
% INPUT
%
%     S             tsdaArray with all spike times
%     Epoch         intervalSet to calculate the FRs into
% 
%  OUTPUT
%
%     firingrates   vector with firingrates of all neurons
%
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 19/01/2021
% github.com/bryzgalovdm

firingrates = nan(length(S),1);
lEpoch = sum(End(Epoch, 's') - Start(Epoch, 's'));
for icell = 1:length(S)
    firingrates(icell) = length(S{icell})/lEpoch;
end

end