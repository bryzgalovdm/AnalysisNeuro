function redectallripples_DB(mice)
%
% This functions redefines and saves the ripples file
%
%  INPUT
%       
%       mice            array with mice numbers that will harvested from 
%                       ERC PathForExperiments
% 
% 
%  OUTPUT
%
%
%       See
%   
%       CreateSleepSignalsSL
% 
% Coded by Dima Bryzgalov and Samuel Laventure, MOBS team, Paris, France
% 08/2020
% github.com/bryzgalovdm

Dir = PathForExperimentsERC_Dima('UMazePAG');
Dir = RestrictPathForExperiment(Dir, 'nMice', mice);

for imouse = 1:length(Dir.path)
    cd(Dir.path{imouse}{1});
    
    CreateSleepSignalsSL('recompute',1, 'scoring','ob','stim',1, ...
    'down',0,'delta',0,'rip',1,'spindle',0, ...
    'ripthresh',[4 6]);
    
end


end