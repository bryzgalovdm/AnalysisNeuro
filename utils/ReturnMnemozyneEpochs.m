function [HabEpoch, UMazeEpoch, CondEpoch, TaskEpoch, AfterConditioningEpoch] = ReturnMnemozyneEpochs(SessionEpoch, varargin)
% 
% This function returns main types of epochs from ERC-MNEMOZYNE experiments
% 
% INPUT
%   
%   SessionEpoch            structure with epochs for all sessions in the experiment
%   SpeeData (optional)     tsd with speed (needed for filtering out rest periods)
%                           (default = [])
%   SpeedThresh (optional)  Threshold for speed (everything below will be filtered out)
%                           (default = 3 cm/s)
%   NumberTests (optional)  Number of tests in Pre, Cond and Post
%                           (default = 4)
% 
% OUTPUT
% 
%   HabEpoch                Habituation (15  min explo before pretests)
%   UMazeEpoch              Habituation + all pretests
%   CondEpoch               All conditioning epochs
%   TaskEpoch               Habituation + pretests + conditioning
%   AfterConditioningEpoch  All posttests
% 
% By Dima Bryzgalov, MOBS team, Paris,
% 06/07/2020
% github.com/bryzgalovdm

%% Default values of optional arguments
speed_data = [];
speed_thresh = 3;
ntests=4;

%% Optional parameters handling
for i=1:2:length(varargin)    
    switch(lower(varargin{i})) 
        case 'speed'
            speed_data = varargin{i+1};
            if ~isa(speed_data,'tsd')
                error('Incorrect value for property ''SpeedData'' (type ''help PlaceField'' for details).');
            end
        case 'speedthresh'
            speed_thresh = varargin{i+1};
            if ~isa(speed_thresh,'numeric')
                error('Incorrect value for property ''SpeedThresh'' (type ''help PlaceField'' for details).');
            end
        case 'numbertests'
            ntests = varargin{i+1};
            if ~isa(ntests,'numeric')
                error('Incorrect value for property ''NumberTests'' (type ''help PlaceField'' for details).');
            end
    end
end

%% Check for speed requirement
if ~isempty(speed_data)
    smoothspeed  = tsd(Range(speed_data),movmedian(Data(CleanVtsd),5));
    LocomotionEpoch = thresholdIntervals(smoothspeed,speed_thresh,'Direction','Above');
end

%% Create epochs

% HabEpoch
HabEpoch = SessionEpoch.Hab;
if ~isempty(speed_data)
    HabEpoch = and(LocomotionEpoch, HabEpoch);
end

% BaselineExplo Epoch
UMazeEpoch = or(SessionEpoch.Hab,SessionEpoch.TestPre1);
for itest = 2:ntests
    UMazeEpoch = or(UMazeEpoch,SessionEpoch.(['TestPre' num2str(itest)]));
    UMazeEpoch = or(UMazeEpoch,SessionEpoch.(['TestPre' num2str(itest)]));
    UMazeEpoch = or(UMazeEpoch,SessionEpoch.(['TestPre' num2str(itest)]));
end
if ~isempty(speed_data)
    UMazeEpoch = and(LocomotionEpoch, UMazeEpoch);
end

% Conditioning Epoch
CondEpoch = or(SessionEpoch.Cond1,SessionEpoch.Cond2);
if ntests > 2
    for itest = 3:ntests
        CondEpoch = or(CondEpoch,SessionEpoch.(['Cond' num2str(itest)]));
        CondEpoch = or(CondEpoch,SessionEpoch.(['Cond' num2str(itest)]));
    end
end
if ~isempty(speed_data)
    CondEpoch = and(LocomotionEpoch, CondEpoch);
end

% Whole task epoch
TaskEpoch = or(UMazeEpoch,CondEpoch);

% After Conditioning
AfterConditioningEpoch = or(SessionEpoch.TestPost1,SessionEpoch.TestPost2);
if ntests > 2
    for itest = 3:ntests
        AfterConditioningEpoch = or(AfterConditioningEpoch,SessionEpoch.(['TestPost' num2str(itest)]));
        AfterConditioningEpoch = or(AfterConditioningEpoch,SessionEpoch.(['TestPost' num2str(itest)]));
    end
end
if ~isempty(speed_data)
    AfterConditioningEpoch = and(LocomotionEpoch, AfterConditioningEpoch);
end

end