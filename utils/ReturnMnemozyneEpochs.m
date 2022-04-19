function [HabEpoch, PreTestEpoch, UMazeEpoch, CondEpoch, TaskEpoch, PostTestEpoch] = ReturnMnemozyneEpochs(SessionEpoch, varargin)
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
% EXAMPLE
% 
%   [HabEpoch, UMazeEpoch, CondEpoch, TaskEpoch, AfterConditioningEpoch] = ReturnMnemozyneEpochs(SessionEpoch)
%   [~, UMazeEpoch, CondEpoch, TaskEpoch, AfterConditioningEpoch] = ReturnMnemozyneEpochs(SessionEpoch, 'Speed', Vtsd, 'SpeedThresh', 4);
% 
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
    smoothspeed  = tsd(Range(speed_data),movmedian(Data(speed_data),5));
    LocomotionEpoch = thresholdIntervals(smoothspeed,speed_thresh,'Direction','Above');
end

%% Create epochs

% HabEpoch
if isfield(SessionEpoch, 'Hab2')
    try
        HabEpoch = or(SessionEpoch.Hab1, SessionEpoch.Hab2);
    catch
        HabEpoch = or(SessionEpoch.Hab, SessionEpoch.Hab2);
    end
else
    HabEpoch = SessionEpoch.Hab;
end
if ~isempty(speed_data)
    HabEpoch = and(LocomotionEpoch, HabEpoch);
end

%PreTestEpoch
PreTestEpoch = or(SessionEpoch.TestPre1, SessionEpoch.TestPre2);
if ntests > 2
    for itest = 3:ntests
        try
            PreTestEpoch = or(PreTestEpoch,SessionEpoch.(['TestPre' num2str(itest)]));
        catch
            warning(['No TestPre' num2str(itest) ' found']);
        end
    end
end
if ~isempty(speed_data)
    PreTestEpoch = and(LocomotionEpoch, PreTestEpoch);
end


% BaselineExplo Epoch
UMazeEpoch = or(HabEpoch,PreTestEpoch);
if ~isempty(speed_data)
    UMazeEpoch = and(LocomotionEpoch, UMazeEpoch);
end

% Conditioning Epoch
CondEpoch = or(SessionEpoch.Cond1,SessionEpoch.Cond2);
if ntests > 2
    for itest = 3:ntests
        try 
            CondEpoch = or(CondEpoch,SessionEpoch.(['Cond' num2str(itest)]));
        catch
            warning(['No Cond' num2str(itest) ' found']);
        end
    end
end
if ~isempty(speed_data)
    CondEpoch = and(LocomotionEpoch, CondEpoch);
end

% Whole task epoch
TaskEpoch = or(UMazeEpoch,CondEpoch);

% After Conditioning
if isfield(SessionEpoch, 'TestPost1')
    PostTestEpoch = or(SessionEpoch.TestPost1,SessionEpoch.TestPost2);
    if ntests > 2
        for itest = 3:ntests
            try
                PostTestEpoch = or(PostTestEpoch,SessionEpoch.(['TestPost' num2str(itest)]));
            catch
                warning(['No TestPost' num2str(itest) ' found']);
            end
        end
    end
    if ~isempty(speed_data)
        PostTestEpoch = and(LocomotionEpoch, PostTestEpoch);
    end
else
    PostTestEpoch = [];
end

end