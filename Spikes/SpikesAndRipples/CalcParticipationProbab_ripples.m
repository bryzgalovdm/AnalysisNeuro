function [ripprob_pc, ripprob_oc] = CalcParticipationProbab_ripples(mice, varargin)
% 
% Function calculates probablility of a given neuron to fire during
% ripples, i.e. number of ripples where a neuron fired at least once
% divided by overall number of ripples for a batch of mice during:
%   - NREM sleep of PreSleep
%   - Conditoning
%   - NREM sleep of PostSleep
% 
%  INPUT
% 
%         mice                  list of mice in in the analysis
%         ExpType (optional)    experiment type to get the data from: 'MFB'
%                               or 'PAG' (default='PAG') 
%         SepArea (optional)    Cell with areas that separate place cells that overlap with them.
%                               Area: 2*2 matrix that delineates area to look for
%                               place fields into. Must be in units if rate
%                               maps. First column is borders on X, second
%                               column is borders on Y
% 
%         
%  OUTPUT
%  
%         ripprob_pc            ripples probability for place cells
%         ripprob_oc            ripples probability for other cells
%                               ripprob_oc.Pre
%                               ripprob_oc.Cond
%                               ripprob_oc.Post
% 
% 
%  EXAMPLE
% 
%           [ripprob_pc, ripprob_oc] = CalcParticipationProbab_ripples([797 798 994]);
%           [ripprob_pc, ripprob_oc] = CalcParticipationProbab_ripples([797 798 994], 'ExpType', 'MFB');
%           [ripprob_pc, ripprob_oc] = CalcParticipationProbab_ripples([797 798 994], 'ExpType', 'MFB', 'SepArea', [7 8; 25 30]);
% 
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% Dec 2020
% github.com/bryzgalovdm

%% Defaults
Type = 'PAG';
sep = [];

%% Optional arguments management
for i=1:2:length(varargin)
    switch(lower(varargin{i}))
        case 'exptype'
            Type = varargin{i+1};
            if ~isa(Type,'char')
                error('Incorrect value for property ''ExpType'' (type ''help CalcParticipationProbab_ripples'' for details).');
            end
        case 'separea'
            sep = varargin{i+1};
            if ~iscell(sep)
                error('Incorrect value for property ''SepArea'' (type ''help CalcParticipationProbab_ripples'' for details).');
            end
            for iarea = 1:length(sep)
                if isdmatrix(sep{iarea}, '<0')
                    error('Incorrect value for property ''SepArea'' (type ''help CalcParticipationProbab_ripples'' for details).');
                end
            end
    end
end

%% Load data
if strcmp(Type, 'PAG')
    if ismac
        [~,comp_name] = system('hostname');
        if contains(comp_name, 'Dmitris-Air') % For debugging
            Dir = PathForExperimentsERC_DimaMAC('UMazePAG'); 
        else
            Dir = PathForExperimentsERC_Dima('UMazePAG');
        end
    else
        Dir = PathForExperimentsERC_Dima('UMazePAG');
    end
    
elseif strcmp(Type, 'MFB')
    Dir = PathForExperimentsERC_SL('StimMFBWake');
end
Dir = RestrictPathForExperiment(Dir,'nMice', mice);

% Allocate arrays
b = cell(length(Dir.path),1); % behavior
r = cell(length(Dir.path),1); % ripples
s = cell(length(Dir.path),1); % spikes
ss = cell(length(Dir.path),1); % sleep scoring

for imouse = 1:length(Dir.path)
    b{imouse} = load([Dir.path{imouse}{1} '/behavResources.mat'], 'SessionEpoch');
    r{imouse} = load([Dir.path{imouse}{1} '/Ripples.mat'], 'ripples');
    s{imouse} = load([Dir.path{imouse}{1} '/SpikeData.mat'], 'S', 'BasicNeuronInfo', 'PlaceCells');
    if strcmp(Type, 'PAG')
        if strcmp(Dir.name{imouse},'Mouse861') || strcmp(Dir.name{imouse},'Mouse906') % bad scoring for 861 and no scoring for 906
            ss{imouse} = load([Dir.path{imouse}{1} 'SleepScoring_Accelero.mat'], 'REMEpoch', 'SWSEpoch', 'Wake', 'Sleep');
        else
            ss{imouse} = load([Dir.path{imouse}{1} 'SleepScoring_OBGamma.mat'], 'REMEpoch', 'SWSEpoch', 'Wake', 'Sleep');
        end
    elseif strcmp(Type, 'MFB')
        ss{imouse} = load([Dir.path{imouse}{1} 'SleepScoring_OBGamma.mat'], 'REMEpoch', 'SWSEpoch', 'Wake', 'Sleep');
    end
end

%% Allocate arrays
% Neurons
PC = cell(length(Dir.path), 1);
OC = cell(length(Dir.path), 1);
% Epochs
ripIS = cell(length(Dir.path), 1);
ripSt = cell(length(Dir.path), 1);
CondEpoch = cell(length(Dir.path), 1);
RipCondIS = cell(length(Dir.path), 1);
SWSPreIS = cell(length(Dir.path), 1);
RipPreIS = cell(length(Dir.path), 1);
SWSPost = cell(length(Dir.path), 1);
RipPostIS = cell(length(Dir.path), 1);
% Prob to fire during ripples
ripprob_pc.Pre = cell(length(Dir.path), 1);
ripprob_pc.Cond = cell(length(Dir.path), 1);
ripprob_pc.Post = cell(length(Dir.path), 1);
ripprob_oc.Pre = cell(length(Dir.path), 1);
ripprob_oc.Cond = cell(length(Dir.path), 1);
ripprob_oc.Post = cell(length(Dir.path), 1);

%% Get cells and epochs
for imouse = 1:length(Dir.path)
    
    % Sort neuros into place cells and other tentative pyr cells
    if ~isempty(s{imouse}.PlaceCells.idx)
        PC{imouse} = s{imouse}.S(s{imouse}.PlaceCells.idx); % Place cells
        OC{imouse} = s{imouse}.S(intersect(setdiff(s{imouse}.BasicNeuronInfo.idx_SUA,s{imouse}.PlaceCells.idx),...
            find(s{imouse}.BasicNeuronInfo.neuroclass>0))); % Other pyramidal cells
    else
        PC{imouse}=[]; % Place cells
        OC{imouse} = []; % Other cells
    end
    
    % Create the most important epochs
    % Ripples
    ripIS{imouse}=intervalSet(r{imouse}.ripples(:,1)*1E4, r{imouse}.ripples(:,3)*1E4);
    ripSt{imouse} = Start(ripIS{imouse});
    % CondEpoch
    CondEpoch{imouse} = or(b{imouse}.SessionEpoch.Cond1,b{imouse}.SessionEpoch.Cond2);
    CondEpoch{imouse} = or(CondEpoch{imouse},b{imouse}.SessionEpoch.Cond3);
    CondEpoch{imouse} = or(CondEpoch{imouse},b{imouse}.SessionEpoch.Cond4);
    
    %%% Create epochs
    SWSPreIS{imouse} = and(ss{imouse}.SWSEpoch,b{imouse}.SessionEpoch.PreSleep);
    RipPreIS{imouse} = and(ripIS{imouse},b{imouse}.SessionEpoch.PreSleep);
    RipCondIS{imouse} = and(ripIS{imouse},CondEpoch{imouse});
    SWSPost{imouse} = and(ss{imouse}.SWSEpoch,b{imouse}.SessionEpoch.PostSleep);
    RipPostIS{imouse} = and(ripIS{imouse},b{imouse}.SessionEpoch.PostSleep);
    
end

%% Calculate probabilities
for imouse = 1:length(Dir.path)
    
    % Allocate arrays
    ripprob_pc.Pre{imouse} = zeros(length(PC{imouse}), 1);
    ripprob_pc.Cond{imouse} = zeros(length(PC{imouse}), 1);
    ripprob_pc.Post{imouse} = zeros(length(PC{imouse}), 1);
    ripprob_oc.Pre{imouse} = zeros(length(OC{imouse}), 1);
    ripprob_oc.Cond{imouse} = zeros(length(OC{imouse}), 1);
    ripprob_oc.Post{imouse} = zeros(length(OC{imouse}), 1);
    
    for icell=1:length(PC{imouse})
        if ~isempty(PC{imouse})      
            ripprob_pc.Pre{imouse}(icell) = RipProbDB(PC{imouse}{icell}, RipPreIS{imouse});
            ripprob_pc.Cond{imouse}(icell) = RipProbDB(PC{imouse}{icell}, RipCondIS{imouse});
            ripprob_pc.Post{imouse}(icell) = RipProbDB(PC{imouse}{icell}, RipPostIS{imouse});
        end
    end
    
    for icell=1:length(OC{imouse})
        if ~isempty(OC{imouse})      
            ripprob_oc.Pre{imouse}(icell) = RipProbDB(OC{imouse}{icell}, RipPreIS{imouse});
            ripprob_oc.Cond{imouse}(icell) = RipProbDB(OC{imouse}{icell}, RipCondIS{imouse});
            ripprob_oc.Post{imouse}(icell) = RipProbDB(OC{imouse}{icell}, RipPostIS{imouse});
        end
    end
    
end

%% Seaprate shcok zone place cells from other place cells (optional)
if ~isempty(sep)
    % Allocate
    temp = cell(length(Dir.path),1);
    id_overlap = cell(length(Dir.path),1);
    id_others = cell(length(Dir.path),1);
    MovingEpoch = cell(length(Dir.path),1);
    
    % Load the data
    for imouse = 1:length(Dir.path)
        temp{imouse} = load([Dir.path{imouse}{1} '/behavResources.mat'], 'CleanAlignedXtsd', 'CleanAlignedYtsd', 'CleanVtsd');
        VtsdSmoothed  = tsd(Range(temp{imouse}.CleanVtsd),movmedian(Data(temp{imouse}.CleanVtsd),5)); % SmoFac = 5
        LocomotionEpoch = thresholdIntervals(VtsdSmoothed,3,'Direction','Above');
        
        UMazeEpoch = or(b{imouse}.SessionEpoch.Hab, b{imouse}.SessionEpoch.TestPre1);
        UMazeEpoch = or(UMazeEpoch, b{imouse}.SessionEpoch.TestPre2);
        UMazeEpoch = or(UMazeEpoch, b{imouse}.SessionEpoch.TestPre3);
        UMazeEpoch = or(UMazeEpoch, b{imouse}.SessionEpoch.TestPre4);
        
        MovingEpoch{imouse} = and(LocomotionEpoch, UMazeEpoch);
    end
    
    % Separate
    for imouse = 1:length(Dir.path)
        if ~isempty(s{imouse}.PlaceCells.idx)
            id_overlap{imouse} = cell(length(sep),1);
            id_others{imouse} = cell(length(sep),1);
            for iarea = 1:length(sep)
                id_overlap{imouse}{iarea} = FindIDofPC_InZone(s{imouse}.S, s{imouse}.PlaceCells.idx,...
                    temp{imouse}.CleanAlignedXtsd, temp{imouse}.CleanAlignedYtsd, MovingEpoch{imouse}, sep{iarea});
                id_overlap{imouse}{iarea} = FindIDofPC_InZone(s{imouse}.S, s{imouse}.PlaceCells.idx,...
                    temp{imouse}.CleanAlignedXtsd, temp{imouse}.CleanAlignedYtsd, MovingEpoch{imouse}, sep{iarea});                
            end
            id_others{imouse} = setdiff(1:length(s{imouse}.PlaceCells.idx), [id_overlap{imouse}{1};...
                    id_overlap{imouse}{2}]);
        end
    end
    
    % Write
    for imouse = 1:length(Dir.path)
        % Pre
        temp = ripprob_pc.Pre{imouse};
        ripprob_pc.Pre{imouse} = [];
        if ~isempty(id_overlap{imouse})
            for iarea = 1:length(sep)
                ripprob_pc.Pre{imouse}{iarea} = temp(id_overlap{imouse}{iarea});
            end
        end
        ripprob_pc.Pre{imouse}{length(sep)+1} = temp(id_others{imouse});
        % Cond
        temp = ripprob_pc.Cond{imouse};
        ripprob_pc.Cond{imouse} = [];
        if ~isempty(id_overlap{imouse})
            for iarea = 1:length(sep)
                ripprob_pc.Cond{imouse}{iarea} = temp(id_overlap{imouse}{iarea});
            end
        end
        ripprob_pc.Cond{imouse}{length(sep)+1} = temp(id_others{imouse});
        % Post
        temp = ripprob_pc.Post{imouse};
        ripprob_pc.Post{imouse} = [];
        if ~isempty(id_overlap{imouse})
            for iarea = 1:length(sep)
                ripprob_pc.Post{imouse}{iarea} = temp(id_overlap{imouse}{iarea});
            end
        end
        ripprob_pc.Post{imouse}{length(sep)+1} = temp(id_others{imouse});
    end
    ripprob_pc.names = cell(length(sep)+1,1);
    for iarea = 1:length(sep)
        ripprob_pc.names{iarea} = ['Area' num2str(iarea)];
    end
    ripprob_pc.names{length(sep)+1} = 'Other place cells';
end
      
end






%% Auxiliary function
function [id_sz, id_am] = GetIDSZneurons(placecells)
    [~, id_sz] = intersect(placecells.idx, placecells.SZ); % shock zone
    id_am = setdiff(1:length(placecells.idx), id_sz); % all maze
end
