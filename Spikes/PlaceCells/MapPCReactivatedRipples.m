function [map_react, allfields, coefs] = MapPCReactivatedRipples(mice, varargin)
% 
% Function builds reactivation maps for place cells active during ripples.
% If a neuron gets active during a ripple, its place field is drawn in
% map_react. map_react is build for:
%   - NREM sleep of PreSleep
%   - Conditoning
%   - NREM sleep of PostSleep
% 
%  INPUT
% 
%         mice                  list of mice in in the analysis
%         ExpType (optional)    experiment type to get the data from: 'MFB'
%                               or 'PAG' (default='PAG') 
%         
%  OUTPUT
%  
%         map_react             map of reactivations
%         allfields             map with all place fields drawn once
%         coefs                 coefficients to multiply allfields if you
%                               want to normalize map_react. Normalization
%                               goes as follows: map_react./(coefs*allfields)
% 
%                               map_react.Pre
%                               map_react.Cond
%                               map_react.Post
% 
% 
%  EXAMPLE
% 
%           [map_react, allfields, coefs] = MapPCReactivatedRipples([797 798 994]);
%           [map_react, allfields, coefs] = MapPCReactivatedRipples([797 798 994], 'ExpType', 'MFB');
% 
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 04/12/2020
% github.com/bryzgalovdm

%% Defaults
Type = 'PAG';

%% Optional arguments management
for i=1:2:length(varargin)
    switch(lower(varargin{i}))
        case 'exptype'
            Type = varargin{i+1};
            if ~isa(Type,'char')
                error('Incorrect value for property ''ExpType'' (type ''help CalcParticipationProbab_ripples'' for details).');
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
    b{imouse} = load([Dir.path{imouse}{1} '/behavResources.mat'], 'SessionEpoch', 'CleanAlignedXtsd', 'CleanAlignedYtsd');
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
AC = cell(length(Dir.path), 1);
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
% Maps
map = cell(length(Dir.path), 1);
stats = cell(length(Dir.path), 1);
FR = cell(length(Dir.path), 1);

%% Get cells and epochs
for imouse = 1:length(Dir.path)
    
    % Sort neuros into place cells and other tentative pyr cells
    if ~isempty(s{imouse}.PlaceCells.idx)
        PC{imouse} = s{imouse}.S(s{imouse}.PlaceCells.idx); % Place cells
    else
        PC{imouse}=[]; % Place cells
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


%% Find all Place fields
for imouse=1:length(Dir.path)
    % Allocate
    map{imouse} = cell(length(PC{imouse}), 1);
    stats{imouse} = cell(length(PC{imouse}), 1);
    FR{imouse} = cell(length(PC{imouse}), 1);
    for icell=1:length(PC{imouse})
        try
            [map{imouse}{icell}, ~, stats{imouse}{icell}, ~, ~, FR{imouse}{icell}] = PlaceField_DB(PC{imouse}{icell},...
                b{imouse}.CleanAlignedXtsd, b{imouse}.CleanAlignedYtsd,...
                'epoch', b{imouse}.SessionEpoch.Hab, 'PlotResults', 0, 'PlotPoisson', 0);
        catch
            map{imouse}{icell} = [];
            stats{imouse}{icell} = [];
            FR{imouse}{icell} = [];
        end
    end
end

%% Create all place fields
%Prepare an array
allfields=zeros(62,62);
for imouse=1:length(Dir.path)
    for icell=1:length(PC{imouse})
        if ~isempty(PC{imouse})
            if isfield(stats{imouse}{icell},'field')
                if iscell(stats{imouse}{icell}.field)
                    for k=1:2
                        allfields = allfields+stats{imouse}{icell}.field{k};
                    end
                else
                    allfields = allfields+stats{imouse}{icell}.field;
                end
            end
        end
    end
end

%% Create spatial reactivation map
map_react.Pre = zeros(62,62);
map_react.Cond = zeros(62,62);
map_react.Post = zeros(62,62);
for imouse=1:length(Dir.path)
    for icell=1:length(PC{imouse})
        if ~isempty(PC{imouse})
            numSpikes_pre = sum(Data(inInterval(RipPreIS{imouse}, PC{imouse}{icell}))); % Count in how many ripples you've got at least a spike
            numSpikes_cond = sum(Data(inInterval(RipCondIS{imouse}, PC{imouse}{icell})));
            numSpikes_post = sum(Data(inInterval(RipPostIS{imouse}, PC{imouse}{icell})));
            if isfield(stats{imouse}{icell},'field')
                if iscell(stats{imouse}{icell}.field)
                    for y = 1:length(stats{imouse}{icell}.field)
                        map_react.Pre = map_react.Pre + numSpikes_pre*stats{imouse}{icell}.field{y};
                        map_react.Cond = map_react.Cond + numSpikes_cond*stats{imouse}{icell}.field{y};
                        map_react.Post = map_react.Post + numSpikes_post*stats{imouse}{icell}.field{y};
                    end
                else
                    map_react.Pre = map_react.Pre + numSpikes_pre*stats{imouse}{icell}.field;
                    map_react.Cond = map_react.Cond + numSpikes_cond*stats{imouse}{icell}.field;
                    map_react.Post = map_react.Post + numSpikes_post*stats{imouse}{icell}.field;
                end
            end
        end
    end
end

%% Calculate coefficients
for imouse = 1:length(Dir.path)
    % Allocate arrays
    ripprob_pc.Pre{imouse} = zeros(length(PC{imouse}), 1);
    ripprob_pc.Cond{imouse} = zeros(length(PC{imouse}), 1);
    ripprob_pc.Post{imouse} = zeros(length(PC{imouse}), 1);
    numRip_pre = 0;
    numRip_cond = 0;
    numRip_post = 0;
    
    for icell=1:length(PC{imouse})
        if ~isempty(PC{imouse})      
            ripprob_pc.Pre{imouse}(icell) = RipProbDB(PC{imouse}{icell}, RipPreIS{imouse});
            numRip_pre = numRip_pre + length(Start(RipPreIS{imouse}));
            ripprob_pc.Cond{imouse}(icell) = RipProbDB(PC{imouse}{icell}, RipCondIS{imouse});
            numRip_cond = numRip_cond + length(Start(RipCondIS{imouse}));
            ripprob_pc.Post{imouse}(icell) = RipProbDB(PC{imouse}{icell}, RipPostIS{imouse});
            numRip_post = numRip_post + length(Start(RipPostIS{imouse}));
        end
    end
    
    prob_pre = flatten_cellarray(ripprob_pc.Pre, 1);
    prob_cond = flatten_cellarray(ripprob_pc.Cond, 1);
    prob_post = flatten_cellarray(ripprob_pc.Post, 1);
    
    coefs.Pre = numRip_pre * nanmean(prob_pre);
    coefs.Cond = numRip_cond * nanmean(prob_cond);
    coefs.Post = numRip_post * nanmean(prob_post);
    
end


end