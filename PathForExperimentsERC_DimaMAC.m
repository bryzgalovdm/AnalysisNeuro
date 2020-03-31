 function Dir=PathForExperimentsERC_DimaMAC(experiment)

% input:
% name of the experiment.
% possible choices:
% 'UMazePAG'
% 'HabBehav', 'TestPreBehav'
% 'Hab' 'PreSleep' 'TestPre' 'Cond' 'PostSleep' 'TestPost' 'TestPost24h'
% 'FindSleep', TestPrePooled', TestPostPooled' 'TestPost24hPooled'
% 'NeutralContextDay1' 'NeutralContextDay2' 'AversiveContextDay1' 'AversiveContextDay2'
% 'Calib', 'FirstExplo'

% output
% Dir = structure containing paths / names / strains / name of the
% experiment (manipe) / correction for amplification (default=1000)
%
%
% example:
% Dir=PathForExperimentsML('EPM');
%
% 	merge two Dir:
% Dir=MergePathForExperiment(Dir1,Dir2);
%lHav
%   restrict Dir to mice or group:
% Dir=RestrictPathForExperiment(Dir,'nMice',[245 246])
% Dir=RestrictPathForExperiment(Dir,'Group',{'OBX','hemiOBX'})
% Dir=RestrictPathForExperiment(Dir,'Group','OBX')
%
% similar functions:
% PathForExperimentFEAR.m
% PathForExperimentsDeltaSleep.m
% PathForExperimentsKB.m PathForExperimentsKBnew.m
% PathForExperimentsML.m

%% strains inputs

% MICEgroups={'Ephys','Behaviour'};
% 
% % Animals that were recorded
% Ephys={'Mouse711' 'Mouse712' 'Mouse714' 'Mouse741' 'Mouse742', 'Mouse753', 'Mouse797', 'Mouse798'};
% 
% % Animals with only behaviour
% Behaviour={'Mouse621' 'Mouse626' 'Mouse627' 'Mouse703' 'Mouse708' 'Mouse743'};

%% Groups

% Concatenated - small zone (starting with 797)
LFP = [1:11]; % Good LFP
Neurons = [1 2 3 5 7 8 9 10 11]; % Good Neurons
ECG = [1 4 5 9 10]; % Good ECG

%% Path
a=0;
I_CA=[];
% 'ERC'
% 'Hab' 'PreSleep' 'TestPre' 'Cond' 'PostSleep' 'TestPost' 'TestPost24h' 'FindSleep' 'UMaze'
%%
if strcmp(experiment,'UMazePAG')
    
    % Mouse797
    a=a+1;Dir.path{a}{1}='/Volumes/Backup Plus/M797/';
    load([Dir.path{a}{1},'ExpeInfo.mat']),Dir.ExpeInfo{a}=ExpeInfo;
    
    % Mouse798
    a=a+1;Dir.path{a}{1}='/Volumes/Backup Plus/M798/';
    load([Dir.path{a}{1},'ExpeInfo.mat']),Dir.ExpeInfo{a}=ExpeInfo;
    
    % Mouse828
    a=a+1;Dir.path{a}{1}='/Volumes/Backup Plus/M828/';
    load([Dir.path{a}{1},'ExpeInfo.mat']),Dir.ExpeInfo{a}=ExpeInfo;
    
    % Mouse861
    a=a+1;Dir.path{a}{1}='/Volumes/Backup Plus/M861/';
    load([Dir.path{a}{1},'ExpeInfo.mat']),Dir.ExpeInfo{a}=ExpeInfo;
    
    % Mouse882
    a=a+1;Dir.path{a}{1}='/Volumes/Backup Plus/M882/';
    load([Dir.path{a}{1},'ExpeInfo.mat']),Dir.ExpeInfo{a}=ExpeInfo;
    
    % Mouse905
    a=a+1;Dir.path{a}{1}='/Volumes/Backup Plus/M905/';
    load([Dir.path{a}{1},'ExpeInfo.mat']),Dir.ExpeInfo{a}=ExpeInfo;
    
    % Mouse906
    a=a+1;Dir.path{a}{1}='/Volumes/Backup Plus/M906/';
    load([Dir.path{a}{1},'ExpeInfo.mat']),Dir.ExpeInfo{a}=ExpeInfo;
    
    % Mouse911
    a=a+1;Dir.path{a}{1}='/Volumes/Backup Plus/M911/';
    load([Dir.path{a}{1},'ExpeInfo.mat']),Dir.ExpeInfo{a}=ExpeInfo;
    
    % Mouse912
    a=a+1;Dir.path{a}{1}='/Volumes/Backup Plus/M912/';
    load([Dir.path{a}{1},'ExpeInfo.mat']),Dir.ExpeInfo{a}=ExpeInfo;
    
    % Mouse977
    a=a+1;Dir.path{a}{1}='/Volumes/Backup Plus/M977/';
    load([Dir.path{a}{1},'ExpeInfo.mat']),Dir.ExpeInfo{a}=ExpeInfo;
    
    % Mouse994
    a=a+1;Dir.path{a}{1}='/Volumes/Backup Plus/M994/20191013/';
%     a=a+1;Dir.path{a}{1}='/Volumes/Backup Plus/M994/20191106/';
load([Dir.path{a}{1},'ExpeInfo.mat']),Dir.ExpeInfo{a}=ExpeInfo;


elseif strcmp(experiment,'FirstExplo')
    
%     % Mouse828
%     a=a+1;Dir.path{a}{1}='/media/nas5/ProjetERC2/Mouse-828/20190301/ExploDay/_Concatenated/';
%     load([Dir.path{a}{1},'ExpeInfo.mat']),Dir.ExpeInfo{a}=ExpeInfo;
%     
%     % Mouse861
%     a=a+1;Dir.path{a}{1}='/media/nas5/ProjetERC2/Mouse-861/20190312/ExploDay/_Concatenated/';
%     load([Dir.path{a}{1},'ExpeInfo.mat']),Dir.ExpeInfo{a}=ExpeInfo;
    
%     % Mouse882
%     a=a+1;Dir.path{a}{1}='/media/nas5/ProjetERC1/M0882/First Exploration/';
%     load([Dir.path{a}{1},'ExpeInfo.mat']),Dir.ExpeInfo{a}=ExpeInfo;
%     
%     % Mouse905
%     a=a+1;Dir.path{a}{1}='/media/DataMOBsRAIDN/ProjetERC2/Mouse-714/27022018/_Concatenated/';
%     load([Dir.path{a}{1},'ExpeInfo.mat']),Dir.ExpeInfo{a}=ExpeInfo;
    
    % Mouse912
    a=a+1;Dir.path{a}{1}='/Volumes/Backup Plus/M912/1/';
    load([Dir.path{a}{1},'ExpeInfo.mat']),Dir.ExpeInfo{a}=ExpeInfo;
    
    % Mouse977
    a=a+1;Dir.path{a}{1}='/Volumes/Backup Plus/M977/1/';
    load([Dir.path{a}{1},'ExpeInfo.mat']),Dir.ExpeInfo{a}=ExpeInfo;
    
    % Mouse979
    a=a+1;Dir.path{a}{1}='/Volumes/Backup Plus/M979/';
    load([Dir.path{a}{1},'ExpeInfo.mat']),Dir.ExpeInfo{a}=ExpeInfo;
    
    % MouseK016
    a=a+1;Dir.path{a}{1}='/Volumes/Backup Plus/M1016/';
    load([Dir.path{a}{1},'ExpeInfo.mat']),Dir.ExpeInfo{a}=ExpeInfo;
    
    
else
    error('Invalid name of experiment')
end


%% Get mice names
for i=1:length(Dir.path)
    Dir.manipe{i}=experiment;
    temp=strfind(Dir.path{i}{1},'M');
    if isempty(temp)
        Dir.name{i}=Dir.path{i}{1}(strfind(Dir.path{i}{1},'Mouse'):strfind(Dir.path{i}{1},'Mouse')+7);
    else
        Dir.name{i}=['Mouse',Dir.path{i}{1}(temp+1:temp+3)];
    end
    %fprintf(Dir.name{i});
end


%% Get Groups

for i=1:length(Dir.path)
    Dir.manipe{i}=experiment;
    if strcmp(Dir.manipe{i},'UMazePAG')
        for j=1:length(LFP)
            Dir.group{1}{LFP(j)} = 'LFP';
        end
        for j=1:length(Neurons)
            Dir.group{2}{Neurons(j)} = 'Neurons';
        end
        for j=1:length(ECG)
            Dir.group{3}{ECG(j)} = 'ECG';
        end
    end
    
    if strcmp(Dir.manipe{i},'Hab') || strcmp(Dir.manipe{i},'TestPrePooled') ||strcmp(Dir.manipe{i},'CondPooled') ||...
            strcmp(Dir.manipe{i},'TestPostPooled')
        for j=1:length(LFP_All)
            Dir.group{1}{LFP_All(j)} = 'LFP';
        end
        for j=1:length(Neurons_All)
            Dir.group{2}{Neurons_All(j)} = 'Neurons';
        end
        for j=1:length(ECG_All)
            Dir.group{3}{ECG_All(j)} = 'ECG';
        end
    end
end

end