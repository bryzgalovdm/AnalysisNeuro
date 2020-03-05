%%%% ExplainedVariance_Kud1999_DB
%
%
% This script calculates explained variance (EV) and reverse EV
% the same manner it was calculated in Kudrimoti et al., 1999
%
% EV = ((R_task,post - R_task,pre*R_post,pre)/sqrt((1-R_task,pre.^2)(1-R_post,pre.^2))).^2;
%
% REV = ((R_task,pre - R_task,post*R_post,pre)/sqrt((1-R_task,post.^2)(1-R_post,pre.^2))).^2;
%
%
% Please go inside the script and check the parameters

%% Parameters
% Mice that go in the analysis
nmouse = [797 798 828 861 882 905 906 911 977 994];
% nmouse = [906 912]; % Had PreMazes
% nmouse = [905 911]; % Did not have PreMazes

% Get paths of each individual mouse
Dir = PathForExperimentsERC_Dima('UMazePAG');
Dir = RestrictPathForExperiment(Dir,'nMice',nmouse);

% Do you want to save the figures?
sav = 0; %%% Does not work now

% Binsize for firing rate histogram (in timestamps!!!)
binsize = 2000; % =100ms (sampling rate = 20000Hz);

% Do you want to split NREM sleep epochs into several consecutive intervals (in min)?
% If not - put []. Note that number of intervals will be floor(40/splitSleep)
% It cannot be more than 40
splitSleep = 20;
if splitSleep > 20
    warning(['There is going to be only one interval: ' num2str(splitSleep) ' min']);
end

%% Allocate memory
if ~isempty(splitSleep)
    if splitSleep <= 20
        % Calculate number of intervals ~~~ !!! 40 min hardcoded !!! ~~~
        numint = floor(40/splitSleep);
        
        PreSleepSWS_Split = cell(1,length(Dir.path)); % interval set
        PostSleepSWS_Split = cell(1,length(Dir.path)); % interval set
        QPRESWS_Split = cell(1,length(Dir.path)); % firing rate histogram
        QPOSTSWS_Split = cell(1,length(Dir.path)); % firing rate histogram
        CorrCo_PostSWS_Split_Task = cell(1,length(Dir.path));
        CorrCo_TaskPreSWS_Split = cell(1,length(Dir.path));
        CorrCo_PostSWSPreSWS_Split = cell(1,length(Dir.path));
        EVSWS_Split = cell(1,numint);
        REVSWS_Split = cell(1,numint);
        
        
        for j = 1:length(length(Dir.path))
            QPRESWS_Split{j} = cell(1,numint);
            QPOSTSWS_Split{j} = cell(1,numint);
            CorrCo_PostSWS_Split_Task{j} = zeros(1,numint);
            CorrCo_TaskPreSWS_Split{j} = zeros(1,numint);
            CorrCo_PostSWSPreSWS_Split{j} = zeros(1,numint);
        end
        for i=1:length(numint)
            EVSWS_Split{i} = zeros(1,length(Dir.path));
            REVSWS_Split{i} = zeros(1,length(Dir.path));
        end
    end
end

% firing rate histograms
QPRE = cell(1,length(Dir.path));
QPRESWS = cell(1,length(Dir.path));
QPREREM = cell(1,length(Dir.path));
QTASK = cell(1,length(Dir.path));
QPOST = cell(1,length(Dir.path));
QPOSTSWS = cell(1,length(Dir.path));
QPOSTREM = cell(1,length(Dir.path));
QPOSTTEST = cell(1,length(Dir.path));

% Structure with correlation matrices
CorrM = cell(1,length(Dir.path));

% Correlation coeficients
CorrCo_PostTask = zeros(1,length(Dir.path));
CorrCo_TaskPre = zeros(1,length(Dir.path));
CorrCo_PostPre = zeros(1,length(Dir.path));

CorrCo_PostSWSTask = zeros(1,length(Dir.path));
CorrCo_PostREMTask = zeros(1,length(Dir.path));
CorrCo_TaskPreSWS = zeros(1,length(Dir.path));
CorrCo_TaskPreREM = zeros(1,length(Dir.path));
CorrCo_PostSWSPreSWS = zeros(1,length(Dir.path));
CorrCo_PostREMPreREM = zeros(1,length(Dir.path));

% Explained variance and reverse explained variance
EV = zeros(1,length(Dir.path));
REV = zeros(1,length(Dir.path));
EVSWS = zeros(1,length(Dir.path));
REVSWS = zeros(1,length(Dir.path));
EVREM = zeros(1,length(Dir.path));
REVREM = zeros(1,length(Dir.path));

%% Load the data
for j=1:length(Dir.path)
    
    cd(Dir.path{j}{1});
    load('SpikeData.mat','S','PlaceCells');
    % If there are less than 2 PCs - don't do
    if isfield(PlaceCells,'idx')
        if length(PlaceCells.idx)>2
           
            load('behavResources.mat','SessionEpoch', 'CleanVtsd');
            load('Ripples.mat','ripples');
            if strcmp(Dir.name{j}, 'Mouse906') || strcmp(Dir.name{j}, 'Mouse977') % Mice with bad OB-based sleep scoring
                load('SleepScoring_Accelero.mat','SWSEpoch','REMEpoch','Sleep'); % Sleep is not used
            else
                load('SleepScoring_OBGamma.mat','SWSEpoch','REMEpoch','Sleep');  % Sleep is not used
            end
            
            %% Create interval sets (epochs)
            
            % Split epochs if necessary
            if ~isempty(splitSleep)
                if splitSleep <= 20
                    PreSleepSWS_Split{j} = SplitIntervals(and(SessionEpoch.PreSleep, SWSEpoch),...
                        splitSleep*60*1e4);
                    PostSleepSWS_Split{j} = SplitIntervals(and(SessionEpoch.PostSleep, SWSEpoch),...
                        splitSleep*60*1e4);
                    
                    % Allocate memory for split intervals 
                    idx_nonexist_PreSWS_Split = cell(1,numint);
                    idx_nonexist_PostSWS_Split = cell(1,numint);
                   
                end
            end
            
            % Create ripples epochs
            ripplesEpoch = intervalSet(ripples(:,2)*1e4-1e4,ripples(:,2)*1e4+1e4);
            PreSleepRipplesEpoch = and(and(SessionEpoch.PreSleep,SWSEpoch),ripplesEpoch);
            PostSleepRipplesEpoch = and(and(SessionEpoch.PostSleep,SWSEpoch),ripplesEpoch);
            
            % BaselineExplo Epoch
            UMazeEpoch = or(SessionEpoch.Hab,SessionEpoch.TestPre1);
            UMazeEpoch = or(UMazeEpoch,SessionEpoch.TestPre2);
            UMazeEpoch = or(UMazeEpoch,SessionEpoch.TestPre3);
            UMazeEpoch = or(UMazeEpoch,SessionEpoch.TestPre4);
            
            % After Conditioning
            AfterConditioningEpoch = or(SessionEpoch.TestPost1,SessionEpoch.TestPost2);
            AfterConditioningEpoch = or(AfterConditioningEpoch,SessionEpoch.TestPost3);
            AfterConditioningEpoch = or(AfterConditioningEpoch,SessionEpoch.TestPost4);
            
            
            % Locomotion threshold
            VtsdSmoothed  = tsd(Range(CleanVtsd),movmedian(Data(CleanVtsd),5));
            LocomotionEpoch = thresholdIntervals(VtsdSmoothed,3,'Direction','Above');
            
            % Get resulting epoch
            UMazeMovingEpoch = and(LocomotionEpoch, UMazeEpoch);
            AfterConditioningMovingEpoch = and(LocomotionEpoch, AfterConditioningEpoch);
            
            %% Create firing rate histograms
            
            % Bin the trains
            Q=MakeQfromS(S,binsize);
            
            % Create epochs for different periods of the day
                        
            % PreSleep full
            QPRE{j}=zscore(full(Data(Restrict(Q,SessionEpoch.PreSleep))));
            
            % PreSleep NREM
            QPRESWS{j}=zscore(full(Data(Restrict(Q,and(SessionEpoch.PreSleep,SWSEpoch))))); % full sleep
%             QPRESWS{j}=zscore(full(Data(Restrict(Q,PreSleepRipplesEpoch)))); % only during nREM ripples

            % PreSleep NREM binned if necessary
            if ~isempty(splitSleep)
                if splitSleep <= 20
                    for i=1:numint
                        QPRESWS_Split{j}{i} = zscore(full(Data(Restrict(Q,PreSleepSWS_Split{j}{i}))));
                    end
                end
            end
            
            % PreSleep REM
            QPREREM{j}=zscore(full(Data(Restrict(Q,and(SessionEpoch.PreSleep,REMEpoch)))));
            
            % Exploration + PreTests (Locomotion only)
            QTASK{j}=zscore(full(Data(Restrict(Q,UMazeMovingEpoch))));
            
            % PostSleep full
            QPOST{j}=full(Data(Restrict(Q,SessionEpoch.PostSleep)));
            
            % PostSleep NREM
            QPOSTSWS{j}=zscore(full(Data(Restrict(Q,and(SessionEpoch.PostSleep,SWSEpoch))))); % full sleep
            %             QPOSTSWS{j}=zscore(full(Data(Restrict(Q,PostSleepRipplesEpoch)))); % only during nREM ripples
            
            % PostSleep NREM binned if necessary
            if ~isempty(splitSleep)
                if splitSleep <= 20
                    for i=1:numint
                        QPOSTSWS_Split{j}{i} = zscore(full(Data(Restrict(Q,PostSleepSWS_Split{j}{i}))));
                    end
                end
            end
            
            % PostSleep REM
            QPOSTREM{j}=zscore(full(Data(Restrict(Q,and(SessionEpoch.PostSleep,REMEpoch)))));
            
            % PostTests (Locomotion only)
            QPOSTTEST{j}=zscore(full(Data(Restrict(Q,AfterConditioningMovingEpoch))));
            
            % Get rid the variable unneccessry in future
            clear Q SWSEpoch REMEpoch SessionEpoch ripples CleanVtsd
            
            %% Calculate the correlation maps and coefficients  
            
            % PreSleep full
            CorrM{j}.Pre = corr(QPRE{j});
            
            % PreSleep NREM
            CorrM{j}.PreSWS = corr(QPRESWS{j});
            
            % PreSleep NREM binned if necessary
            if ~isempty(splitSleep)
                if splitSleep <= 20
                    for i=1:numint
                        CorrM{j}.PreSWS_Split{i}=corr(QPRESWS_Split{j}{i});
                    end
                end
            end
            
            % PreSleep REM
            CorrM{j}.PreREM = corr(QPREREM{j});
            
            % Exploration + PreTests (Locomotion only)
            CorrM{j}.Task_full = corr(QTASK{j});
            CorrM{j}.Task_stages = corr(QTASK{j});
            if ~isempty(splitSleep)
                if splitSleep <= 20
                    CorrM{j}.Task_splitStages = corr(QTASK{j});
                end
            end
            
            % PostSleep full
            CorrM{j}.Post = corr(QPOST{j});
            
            % PreSleep REM
            CorrM{j}.PreREM = corr(QPREREM{j});
          
            % PostSleep NREM            
            CorrM{j}.PostSWS = corr(QPOSTSWS{j});
            
            % PostSleep NREM binned if necessary
            if ~isempty(splitSleep)
                if splitSleep <= 20
                    for i=1:numint
                        CorrM{j}.PostSWS_Split{i}=corr(QPOSTSWS_Split{j}{i});
                    end
                end
            end
            
            % PostSleep REM
            CorrM{j}.PostREM = corr(QPOSTREM{j});
                        
            % PostTests (Locomotion only)
            CorrM{j}.PostTest_full = corr(QPOSTTEST{j});
            CorrM{j}.PostTest_stages = corr(QPOSTTEST{j});
            if ~isempty(splitSleep)
                if splitSleep <= 20
                    CorrM{j}.PostTest_stagesSplit = corr(QPOSTTEST{j});
                end
            end
            
            
            %% Find non-firing neurons in each epoch
            
            idx_nonexist_Pre = find(isnan(CorrM{j}.Pre(:,1)));
            idx_nonexist_PreSWS = find(isnan(CorrM{j}.PreSWS(:,1)));
            
            if ~isempty(splitSleep)
                if splitSleep <= 20
                    for i=1:numint
                        idx_nonexist_PreSWS_Split{i} = find(isnan(CorrM{j}.PreSWS_Split{i}(:,1)));
                    end
                    idx_nonexist_PreSWS_Split_Final = idx_nonexist_PreSWS_Split{1};
                    for i = 2:numint
                        idx_nonexist_PreSWS_Split_Final = [idx_nonexist_PreSWS_Split_Final idx_nonexist_PreSWS_Split{i}];
                    end
                    
                    idx_nonexist_PreSWS_Split_Final = unique(idx_nonexist_PreSWS_Split_Final);
                end
            end
            
            idx_nonexist_PreREM = find(isnan(CorrM{j}.PreREM(:,1)));
            
            idx_nonexist_Task = find(isnan(CorrM{j}.Task_full(:,1)));
            
            idx_nonexist_Post = find(isnan(CorrM{j}.Post(:,1)));
            idx_nonexist_PostSWS = find(isnan(CorrM{j}.PostSWS(:,1)));
            
            if ~isempty(splitSleep)
                if splitSleep <= 20
                    for i=1:numint
                        idx_nonexist_PostSWS_Split{i} = find(isnan(CorrM{j}.PostSWS_Split{i}(:,1)));
                    end
                    idx_nonexist_PostSWS_Split_Final = idx_nonexist_PostSWS_Split{1};
                    for i = 2:numint
                        idx_nonexist_PostSWS_Split_Final = [idx_nonexist_PostSWS_Split_Final idx_nonexist_PostSWS_Split{i}];
                    end
                    
                    idx_nonexist_PostSWS_Split_Final = unique(idx_nonexist_PostSWS_Split_Final);
                end
            end

            idx_nonexist_PostREM = find(isnan(CorrM{j}.PostREM(:,1)));
            
            idx_nonexist_PostTest = find(isnan(CorrM{j}.PostTest_full(:,1)));
            
            % Three sets of non-firing indexes:
            % - Sleep - Run - Sleep - Run (_full)
            % - NREM - REM - Run - NREM - REM - Run (_stages)
            % - NREM_Split - REM_Split - Run - NREM_Split - REM_Split - Run (_splitStages)
            
            idx_toremove_full = unique(([idx_nonexist_Pre; idx_nonexist_Task; idx_nonexist_Post; idx_nonexist_PostTest])');
            idx_toremove_stages = unique(([idx_nonexist_PreSWS; idx_nonexist_PreREM;...
                idx_nonexist_Task; idx_nonexist_PostSWS; idx_nonexist_PostREM; idx_nonexist_PostTest])');
            if ~isempty(splitSleep)
                if splitSleep <= 20
                    idx_toremove_splitStages = unique(([idx_nonexist_PreSWS_Split_Final; ...
                        idx_nonexist_Task; idx_nonexist_PostSWS_Split_Final; idx_nonexist_PostTest])');
                end
            end
            
            % Get rid the variable unneccessry in future
            clear idx_nonexist_Pre idx_nonexist_PreSWS idx_nonexist_PreSWS_Split idx_nonexist_PreSWS_Split_Final
            clear idx_nonexist_PreREM idx_nonexist_Task idx_nonexist_Post idx_nonexist_PostSWS
            clear idx_nonexist_PostSWS_Split idx_nonexist_PostSWS_Split_Final idx_nonexist_PostREM idx_nonexist_PostTest
            
            %% Remove neurons that do not fire in some epochs
            
            %%%%%%%% Intance _full %%%%%%%%
            % PreSleep full
            CorrM{j}.Pre(idx_toremove_full,:) = [];
            CorrM{j}.Pre(:,idx_toremove_full) = [];
            % Exploration + PreTests (Locomotion only)
            CorrM{j}.Task_full(:,idx_toremove_full) = [];
            CorrM{j}.Task_full(idx_toremove_full,:) = [];
            % PostSleep full
            CorrM{j}.Post(:,idx_toremove_full) = [];
            CorrM{j}.Post(idx_toremove_full,:) = [];
            % PostTests (Locomotion only)
            CorrM{j}.PostTest_full(:,idx_toremove_full) = [];
            CorrM{j}.PostTest_full(idx_toremove_full,:) = [];
           
            %%%%%%%% Intance _stages %%%%%%%%
            % PreSleep NREM
            CorrM{j}.PreSWS(idx_toremove_stages,:) = [];
            CorrM{j}.PreSWS(:,idx_toremove_stages) = [];
            % PreSleep REM
            CorrM{j}.PreREM(idx_toremove_stages,:) = [];
            CorrM{j}.PreREM(:,idx_toremove_stages) = [];
            % Exploration + PreTests (Locomotion only)
            CorrM{j}.Task_stages(:,idx_toremove_stages) = [];
            CorrM{j}.Task_stages(idx_toremove_stages,:) = [];
            % PostSleep NREM
            CorrM{j}.PostSWS(idx_toremove_stages,:) = [];
            CorrM{j}.PostSWS(:,idx_toremove_stages) = [];
            % PostSleep REM
            CorrM{j}.PostREM(idx_toremove_stages,:) = [];
            CorrM{j}.PostREM(:,idx_toremove_stages) = [];
            % PostTests (Locomotion only)
            CorrM{j}.PostTest_stages(:,idx_toremove_stages) = [];
            CorrM{j}.PostTest_stages(idx_toremove_stages,:) = [];
            
            %%%%%%%% Intance _splitStages %%%%%%%%
            if ~isempty(splitSleep)
                if splitSleep <= 20
                    % PreSleep NREM split
                    for i=1:numint
                        CorrM{j}.PreSWS_Split{i}(:,idx_toremove_splitStages) = [];
                        CorrM{j}.PreSWS_Split{i}(idx_toremove_splitStages,:) = [];
                    end
                    % Exploration + PreTests (Locomotion only)
                    CorrM{j}.Task_splitStages(:,idx_toremove_splitStages) = [];
                    CorrM{j}.Task_splitStages(idx_toremove_splitStages,:) = [];
                    % PostSleep NREM split
                    for i=1:numint
                        CorrM{j}.PostSWS_Split{i}(idx_toremove_splitStages,:) = [];
                        CorrM{j}.PostSWS_Split{i}(:,idx_toremove_splitStages) = [];
                    end
                    % PostTests (Locomotion only)
                    CorrM{j}.PostTest_stagesSplit(:,idx_toremove_splitStages) = [];
                    CorrM{j}.PostTest_stagesSplit(idx_toremove_splitStages,:) = [];
                end
            end
            
            
            %% Calculate correlation coefficients of correlation matrices
            %%%%%%%% Intance _full %%%%%%%%
            temp = corrcoef(CorrM{j}.Post,CorrM{j}.Task_full); CorrCo_PostTask(j) = temp(1,2); clear temp
            temp = corrcoef(CorrM{j}.Task_full,CorrM{j}.Pre); CorrCo_TaskPre(j) = temp(1,2); clear temp
            temp = corrcoef(CorrM{j}.Post,CorrM{j}.Pre); CorrCo_PostPre(j) = temp(1,2); clear temp
            
            %%%%%%%% Intance _stages %%%%%%%%
            temp = corrcoef(CorrM{j}.PostSWS,CorrM{j}.Task_stages); CorrCo_PostSWSTask(j) = temp(1,2); clear temp
            temp = corrcoef(CorrM{j}.PostREM,CorrM{j}.Task_stages); CorrCo_PostREMTask(j) = temp(1,2); clear temp
            
            temp = corrcoef(CorrM{j}.Task_stages,CorrM{j}.PreSWS); CorrCo_TaskPreSWS(j) = temp(1,2); clear temp
            temp = corrcoef(CorrM{j}.Task_stages,CorrM{j}.PreREM); CorrCo_TaskPreREM(j) = temp(1,2); clear temp
            
            temp = corrcoef(CorrM{j}.PostSWS,CorrM{j}.PreSWS); CorrCo_PostSWSPreSWS(j) = temp(1,2); clear temp
            temp = corrcoef(CorrM{j}.PostREM,CorrM{j}.PreREM); CorrCo_PostREMPreREM(j) = temp(1,2); clear temp
            
            %%%%%%%% Intance _splitStages %%%%%%%%
            if ~isempty(splitSleep)
                if splitSleep <= 20
                    for i=1:numint
                        temp = corrcoef(CorrM{j}.PostSWS_Split{i},CorrM{j}.Task_splitStages);
                        CorrCo_PostSWS_Split_Task{j}(i) = temp(1,2);
                        clear temp
                    end
                    for i=1:numint
                        temp = corrcoef(CorrM{j}.Task_splitStages,CorrM{j}.PreSWS_Split{i});
                        CorrCo_TaskPreSWS_Split{j}(i) = temp(1,2);
                        clear temp
                    end
                    for i=1:numint
                        temp = corrcoef(CorrM{j}.PostSWS_Split{i},CorrM{j}.PreSWS_Split{i}); 
                        CorrCo_PostSWSPreSWS_Split{j}(i) = temp(1,2);
                        clear temp
                    end
                end
            end
            
            %% Calculate EV and REV
            %%%%%%%% Intance _full %%%%%%%%
            EV(j) = ((CorrCo_PostTask(j)-CorrCo_TaskPre(j)*CorrCo_PostPre(j))/...
                sqrt((1-(CorrCo_TaskPre(j))^2)*(1-(CorrCo_PostPre(j))^2)))^2;
            REV(j) = ((CorrCo_TaskPre(j)-CorrCo_PostTask(j)*CorrCo_PostPre(j))/...
                sqrt((1-(CorrCo_PostTask(j))^2)*(1-(CorrCo_PostPre(j))^2)))^2;
            
            %%%%%%%% Intance _stages %%%%%%%%
            EVSWS(j) = ((CorrCo_PostSWSTask(j)-CorrCo_TaskPreSWS(j)*CorrCo_PostSWSPreSWS(j))/...
                sqrt((1-(CorrCo_TaskPreSWS(j))^2)*(1-(CorrCo_PostSWSPreSWS(j))^2)))^2;
            REVSWS(j) = ((CorrCo_TaskPreSWS(j)-CorrCo_PostSWSTask(j)*CorrCo_PostSWSPreSWS(j))/...
                sqrt((1-(CorrCo_PostSWSTask(j))^2)*(1-(CorrCo_PostSWSPreSWS(j))^2)))^2;
            
            EVREM(j) = ((CorrCo_PostREMTask(j)-CorrCo_TaskPreREM(j)*CorrCo_PostREMPreREM(j))/...
                sqrt((1-(CorrCo_TaskPreREM(j))^2)*(1-(CorrCo_PostREMPreREM(j))^2)))^2;
            REVREM(j) = ((CorrCo_TaskPreREM(j)-CorrCo_PostREMTask(j)*CorrCo_PostREMPreREM(j))/...
                sqrt((1-(CorrCo_PostREMTask(j))^2)*(1-(CorrCo_PostREMPreREM(j))^2)))^2;
            
            %%%%%%%% Intance _splitStages %%%%%%%%
            if ~isempty(splitSleep)
                if splitSleep <= 20
                    for i=1:numint
                        EVSWS_Split{i}(j) = ((CorrCo_PostSWS_Split_Task{j}(i)-CorrCo_TaskPreSWS_Split{j}(i)*CorrCo_PostSWSPreSWS_Split{j}(i))/...
                            sqrt((1-(CorrCo_TaskPreSWS_Split{j}(i))^2)*(1-(CorrCo_PostSWSPreSWS_Split{j}(i))^2)))^2;
                        REVSWS_Split{i}(j) = ((CorrCo_TaskPreSWS_Split{j}(i)-CorrCo_PostSWS_Split_Task{j}(i)*CorrCo_PostSWSPreSWS_Split{j}(i))/...
                            sqrt((1-(CorrCo_PostSWS_Split_Task{j}(i))^2)*(1-(CorrCo_PostSWSPreSWS_Split{j}(i))^2)))^2;
                    end
                end
            end
            clear idx_nonexist_Pre idx_nonexist_Task idx_nonexist_Post idx_nonexist_PostTest SWSEpoch REMEpoch
            
        end
    end
end

%% Remove mice without data

QPRE = QPRE(~cellfun('isempty',QPRE));
QPRESWS = QPRESWS(~cellfun('isempty',QPRESWS));
QPREREM = QPREREM(~cellfun('isempty',QPREREM));
QTASK = QTASK(~cellfun('isempty',QTASK));
QPOST = QPOST(~cellfun('isempty',QPOST));
QPOSTSWS = QPOSTSWS(~cellfun('isempty',QPOSTSWS));
QPOSTREM = QPOSTREM(~cellfun('isempty',QPOSTREM));
QPOSTTEST = QPOSTTEST(~cellfun('isempty',QPOSTTEST));

if ~isempty(splitSleep)
    if splitSleep <= 20
        QPRESWS_Split = QPRESWS_Split(~cellfun('isempty',QPRESWS_Split));
        QPOSTSWS_Split = QPOSTSWS_Split(~cellfun('isempty',QPOSTSWS_Split));
    end
end

CorrM = CorrM(~cellfun('isempty',CorrM));

CorrCo_PostTask = nonzeros(CorrCo_PostTask);
CorrCo_TaskPre = nonzeros(CorrCo_PostTask);
CorrCo_PostPre = nonzeros(CorrCo_PostTask);

CorrCo_PostSWSTask = nonzeros(CorrCo_PostSWSTask);
CorrCo_PostREMTask = nonzeros(CorrCo_PostREMTask);
CorrCo_TaskPreSWS = nonzeros(CorrCo_TaskPreSWS);
CorrCo_TaskPreREM = nonzeros(CorrCo_TaskPreREM);
CorrCo_PostSWSPreSWS = nonzeros(CorrCo_PostSWSPreSWS);
CorrCo_PostREMPreREM = nonzeros(CorrCo_PostREMPreREM);

if ~isempty(splitSleep)
    if splitSleep <= 20
        CorrCo_PostSWS_Split_Task = CorrCo_PostSWS_Split_Task(~cellfun('isempty',CorrCo_PostSWS_Split_Task));
        CorrCo_TaskPreSWS_Split = CorrCo_TaskPreSWS_Split(~cellfun('isempty',CorrCo_TaskPreSWS_Split));
        CorrCo_PostSWSPreSWS_Split = CorrCo_PostSWSPreSWS_Split(~cellfun('isempty',CorrCo_PostSWSPreSWS_Split));
    end
end

EV = nonzeros(EV);
REV = nonzeros(REV);
EVSWS = nonzeros(EVSWS);
REVSWS = nonzeros(REVSWS);
EVREM = nonzeros(EVREM);
REVREM = nonzeros(REVREM);

if ~isempty(splitSleep)
    if splitSleep <= 20
        for i=1:numint
            EVSWS_Split{i} = nonzeros(EVSWS_Split{i});
            REVSWS_Split{i} = nonzeros(REVSWS_Split{i});
        end
    end
end

%% Some message for the public
fprintf(['~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n',...
    '         ' num2str(length(EV)) ' mice are in the analysis           \n' ,...
    '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n']);

%% Plot

% Figure
fh = figure('units', 'normalized', 'outerposition', [0 0 0.8 0.5]);

% general explained variance, regardless of sleep state
subplot(1,3,1)
[p_occ,h_occ, her_occ] = PlotErrorBarN_DB([EV*100 REV*100],...
    'barcolors', [1 1 1], 'barwidth', 0.6, 'newfig', 0, 'showpoints',1);
h_occ.FaceColor = 'flat';
h_occ.CData(2,:) = [0 0 0];
set(gca,'Xtick',[1:2],'XtickLabel',{'EV', 'REV'});
set(gca, 'FontSize', 14, 'FontWeight',  'bold');
set(gca, 'LineWidth', 3);
set(h_occ, 'LineWidth', 3);
set(her_occ, 'LineWidth', 3);
ylabel('% explained');
title('Pre-Post sleep EV', 'FontSize', 14);

% general explained variance in nonREM sleep
subplot(1,3,2)
[p_occ,h_occ, her_occ] = PlotErrorBarN_DB([EVSWS*100 REVSWS*100],...
    'barcolors', [1 1 1], 'barwidth', 0.6, 'newfig', 0, 'showpoints',1);
h_occ.FaceColor = 'flat';
h_occ.CData(2,:) = [0 0 0];
set(gca,'Xtick',[1:2],'XtickLabel',{'EV', 'REV'});
set(gca, 'FontSize', 14, 'FontWeight',  'bold');
set(gca, 'LineWidth', 3);
set(h_occ, 'LineWidth', 3);
set(her_occ, 'LineWidth', 3);
ylabel('% explained');
title('Pre-Post NREM sleep EV', 'FontSize', 14);

% general explained variance in REM sleep
subplot(1,3,3)
[p_occ,h_occ, her_occ] = PlotErrorBarN_DB([EVREM*100 REVREM*100],...
    'barcolors', [1 1 1], 'barwidth', 0.6, 'newfig', 0, 'showpoints',1);
h_occ.FaceColor = 'flat';
h_occ.CData(2,:) = [0 0 0];
set(gca,'Xtick',[1:2],'XtickLabel',{'EV', 'REV'});
set(gca, 'FontSize', 14, 'FontWeight',  'bold');
set(gca, 'LineWidth', 3);
set(h_occ, 'LineWidth', 3);
set(her_occ, 'LineWidth', 3);
ylabel('% explained');
title('Pre-Post REM sleep EV', 'FontSize', 14);


%% Plot the figure with dynamics
if ~isempty(splitSleep)
    if splitSleep <= 20
        
        fi = figure('units', 'normalized', 'outerposition', [0 0 0.8 0.5]);
        
        % Prepare the data
        % NREM Sleep
        dat_SWS =  EVSWS_Split{1};
        for i=1:numint
            if i == 1
                dat_SWS = [dat_SWS REVSWS_Split{i}];
                dat_SWS = [dat_SWS zeros(size(dat_SWS,1),1)];
            elseif i == numint
                dat_SWS = [dat_SWS EVSWS_Split{i}];
                dat_SWS = [dat_SWS REVSWS_Split{i}];
            else
                dat_SWS = [dat_SWS EVSWS_Split{i}];
                dat_SWS = [dat_SWS REVSWS_Split{i}];
                dat_SWS = [dat_SWS zeros(size(dat_SWS,1),1)];
            end
        end

        % Plot
        
        [p_occ,h_occ, her_occ] = PlotErrorBarN_DB(dat_SWS*100,...
            'barcolors', [0 0 0], 'barwidth', 0.6, 'newfig', 0, 'showpoints',1);
        h_occ.FaceColor = 'flat';
        h_occ.CData(2,:) = [1 1 1];
        for i=2:3:size(dat_SWS,2)
            h_occ.CData(i,:) = [1 1 1];
        end
        labels = ['0-' num2str(splitSleep) ' min'];
        for i=2:numint
            labels = {labels,[num2str((i-1)*splitSleep) '-' num2str(i*splitSleep) ' min']};
        end
        
        set(gca,'Xtick',[1:3:size(dat_SWS,2)],'XtickLabel',labels);
        set(gca, 'FontSize', 14, 'FontWeight',  'bold');
        set(gca, 'LineWidth', 3);
        set(h_occ, 'LineWidth', 3);
        set(her_occ, 'LineWidth', 3);
        ylabel('% explained');
        title('EV (black) and REV (white) split in intervals', 'FontSize', 14);
    end
end

%% How many mice have EV>REV
idx = EVSWS>REVSWS;
disp(['Number of mice with EV>REV is ' num2str(length(find(idx))) '/' num2str(length(idx))])
%
% % Plot only "working" mice (horrible, I know)
% EV_new = EV(idx);
% REV_new = REV(idx);
%
% Pl = {EV_new*100; REV_new*100};
%
% Cols = {[0.7 0.7 0.7], [0.2 0.2 0.2]};
%
% addpath(genpath('/home/mobsrick/Dima/MatlabToolbox-master/'));
%
% fh = figure('units', 'normalized', 'outerposition', [0 0 0.7 0.9]);
% MakeSpreadAndBoxPlot_SB(Pl,Cols,[1,2]);
% [p,h5,stats] = signrank(Pl{1},Pl{2});
% if p < 0.05
%     sigstar_DB({{1,2}},p,0, 'StarSize',14);
% end
% set(gca,'LineWidth',3,'FontWeight','bold','FontSize',16,'XTick',1:2,'XTickLabel',...
%     {'Explained variance','Reverse explained variance'})
% % ylim([0.15 0.9])
% ylabel('EV or REV in %')
% % title('Place cell stability after conditioning')
%
% rmpath(genpath('/home/mobsrick/Dima/MatlabToolbox-master/'));
%
%
% if sav
%     saveas(gcf,[dropbox '/MOBS_workingON/Dima/Ongoing results/PlaceField_Final/ExplainedVariance.fig']);
%     saveFigure(gcf,'ExplainedVariance',...
%         [dropbox '/MOBS_workingON/Dima/Ongoing results/PlaceField_Final/']);
% end