%%%%% PlaceCellStability_afterCond
% 
% This script calculates correlation between the rate maps
% of pre-exploration (15 min long) and 4 post-tests concatenated together
% (4 * 4 min = 16 min overall). Control is correlation between the first half
% and the second half of pre-exploration (7,5 min each)
% 
% It plots out only cells with within-session stability <= 0.5
% 
% Last touched on 23/04/2020
% 
% By Dima Bryzgalov, MOBS team, Paris
% 

%% Parameters

% Mice in the analysis
nmouse = [797 798 828 861 882 905 906 911 912 977 994];
% nmouse = [906 912]; % Had PreMazes
% nmouse = [905 911]; % Did not have PreMazes

% Paths retrieved
% Dir = PathForExperimentsERC_Dima('UMazePAG'); 
Dir = PathForExperimentsERC_DimaMAC('UMazePAG');
Dir = RestrictPathForExperiment(Dir,'nMice',nmouse);

% Threshold on speed (im cm/s, epoch with lower are not considered)
speed_thresh = 3;

% Do you want to save a figure?
sav = false;

% Paths and names to save
pathfig = '/MOBS_workingON/Dima/Ongoing_results/PlaceField_Final/Stability/'; % without dropbox path
figbox = 'Stability_PrePost_box';
figbar = 'Stability_PrePost_bar';

%% How many cells do you have
numPCs = CountPlaceCells(Dir);

%% Main loop

% Loop counter (who knows better practice?)
cnt=0;

% Allocate memory
CorrStabSameCell = nan(1,numPCs);
CorrStabAfterCondCell = nan(1,numPCs);
idPC = cell(1, numPCs);
idxSZ = false(1, numPCs);

for j=1:length(Dir.path)
    
    % Get the data
    cd(Dir.path{j}{1});
    spikes = load('SpikeData.mat','S','PlaceCells');
    beh = load('behavResources.mat','SessionEpoch','CleanVtsd','CleanAlignedXtsd','CleanAlignedYtsd');
    
    % Split the exploration in two
    st = Start(beh.SessionEpoch.Hab);
    en = End(beh.SessionEpoch.Hab);
    Hab1Half = intervalSet(st,((en-st)/2)+st);
    Hab2Half = intervalSet(((en-st)/2)+st,en);
    
    % Create BaselineExplo Epoch (NOT TO BE USED! SO FAR)
    BigMazeEpoch = or(beh.SessionEpoch.Hab,beh.SessionEpoch.TestPre1);
    BigMazeEpoch = or(BigMazeEpoch,beh.SessionEpoch.TestPre2);
    BigMazeEpoch = or(BigMazeEpoch,beh.SessionEpoch.TestPre3);
    BigMazeEpoch = or(BigMazeEpoch,beh.SessionEpoch.TestPre4);
    
    % After Conditioning
    AfterConditioningEpoch = or(beh.SessionEpoch.TestPost1,beh.SessionEpoch.TestPost2);
    AfterConditioningEpoch = or(AfterConditioningEpoch,beh.SessionEpoch.TestPost3);
    AfterConditioningEpoch = or(AfterConditioningEpoch,beh.SessionEpoch.TestPost4);
    
    % Locomotion threshold
    LocomotionEpoch = thresholdIntervals(tsd(Range(beh.CleanVtsd),movmedian(Data(beh.CleanVtsd),5)),speed_thresh,'Direction','Above'); % smoothing = 5
    
    % Get resulting epochs
    BigMazeMovingEpoch = and(LocomotionEpoch, BigMazeEpoch); % NOT TO BE USED! SO FAR
    UMazeMovingEpoch = and(beh.SessionEpoch.Hab, LocomotionEpoch);
    AfterConditioningMovingEpoch = and(LocomotionEpoch, AfterConditioningEpoch);
    MovingHab1Half = and(LocomotionEpoch,Hab1Half);
    MovingHab2Half = and(LocomotionEpoch,Hab2Half);
    
    % Calculate same cell stabity
    if isfield(spikes.PlaceCells,'idx')
        for i=1:length(spikes.PlaceCells.idx)
            cnt=cnt+1;
            idPC{cnt} = [j spikes.PlaceCells.idx(i)];
            % Grab shock zone place cells
            if sum(idPC{cnt}(2) == spikes.PlaceCells.SZ)
                idxSZ(cnt) = true;
            end
            try % this try-catch syntax avoids spitting out an error in a situation when there are too less spikes for analysis on one half of the recording
                map1=PlaceField_DB(spikes.S{spikes.PlaceCells.idx(i)},beh.CleanAlignedXtsd, beh.CleanAlignedYtsd,... % 1st half
                    'Epoch', MovingHab1Half, 'smoothing', 2.5, 'size', 50, 'plotresults',0, 'plotpoisson', 0);
                map2=PlaceField_DB(spikes.S{spikes.PlaceCells.idx(i)},beh.CleanAlignedXtsd, beh.CleanAlignedYtsd,... % 2nd half
                    'Epoch', MovingHab2Half, 'smoothing', 2.5, 'size', 50, 'plotresults',0, 'plotpoisson', 0);
            catch
                map1 = [];
                map2 = [];
            end
            % Correlation coefficient
            if ~isempty(map1) && ~isempty(map2)
                de = corrcoef(map1.rate,map2.rate);
                CorrStabSameCell(cnt) = de(2,1);
            end
        
        clear map1 map2
        
        % Inter-cell correlation
            try % this try-catch syntax avoids spitting out an error in a situation when there are too less spikes for analysis on one half of the recording
                map1=PlaceField_DB(spikes.S{spikes.PlaceCells.idx(i)},beh.CleanAlignedXtsd, beh.CleanAlignedYtsd,... % Exploration before aversive experience
                    'Epoch', UMazeMovingEpoch, 'smoothing', 2.5, 'size', 50, 'plotresults',0, 'plotpoisson', 0);
                map2=PlaceField_DB(spikes.S{spikes.PlaceCells.idx(i)},beh.CleanAlignedXtsd, beh.CleanAlignedYtsd,... % Exploration after aversive experience
                    'Epoch', AfterConditioningMovingEpoch, 'smoothing', 2.5, 'size', 50, 'plotresults',0, 'plotpoisson', 0);
            catch
                map1 = [];
                map2 = [];
            end
            % Correlation coefficent
            if ~isempty(map1) && ~isempty(map2)
                de = corrcoef(map1.rate,map2.rate);
                CorrStabAfterCondCell(cnt) = de(2,1);
            end
            
            clear map1 map2
        end
    end
    clear spikes beh st en LocomotionEpoch MovingHab1Half MovingHab2Half AfterConditioningMovingEpoch BigMazeMovingEpoch UMazeMovingEpoch
end

%% Create two types of data to plot

% First, pairwise matching of within-session and after-conditioning correlation coefficients
StabilityAcrossExperiment = [CorrStabSameCell' CorrStabAfterCondCell']; % first column - within session, second - across aversive learning

idxnan = isnan(StabilityAcrossExperiment);
idxdel_even = false(1,numPCs);
for i = 1:length(idxnan)
    if idxnan(i,1) || idxnan(i,2)
        idxdel_even(i) = true;
    end
end
StabilityAcrossExperiment_even = StabilityAcrossExperiment;
StabilityAcrossExperiment_even(idxdel_even,:) = [];

% Second, all correlations calculated (arrays are not equal in length)
CorrStabSameCell_all = CorrStabSameCell(~isnan(CorrStabSameCell));
CorrStabAfterCondCell_all = CorrStabAfterCondCell(~isnan(CorrStabAfterCondCell));

%% Grab instances in which place fields are not found or within session stability < 0.5
badcases = cell(1,numPCs);
idxdel_bad = false(1,numPCs);
for i = 1:numPCs
    if isnan(CorrStabSameCell(i)) || CorrStabSameCell(i) < 0.5
        badcases{i} = idPC{i};
        idxdel_bad(i) = true;
    end
end
badcases = badcases(~cellfun('isempty',badcases));

% Remove all bad cases
idxdel = or(idxdel_even, idxdel_bad);
idxdelSZ = or(or(idxdel_even,idxdel_bad),idxSZ);
idxbigSZ = and(not(or(idxdel_even,idxdel_bad)),idxSZ);

StabilityAcrossExperiment_clean = StabilityAcrossExperiment;
StabilityAcrossExperiment_clean(idxdel,:) = [];

% NoShockZones
StabilityAcrossExperiment_clean_NSZ = StabilityAcrossExperiment;
StabilityAcrossExperiment_clean_NSZ(idxdelSZ,:) = [];

% ShockZones
StabilityAcrossExperiment_clean_SZ = StabilityAcrossExperiment;
StabilityAcrossExperiment_clean_SZ = StabilityAcrossExperiment_clean_SZ(idxbigSZ,:);


%% Grab shock zone cells
% SZcells = cell(1,numPCs);
% idxdel_SZ = false(1,numPCs);
% for i = 1:numPCs
%     if isnan(CorrStabSameCell(i)) || CorrStabSameCell(i) < 0.5
%         badcases{i} = idPC{i};
%         idxdel_bad(i) = true;
%     end
% end
% badcases = badcases(~cellfun('isempty',badcases));

%% Plot bad cases
% for i=1:length(badcases)
%     cd(Dir.path{badcases{i}(1)}{1});
%     load ('SpikeData.mat', 'S', 'PlaceCells');
%     load('behavResources.mat','SessionEpoch','CleanVtsd','CleanAlignedXtsd','CleanAlignedYtsd');
%     % Locomotion threshold
%     LocomotionEpoch = thresholdIntervals(tsd(Range(CleanVtsd),movmedian(Data(CleanVtsd),5)),2.5,'Direction','Above'); % smoothing = 5   
%     % Get resulting epochs
%     UMazeMovingEpoch = and(SessionEpoch.Hab, LocomotionEpoch);
%     map = PlaceField_DB(Restrict(S{badcases{i}(2)},UMazeMovingEpoch),... % 1st half
%                     Restrict(CleanAlignedXtsd, UMazeMovingEpoch),...
%                     Restrict(CleanAlignedYtsd, UMazeMovingEpoch), 'smoothing', 2.5, 'size', 50, 'plotresults',0, 'plotpoisson', 0);
% %     map = PlaceField_DB(Restrict(S{badcases{i}(2)},SessionEpoch.Hab),... % 1st half
% %                     Restrict(CleanAlignedXtsd, SessionEpoch.Hab),...
% %                     Restrict(CleanAlignedYtsd, SessionEpoch.Hab), 'smoothing', 2.5, 'size', 50, 'plotresults',1, 'plotpoisson', 1);
%     figure
%     imagesc(map.rate);
%     colormap jet
%     title([Dir.name{badcases{i}(1)} '  PC#' num2str(badcases{i}(2))]);
%     
%     clear S PlaceCells SessionEpoch CleanVtsd CleanAlignedXtsd CleanAlignedYtsd LocomotionEpoch UMazeMovingEpoch map
%     
% end
    

%% Plot boxplot
% Pl = {CorrStabAfterCondCell_all; CorrStabSameCell_all};
% Pl = {StabilityAcrossExperiment_clean(:,2); StabilityAcrossExperiment_clean(:,1)};
% 
% Cols = {[0.2 0.2 0.2], [0.7 0.7 0.7]};
% 
% fh = figure('units', 'normalized', 'outerposition', [0 0 0.5 0.7]);
% MakeSpreadAndBoxPlot_SB(Pl,Cols,[1,2],{},1);
% set(gca,'LineWidth',3,'FontWeight','bold','FontSize',16,'XTick',1:2,'XTickLabel',{'PreVsPost','WithinBaselineSession'})
% ylim([0 0.9])
% ylabel('Correlation Coef.')
% title('Place cell stability after conditioning')

% if sav
%     saveas(fh,[pathfig figbox '.fig']);
%     saveFigure(fh,figbox,pathfig);
% end

%% Plot barplot

fb = figure('units', 'normalized', 'outerposition', [0 0 0.5 0.7]);
[p_s,h_s, her_s] = PlotErrorBarN_DB(StabilityAcrossExperiment_clean, 'barcolors', [0 0 0], 'barwidth', 0.6, 'newfig', 0, 'showpoints',1);
h_s.FaceColor = 'flat'; 
h_s.CData(2,:) = [1 1 1]; % It does not work in 2016b
set(gca,'Xtick',[1:2],'XtickLabel',{'WithinSession', 'PreVsPost'});
set(gca, 'FontSize', 14, 'FontWeight',  'bold');
set(gca, 'LineWidth', 3);
set(h_s, 'LineWidth', 3);
set(her_s, 'LineWidth', 3);
ylabel('Correlation coeff.');
title('Place cells stability', 'FontSize', 14);
ylim([-0.15 1.1])

if sav
    saveas(fb,[dropbox pathfig figbar '.fig']);
    saveFigure(fb,figbar,[dropbox pathfig]);
end

%% Plot distribution of within-session stability

% figure
% hist(CorrStabSameCell_all,10)