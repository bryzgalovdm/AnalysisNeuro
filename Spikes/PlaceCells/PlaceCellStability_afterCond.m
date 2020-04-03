%%%%% PlaceCellStability_afterCond
% 
% This script calculates correlation between the rate maps
% of pre-exploration (15 min long) and 4 post-tests concatenated together
% (4 * 4 min = 16 min overall). Control is correlation between the first half
% and the second half of pre-exploration (7,5 min each)
% 
% Last touched on 03/04/2020
% 
% By Dima Bryzgalov, MOBS team, Paris
%
% TODO: you need to track each place cell
% get the cells with the same cell stability <0.3

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
sav = 0;

%% How many cells do you have
numPCs = CountPlaceCells(Dir);

%% Main loop

% Loop counter (who knows better practice?)
cnt=0;

% Allocate memory
CorrStabSameCell = nan(1,numPCs);
CorrStabAfterCondCell = nan(1,numPCs);

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
            try % this try-catch syntax avoids spitting out an error in a situation when there are too less spikes for analysis on one half of the recording
                map1 = PlaceField_DB(Restrict(spikes.S{spikes.PlaceCells.idx(i)},MovingHab1Half),... % 1st half
                    Restrict(beh.CleanAlignedXtsd, MovingHab1Half),...
                    Restrict(beh.CleanAlignedYtsd, MovingHab1Half), 'smoothing', 2.5, 'size', 50, 'plotresults',0, 'plotpoisson', 0);
                map2 = PlaceField_DB(Restrict(spikes.S{spikes.PlaceCells.idx(i)},MovingHab2Half),... % 2nd half
                    Restrict(beh.CleanAlignedXtsd, MovingHab2Half),...
                    Restrict(beh.CleanAlignedYtsd, MovingHab2Half), 'smoothing', 2.5, 'size', 50, 'plotresults',0, 'plotpoisson', 0);
            catch
                map1 = [];
                map2 = [];
            end
            % Correlation coefficient
            if ~isempty(map1) && ~isempty(map2)
                cnt=cnt+1;
                de = corrcoef(map1.rate,map2.rate);
                CorrStabSameCell(cnt) = de(2,1);
            else
                cnt=cnt+1;
            end
        
        clear map1 map2
        
        % Inter-cell correlation
            try % this try-catch syntax avoids spitting out an error in a situation when there are too less spikes for analysis on one half of the recording
                map1 = PlaceField_DB(Restrict(spikes.S{spikes.PlaceCells.idx(i)},UMazeMovingEpoch),... % Exploration before aversive experience
                    Restrict(beh.CleanAlignedXtsd, UMazeMovingEpoch),...
                    Restrict(beh.CleanAlignedYtsd, UMazeMovingEpoch), 'smoothing', 2.5, 'size', 50, 'plotresults',0, 'plotpoisson', 0);
                map2 = PlaceField_DB(Restrict(spikes.S{spikes.PlaceCells.idx(i)},AfterConditioningMovingEpoch),... % Exploration after aversive experience
                    Restrict(beh.CleanAlignedXtsd, AfterConditioningMovingEpoch),...
                    Restrict(beh.CleanAlignedYtsd, AfterConditioningMovingEpoch), 'smoothing', 2.5, 'size', 50, 'plotresults', 0, 'plotpoisson', 0);
            catch
                map1 = [];
                map2 = [];
            end
            % Correlation coefficent
            if ~isempty(map1) && ~isempty(map2)
                de = corrcoef(map1.rate,map2.rate);
                CorrStabAfterCondCell(cnt) = de(2,1);
            end
        end
    end
    clear spikes beh st en LocomotionEpoch MovingHab1Half MovingHab2Half AfterConditioningMovingEpoch BigMazeMovingEpoch UMazeMovingEpoch
end

%% Create two types of data to plot

% First, pairwise matching of within-session and after-conditioning correlation coefficients
StabilityAcrossExperiment = [CorrStabSameCell' CorrStabAfterCondCell']; % first column - within session, second - across aversive learning

idxnan = isnan(StabilityAcrossExperiment);
idxdel = false(1,numPCs);
for i = 1:length(col1nan)
    if idxnan(i,1) || idxnan(i,2)
        idxdel(i) = true;
    end
end
StabilityAcrossExperiment(idxdel,:) = [];

% Second, all correlations calculated (arrays are not equal in length)
CorrStabSameCell_all = CorrStabSameCell(~isnan(CorrStabSameCell));
CorrStabAfterCondCell_all = CorrStabAfterCondCell(~isnan(CorrStabAfterCondCell));

%% Plot boxplot
Pl = {CorrStabAfterCondCell_all; CorrStabSameCell_all};

Cols = {[0.7 0.7 0.7], [0.2 0.2 0.2]};

fh = figure('units', 'normalized', 'outerposition', [0 0 0.5 0.7]);
MakeSpreadAndBoxPlot_SB(Pl,Cols,[1,2],{},1);
set(gca,'LineWidth',3,'FontWeight','bold','FontSize',16,'XTick',1:2,'XTickLabel',{'PreVsPost','WithinBaselineSession'})
ylim([0 0.9])
ylabel('Correlation Coef.')
title('Place cell stability after conditioning')

% saveas(fh,'/home/mobsrick/Dropbox/MOBS_workingON/Dima/Ongoing results/PlaceField_Final/Stability_PrePostAllCells.fig');
% saveFigure(fh,'Stability_PrePostAllCells','/home/mobsrick/Dropbox/MOBS_workingON/Dima/Ongoing results/PlaceField_Final/');

%% Plot barplot

% PlotErrorBarN_DB(StabilityAcrossExperiment)