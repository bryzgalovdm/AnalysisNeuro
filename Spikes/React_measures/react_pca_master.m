function react_pca_master(nmice, experiment, TemplatePeriod, numtemplates, varargin)
%
% The function calculates reactivation strength for the ERC experiment
% in the same manner it was calculated in Peyrache et al., 2010
%
%
% INPUT
%
%   nmice               array with mice numbers that will harvested from ERC
%                       PathForExperiments. Each mouse should contain
%                       PlaceCells structure in its SpikeData.mat
%   experiment          type of experiment: 'PAG', 'MFB', 'Novel'
%   TemplatePeriod      type of period to be taken as a template. The
%                       following periods supported: 'wake', 'cond', 'condFree', 'postRip'
%
%
%
%   SpeedThresh         threshold (in cm/s) to calculate place fields
%                       (default = 4) (optional)
%   PlotResults         whether to plot results or not (optional - default=true)
%   BinSize             binsize (in tsd units) used to creat spike-time
%                       histograms (default = 0.1*1e4). In RepeatOPostRipplesEpochriginalPaper
%                       mode equals 10 (1 ms) (optional)
%   SplitSleep          Splits 40 min of NREM sleep into n intervals, each
%                       SplitSleep min long. If empty, no split happens (default=20)
%   IncludeInterneurons if false, removes all interneurons (default = true)
%   IsSaveFig           if true, saves figures in dropbox (default = false)(optional)
%
% OUTPUT
%
%  None
%
% EXAMPLE
%
%
% By Dima Bryzgalov, MOBS team, Paris,
% 04/10/2021 based on the codes of Karim Benchenane
% github.com/bryzgalovdm


%% Default optional arguments
binsize = 0.04*1e4;
IsPlotInd = false;

%% Parse parameters
for i=1:2:length(varargin)
    switch(lower(varargin{i}))
        case 'binsize'
            binsize = varargin{i+1};
            if ~isa(binsize, 'numeric')
                error('Incorrect value for property ''BinSize'' (type ''help react_pca_master'' for details).');
            end
        case 'issavefig'
            IsSaveFig = varargin{i+1};
            if IsSaveFig ~= 1 && IsSaveFig ~= 0
                error('Incorrect value for property ''IsSaveFig'' (type ''help ExplainedVariance_master_DB'' for details).');
            end
        case 'plotresults'
            IsPlot = varargin{i+1};
            if IsPlot ~= 1 && IsPlot ~= 0
                error('Incorrect value for property ''PlotResults'' (type ''help ExplainedVariance_master_DB'' for details).');
            end
        case 'plotindiv'
            IsPlotInd = varargin{i+1};
            if IsPlotInd ~= 1 && IsPlotInd ~= 0
                error('Incorrect value for property ''PlotIndiv'' (type ''help react_pca_master'' for details).');
            end
        case 'speedthresh'
            speed_thresh = varargin{i+1};
            if ~isa(speed_thresh, 'numeric') && speed_thresh <= 0
                error('Incorrect value for property ''SpeedThresh'' (type ''help ExplainedVariance_master_DB'' for details).');
            end
        case 'includeinterneurons'
            IsII = varargin{i+1};
            if ~isa(IsII, 'numeric')
                if IsII ~= 1 && IsII ~= 0
                    error('Incorrect value for property ''IncludeInterneurons'' (type ''help ExplainedVariance_master_DB'' for details).');
                end
            end
    end
end

%% Manage inputs
% Experiment
if strcmp(experiment, 'PAG')
    fetchpaths = 'UMazePAG';
elseif strcmp(experiment, 'MFB')
    fetchpaths = 'StimMFBWake';
elseif strcmp(experiment, 'Novel')
    fetchpaths = 'Novel';
end

% Analysis period
if strcmp(TemplatePeriod, 'wake')
    template = 'FullTaskEpoch';
elseif strcmp(TemplatePeriod, 'cond')
    template = 'CondEpoch';
elseif strcmp(TemplatePeriod, 'condFree')
    template = 'CondFreezeEpoch';
elseif strcmp(TemplatePeriod, 'postRip')
    template = 'PostRipplesEpoch';
elseif strcmp(TemplatePeriod, 'condRip')
    template = 'CondRipplesEpoch';
else
    error('TemplatePeriod has the wrong value. Please see help for more info')
end


%% Get paths of each individual mouse
Dir = PathForExperimentsERC(fetchpaths);
Dir = RestrictPathForExperiment(Dir,'nMice',nmice);
[numsessions, micenames] = CountNumSesionsERC(Dir);

%% Allocate data
s = cell(numsessions, 1);
b = cell(numsessions, 1);
f = cell(numsessions, 1);
r = cell(numsessions, 1);
sleep = cell(numsessions, 1);

RipEpoch = cell(numsessions, 1);

Wake = cell(numsessions, 1);
REMEpoch = cell(numsessions, 1);
SWSEpoch = cell(numsessions, 1);

HabEpoch = cell(numsessions, 1);
TestPreEpoch = cell(numsessions, 1);
TestPostEpoch = cell(numsessions, 1);
CondEpoch = cell(numsessions, 1);
FullTaskEpoch = cell(numsessions, 1);
CondFreezeEpoch = cell(numsessions, 1);
PreRipplesEpoch = cell(numsessions, 1);
PostRipplesEpoch = cell(numsessions, 1);
CondRipplesEpoch = cell(numsessions, 1);

Q = cell(numsessions, 1);
PCtemplates = cell(numsessions, 1);
replayGtsd = cell(numsessions, 1);
EpochXY = cell(numsessions, 1);
LocPCHab = cell(numsessions, 1);
LocPCCond = cell(numsessions, 1);
tRipples = cell(numsessions, 1);

SWRrate = nan(numsessions, 4); % 4 columns: PreSleep, Hab, Cond, PostSleep

%% Load data
cnt=1;
for imouse=1:length(Dir.path)
    for isession = 1:length(Dir.path{imouse})
        s{cnt} = load([Dir.path{imouse}{isession} 'SpikeData.mat'], 'S');
        b{cnt} = load([Dir.path{imouse}{isession} 'behavResources.mat'],'SessionEpoch', 'Vtsd',...
            'AlignedXtsd', 'AlignedYtsd', 'ZoneEpoch', 'TTLInfo');
        if strcmp(experiment, 'PAG')
            f{cnt} = load([Dir.path{imouse}{isession} 'behavResources.mat'], 'FreezeAccEpoch');
        end
        try
            r{cnt} = load([Dir.path{imouse}{isession} 'SWR.mat'],'ripples');
        catch
            r{cnt} = load([Dir.path{imouse}{isession} 'Ripples.mat'],'ripples');
        end
        try
            sleep{cnt} = load([Dir.path{imouse}{isession} 'SleepScoring_OBGamma.mat'],'SWSEpoch','REMEpoch','Wake', 'TotalNoiseEpoch');
        catch
            sleep{cnt} = load([Dir.path{imouse}{isession} 'SleepScoring_Accelero.mat'],'SWSEpoch','REMEpoch','Wake', 'TotalNoiseEpoch');
        end
        cnt=cnt+1;
    end
end

%% Prepare data
for isession = 1:numsessions
    RipEpoch{isession} = intervalSet(r{isession}.ripples(:,2)*1E4-0.2E4, r{isession}.ripples(:,2)*1E4+0.2E4);
    
    Wake{isession} = sleep{isession}.Wake - sleep{isession}.TotalNoiseEpoch;
    REMEpoch{isession} = sleep{isession}.REMEpoch - sleep{isession}.TotalNoiseEpoch;
    SWSEpoch{isession} = sleep{isession}.SWSEpoch - sleep{isession}.TotalNoiseEpoch;
    
    if strcmp(experiment, 'PAG')
        [Hab, TestPre, ~, Cond, FullTask, TestPost] = ReturnMnemozyneEpochs(b{isession}.SessionEpoch);
    elseif strcmp(experiment, 'MFB')
        [Hab, TestPre, ~, Cond, FullTask, TestPost] = ReturnMnemozyneEpochs(b{isession}.SessionEpoch, 'NumberTests', 8);
    elseif strcmp(experiment, 'Novel')
        [Hab, TestPre, ~, Cond, FullTask, TestPost] = ReturnMnemozyneEpochs(b{isession}.SessionEpoch, 'NumberTests', 4);
    end
    
    % Remove noise and create possible template epochs
    HabEpoch{isession} = Hab - sleep{isession}.TotalNoiseEpoch;
    TestPreEpoch{isession} = TestPre - sleep{isession}.TotalNoiseEpoch;
    TestPostEpoch{isession} = TestPost - sleep{isession}.TotalNoiseEpoch;
    CondEpoch{isession} = Cond - sleep{isession}.TotalNoiseEpoch;
    FullTaskEpoch{isession} = FullTask - sleep{isession}.TotalNoiseEpoch;
    if strcmp(experiment, 'PAG')
        CondFreezeEpoch{isession} = and(Cond, f{isession}.FreezeAccEpoch) - sleep{isession}.TotalNoiseEpoch;
    end
    PreRipplesEpoch{isession} = and(and(b{isession}.SessionEpoch.PreSleep, sleep{isession}.SWSEpoch), RipEpoch{isession}) - ...
        sleep{isession}.TotalNoiseEpoch;
    PostRipplesEpoch{isession} = and(and(b{isession}.SessionEpoch.PostSleep, sleep{isession}.SWSEpoch), RipEpoch{isession}) - ...
        sleep{isession}.TotalNoiseEpoch;
    CondRipplesEpoch{isession} = and(Cond, RipEpoch{isession}) - sleep{isession}.TotalNoiseEpoch;
    
end

%% Create templates
for isession = 1:numsessions
    Q{isession} = MakeQfromS(s{isession}.S, binsize);
    % Create tempate epoch and get pcs
    eval(['Qtemplate = full(Data(Restrict(Q{isession}, ' template '{isession})));']);
    C = corrcoef(Qtemplate);
    C(isnan(C))=0;
    C=C-diag(diag(C));
    [PCtemplates{isession},L] = pcacov(C);
    
    % Elbow rule?
    % figure, plot(L,'o-');
    % xlabel('Number of PC');
    % ylabel('Eigenvalue')
    % makepretty
end

%% Match it (to the full experiment)
for isession = 1:numsessions
    Qf=full(Data(Q{isession}));
    for itemplate = 1:numtemplates
        scoreG = zscore(Qf)*PCtemplates{isession}(:,itemplate);
        singleNG = (zscore(Qf).^2)*(PCtemplates{isession}(:,itemplate).^2);
        replayG = scoreG.^2 - singleNG;
        replayGtsd{isession}{itemplate}=tsd(Range(Q{isession}),replayG);
    end
end

%% Location of RS
% Bin the maze
for isession = 1:numsessions
    i=1;
    for x=0:0.1:1
        EpochXd=thresholdIntervals(b{isession}.AlignedXtsd,x,'Direction','Above');
        EpochXu=thresholdIntervals(b{isession}.AlignedXtsd,x+0.05,'Direction','Below');
        EpochX=and(EpochXd,EpochXu);
        j=1;
        for y=0:0.1:1
            EpochYd=thresholdIntervals(b{isession}.AlignedYtsd,y,'Direction','Above');
            EpochYu=thresholdIntervals(b{isession}.AlignedYtsd,y+0.05,'Direction','Below');
            EpochY=and(EpochYd,EpochYu);
            EpochXY{isession}{i,j}=and(EpochX,EpochY);
            j=j+1;
        end
        i=i+1;
    end
end

% Put RS on the UMaze
for isession = 1:numsessions
    for itemplate = 1:numtemplates
        for i=1:size(EpochXY{isession},1)
            for j=1:size(EpochXY{isession},2)
                LocPCHab{isession}{itemplate}(i,j)=nanmean(Data(Restrict(replayGtsd{isession}{itemplate},...
                    and(HabEpoch{isession},EpochXY{isession}{i,j}))));
                LocPCCond{isession}{itemplate}(i,j)=nanmean(Data(Restrict(replayGtsd{isession}{itemplate},...
                    and(CondEpoch{isession},EpochXY{isession}{i,j}))));
            end
        end
    end
end

%% Gather all-mice arrays
for isession = 1:numsessions
    tRipples{isession} = ts(r{isession}.ripples(:,2)*1E4);
    nrip_presleep = Range(Restrict(tRipples{isession}, b{isession}.SessionEpoch.PreSleep));
    nrip_hab = Range(Restrict(tRipples{isession}, HabEpoch{isession}));
    nrip_cond = Range(Restrict(tRipples{isession}, CondEpoch{isession}));
    nrip_postsleep = Range(Restrict(tRipples{isession}, b{isession}.SessionEpoch.PostSleep));
    
    SWRrate(isession, 1)=length(nrip_presleep)/sum(tot_length(and(b{isession}.SessionEpoch.PreSleep, SWSEpoch{isession}),'s'));
    SWRrate(isession, 2)=length(nrip_hab)/sum(tot_length(HabEpoch{isession},'s'));
    SWRrate(isession, 3)=length(nrip_cond)/sum(tot_length(CondEpoch{isession},'s'));
    SWRrate(isession, 4)=length(nrip_postsleep)/sum(tot_length(and(b{isession}.SessionEpoch.PostSleep, SWSEpoch{isession}),'s'));
    
    for itemplate = 1:numtemplates
        % Ripples-triggered RS
        [meanRSripPr(isession, itemplate, :),~,tPr]=ETAverage(nrip_presleep, Range(replayGtsd{isession}{itemplate}),...
            Data(replayGtsd{isession}{itemplate}), 100,50);
        [meanRSripCo(isession, itemplate, :),~, tCo]=ETAverage(nrip_cond, Range(replayGtsd{isession}{itemplate}),...
            Data(replayGtsd{isession}{itemplate}), 100,50);
        [meanRSripPo(isession, itemplate, :),~,tPo]=ETAverage(nrip_postsleep, Range(replayGtsd{isession}{itemplate}),...
            Data(replayGtsd{isession}{itemplate}), 100,50);
        % Mean RS
        meanRsPr(isession, itemplate) = mean(Data(Restrict(replayGtsd{isession}{itemplate},b{isession}.SessionEpoch.PreSleep)));
        meanRsHab(isession, itemplate) = mean(Data(Restrict(replayGtsd{isession}{itemplate},HabEpoch{isession})));
        meanRsCo(isession, itemplate) = mean(Data(Restrict(replayGtsd{isession}{itemplate},CondEpoch{isession})));
        meanRsPo(isession, itemplate) = mean(Data(Restrict(replayGtsd{isession}{itemplate},b{isession}.SessionEpoch.PostSleep)));
        
        meanRSCoFreeze(isession, itemplate) = mean(Data(Restrict(replayGtsd{isession}{itemplate},CondFreezeEpoch{isession})));
        meanRSCoRipples(isession, itemplate) = mean(Data(Restrict(replayGtsd{isession}{itemplate},CondRipplesEpoch{isession})));
        meanRSCoNoRipples(isession, itemplate) = mean(Data(Restrict(replayGtsd{isession}{itemplate},...
            CondEpoch{isession} - CondRipplesEpoch{isession})));
        
    end
end

%% Plot individual figure
if IsPlotInd
    for isession = 1:numsessions
        for itemplate = 1:numtemplates
            finv{isession}{itemplate} = plotIndividFigure_react_pca(tRipples{isession}, b{isession}.SessionEpoch,...
                HabEpoch{isession}, CondEpoch{isession}, sleep{isession}.SWSEpoch, ...
                replayGtsd{isession}{itemplate}, LocPCHab{isession}{itemplate}, LocPCCond{isession}{itemplate},...
                SWRrate(isession, :), squeeze(meanRSripPr(isession, itemplate, :)), squeeze(meanRSripCo(isession, itemplate, :)),...
                squeeze(meanRSripPo(isession, itemplate, :)), tPr, tCo, tPo, micenames{isession}, itemplate);
        end
    end
end


%% General figures

% RS triggered on ripples
fh1 = figure('color',[1 1 1], 'units', 'normalized', 'outerposition', [0 0 0.8 0.65]);
subplot(1,3,1), hold on
for itemplate = 1:size(meanRSripPr,2)
    plot(tPr,squeeze(meanRSripPr(:,itemplate,:)), 'Color', [.8 .8 .8]), ylim([-2 20]), xlim([tPr(1) tPr(end)]);
end
xlabel('Time (ms)')
ylabel('Reactivation strength')
line([0 0],ylim,'color','k','linestyle',':')
makepretty_DB
subplot(1,3,2), hold on
for itemplate = 1:size(meanRSripPo,2)
    plot(tPo,squeeze(meanRSripPo(:,itemplate,:)),'k'), ylim([-2 20]), xlim([tPo(1) tPo(end)]);
end
xlabel('Time (ms)')
ylabel('Reactivation strength')
line([0 0],ylim,'color','k','linestyle',':')
makepretty_DB
subplot(1,3,3), hold on, 
plot(tPr,squeeze(nanmean(nanmean(meanRSripPr, 1), 2)),'Color', [.8 .8 .8]);
plot(tCo,squeeze(nanmean(nanmean(meanRSripCo, 1), 2)),'r');
plot(tPo,squeeze(nanmean(nanmean(meanRSripPo, 1), 2)),'k'); 
xlabel('Time (ms)')
ylabel('Reactivation strength')
line([0 0],ylim,'color','k','linestyle',':'), xlim([tPo(1) tPo(end)])
legend('PreSleep', 'PostSleep', 'Learning')
makepretty_DB

% Scatters RS in different periods
fh2 = figure('color',[1 1 1], 'units', 'normalized', 'outerposition', [0 0 0.8 0.8]);
% Prepare data
meanPr_plot = reshape(meanRsPr, [size(meanRsPr,1)*size(meanRsPr,2), 1]);
meanPo_plot = reshape(meanRsPo, [size(meanRsPo,1)*size(meanRsPo,2), 1]);
meanCo_plot = reshape(meanRsCo, [size(meanRsCo,1)*size(meanRsCo,2), 1]);
meanHab_plot = reshape(meanRsHab, [size(meanRsHab,1)*size(meanRsHab,2), 1]);
meanCoFreeze_plot = reshape(meanRSCoFreeze, [size(meanRSCoFreeze,1)*size(meanRSCoFreeze,2), 1]);
meanCoRipples_plot = reshape(meanRSCoRipples, [size(meanRSCoRipples,1)*size(meanRSCoRipples,2), 1]);
meanCoNoRipples_plot = reshape(meanRSCoNoRipples, [size(meanRSCoNoRipples,1)*size(meanRSCoNoRipples,2), 1]);
% plot scatter
subplot(2,3,1)
plot(meanPr_plot,meanPo_plot,'k.','markersize',15), line([-1 8],[-1 8],'color','k','linestyle',':');
[r,p]=corrcoef(meanPr_plot,meanPo_plot);
title(['r=',num2str(r(2,1)),', p=',num2str(p(2,1))])
xlabel('PC, Pre Sleep'),ylabel('PC, Post Sleep')
makepretty_DB
subplot(2,3,2)
plot(meanCo_plot,meanPo_plot,'k.','markersize',15)
[r,p]=corrcoef(meanCo_plot,meanPo_plot);
title(['r=',num2str(r(2,1)),', p=',num2str(p(2,1))])
xlabel('PC, Cond'),ylabel('PC, Post Sleep')
makepretty_DB
subplot(2,3,3)
plot(meanCo_plot,meanPo_plot-meanPr_plot,'k.','markersize',15)
[r,p]=corrcoef(meanCo_plot,meanPo_plot-meanPr_plot);
title(['r=',num2str(r(2,1)),', p=',num2str(p(2,1))])
xlabel('PC, Cond'),ylabel('PC, Post-Pre Sleep')
makepretty_DB
subplot(2,3,4)
plot(meanHab_plot,meanPo_plot-meanPr_plot,'k.','markersize',15)
[r,p]=corrcoef(meanHab_plot,meanPo_plot-meanPr_plot);
title(['r=',num2str(r(2,1)),', p=',num2str(p(2,1))])
xlabel('PC, Hab'),ylabel('PC, Post-Pre Sleep')
makepretty_DB
subplot(2,3,5)
plot(meanCo_plot-meanHab_plot,meanPo_plot-meanPr_plot,'k.','markersize',15)
[r,p]=corrcoef(meanCo_plot-meanHab_plot,meanPo_plot-meanPr_plot);
title(['r=',num2str(r(2,1)),', p=',num2str(p(2,1))])
xlabel('PC, Cond-Hab'),ylabel('PC, Post-Pre Sleep')
makepretty_DB
subplot(2,3,6)
scatter(meanPr_plot,meanPo_plot,40,meanCo_plot,'filled')
xlabel('PC, Pre Sleep'),ylabel('PC, Post Sleep')
line([-1 8],[-1 8],'color','k','linestyle',':');
makepretty_DB

%%% Bar plots
% General
fh3 = figure('color',[1 1 1], 'units', 'normalized', 'outerposition', [0 0 0.5 0.65]);
[~,h1] = PlotErrorBarN_DB([meanPr_plot meanHab_plot meanCo_plot meanPo_plot],...
    'barcolors', [1 1 1], 'barwidth', 0.6, 'newfig', 0, 'showpoints', 1);
h1.FaceColor = 'flat';
h1.FaceAlpha = .6;
h1.CData(2,:) = [0.793 .902 0.184];
h1.CData(3,:) = [0.9 0 0];
h1.CData(4,:) = [0 0 0];
set(gca,'Xtick',[1:4],'XtickLabel',{'PreSleep', 'FreeExplo', 'Learning', 'PostSleep'});
ylabel('PC score');
makepretty_DB
% Detailed conditioning
fh4 = figure('color',[1 1 1], 'units', 'normalized', 'outerposition', [0 0 0.7 0.65]);
[~,h2] = PlotErrorBarN_DB([meanPr_plot meanCo_plot meanCoFreeze_plot meanCoRipples_plot meanCoNoRipples_plot...
    meanPo_plot], 'barcolors', [1 1 1], 'barwidth', 0.6, 'newfig', 0, 'showpoints', 1);
h2.FaceColor = 'flat';
h2.FaceAlpha = .6;
h2.CData(2,:) = [0.9 0 0];
h2.CData(3,:) = [0.9 0 0];
h2.CData(4,:) = [0.9 0 0];
h2.CData(5,:) = [0.9 0 0];
h2.CData(6,:) = [0 0 0];
set(gca,'Xtick',[1:6],'XtickLabel',{'PreSleep', 'Full learning',  '\begin{tabular}{c} Learning \\ freezing\end{tabular}',...
    '\begin{tabular}{c} Learning \\ SWR\end{tabular}', '\begin{tabular}{c} Learning \\ no SWR\end{tabular}',...
    'PostSleep'}, 'TickLabelInterpreter', 'latex');
ylabel('PC score');
makepretty_DB
% Only with eignevalues of more than 1
% fh5 = figure('color',[1 1 1], 'units', 'normalized', 'outerposition', [0 0 0.9 0.65]);
% [~,h3] = PlotErrorBarN_DB([meanPr_plot(meanPr_plot>1) meanHab_plot(meanHab_plot>1) meanCo_plot(meanCo_plot>1) ...
%     meanCoFreeze_plot(meanCoFreeze_plot>1) meanCoRipples_plot(meanCoRipples_plot>1) meanCoNoRipples_plot(meanCoNoRipples_plot>1)...
%     meanPo_plot(meanPo_plot>1)], 'barcolors', [1 1 1], 'barwidth', 0.6, 'newfig', 0, 'showpoints', 1);
% h3.FaceColor = 'flat';
% h3.FaceAlpha = .6;
% h3.CData(2,:) = [0.793 .902 0.184];
% h3.CData(3,:) = [0.9 0 0];
% h3.CData(4,:) = [0.9 0 0];
% h3.CData(5,:) = [0.9 0 0];
% h3.CData(6,:) = [0.9 0 0];
% h3.CData(7,:) = [0 0 0];
% set(gca,'Xtick',[1:7],'XtickLabel',{'PreSleep', 'Free explo', 'Full learning',...
%     '\begin{tabular}{c} Learning \\ freezing\end{tabular}', '\begin{tabular}{c} Learning \\ SWR\end{tabular}',...
%     '\begin{tabular}{c} Learning \\ no SWR\end{tabular}', 'PostSleep'}, 'TickLabelInterpreter', 'latex');
% ylabel('PC score');
% makepretty_DB

%%% Location figures
% Prepare data
cnt=1;
for isession = 1:numsessions
    for itemplates = 1:numtemplates
        XYHa(:,:,cnt)=smooth2a(LocPCHab{isession}{itemplate}, 1, 1);
        XYCo(:,:,cnt)=smooth2a(LocPCCond{isession}{itemplate}, 1, 1);
        cnt=cnt+1;
    end
end
% Do a statistical test
for i=1:size(XYHa,1)
    for j=1:size(XYHa,2)
        [~,p(i,j)]=ttest(XYCo(i,j,:),XYHa(i,j,:));
        try
            [pr(i,j),hr]=ranksum(squeeze(XYCo(i,j,:)),squeeze(XYHa(i,j,:)));
        catch
            pr(i,j)=nan;
        end
    end
end
p = p';

% Mean PC Hab, Cond
fh6 = figure('color',[1 1 1], 'units', 'normalized', 'outerposition', [0 0 0.9 0.65]);
subplot(1,3,1)
imagesc(sum(XYHa,3,'omitnan')'/30)
axis xy
xlim([1.5 10.5]),ylim([1.5 10.5])
title('Mean PC score Hab')
caxis([0 1.5]), colorbar, set(gca,'XTick',[])
set(gca,'YTick',[])
makepretty_DB
subplot(1,3,2)
imagesc(sum(XYCo,3,'omitnan')'/30)
axis xy,
xlim([1.5 10.5]),ylim([1.5 10.5])
title('Mean PC score Cond')
caxis([0 1.5]), colorbar
set(gca,'XTick',[]),set(gca,'YTick',[])
colormap(hot)
makepretty_DB
% Mean PC Cond-Hab
subplot(1,3,3)
imagesc(sum(XYCo-XYHa,3,'omitnan')'/30)
axis xy
xlim([1.5 10.5]),ylim([1.5 10.5])
title('Mean PC score Cond-Hab')
colorbar
set(gca,'XTick',[]),set(gca,'YTick',[])
caxis([-0.2 0.85])


end



%% Auxiliary functions
function fh = plotIndividFigure_react_pca(tRipples, SessionEpoch, HabEpoch, CondEpoch, SWSEpoch,...
    replayGtsd, LocPCHab, LocPCCond, SWRrate, mPr, mCo, mPo, tPr, tCo, tPo, mousename, itemplate)

const = 0.3 * max(Data(replayGtsd));

% Figure
fh = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
subplot(2,4,1:3), hold on,
a1 = plot(Range(Restrict(replayGtsd, SessionEpoch.PreSleep),'s'), ...
    Data(Restrict(replayGtsd, SessionEpoch.PreSleep)),'.','color',[0.8 0.8 0.8]);
a2 = plot(Range(Restrict(replayGtsd, and(SessionEpoch.PreSleep, SWSEpoch)),'s'),...
    Data(Restrict(replayGtsd, and(SessionEpoch.PreSleep, SWSEpoch))), '.','color',[0 0 0.8]);  %% SWS Sleep
a3 = plot(Range(Restrict(replayGtsd, HabEpoch),'s'), Data(Restrict(replayGtsd, HabEpoch)),'.', 'Color', [0.793 .902 0.184]);
a4 = plot(Range(Restrict(replayGtsd, CondEpoch),'s'), Data(Restrict(replayGtsd, CondEpoch)),'r.');
plot(Range(Restrict(replayGtsd, SessionEpoch.PostSleep),'s'), ...
    Data(Restrict(replayGtsd, SessionEpoch.PostSleep)),'.','color',[0.8 0.8 0.8]);
a5 = plot(Range(Restrict(replayGtsd, and(SessionEpoch.PostSleep, SWSEpoch)),'s'),...
    Data(Restrict(replayGtsd, and(SessionEpoch.PostSleep, SWSEpoch))), '.','color', [0 0 0.2]);  %% SWS Sleep
EpochR=or(or(CondEpoch, HabEpoch),or(SessionEpoch.PostSleep,SessionEpoch.PreSleep));
a6 = plot(Range(Restrict(tRipples,EpochR),'s'), zeros(length(Restrict(tRipples,EpochR)), 1)+const, 'm.');
title('PC score over time')
legend([a1, a2, a3, a4, a5, a6],'Wake','Pre-NREM','Hab', 'Cond', 'Post-NREM', 'Ripples');

subplot(2,4,5)
PlotErrorBar4(Data(Restrict(replayGtsd, SessionEpoch.PreSleep)),Data(Restrict(replayGtsd, HabEpoch)),...
    Data(Restrict(replayGtsd, CondEpoch)), Data(Restrict(replayGtsd, SessionEpoch.PostSleep)),0,0);
set(gca,'XTick',1:4,'XTickLabel',{'PreSleep','Hab','Cond','PostSleep'});
ylabel('Mean PC score')

subplot(2,4,6), hold on,
bar(SWRrate,'k')
set(gca,'XTick',1:4,'XTickLabel',{'PreSleep','Hab','Cond','PostSleep'});
ylabel('ripples/s')
title('Ripples rate')


subplot(2,4,7), hold on,
p1 = plot(tPr,mPr,'linewidth',2, 'color', [0 0 .8]);
p2 = plot(tCo,mCo,'r','linewidth',2);
p3 = plot(tPo,mPo, 'linewidth',2, 'color', [0 0 .2]);
line([0 0],ylim,'color','k','linestyle',':')
xlim([tPr(1) tPr(end)])
xlabel('Time to ripple (s)')
ylabel('PC score')
legend([p1, p2, p3],'PreSleep','Cond','PostSleep');

subplot(2,4,4), imagesc(LocPCHab'), axis xy, ca1=caxis;
title('Hab')
subplot(2,4,8), imagesc(LocPCCond'), axis xy, ca2=caxis;
title('Cond')
colormap(hot)
if min([ca1(2),ca2(2)]) > 0
    subplot(2,4,4),caxis([0 min([ca1(2),ca2(2)])]), xlim([1.5 10.5]), ylim([1.5 10.5])
    subplot(2,4,8),caxis([0 min([ca1(2),ca2(2)])]), xlim([1.5 10.5]), ylim([1.5 10.5])
else
    subplot(2,4,4),caxis([min([ca1(2),ca2(2)]) 0]), xlim([1.5 10.5]), ylim([1.5 10.5])
    subplot(2,4,8),caxis([min([ca1(2),ca2(2)]) 0]), xlim([1.5 10.5]), ylim([1.5 10.5])
end

mtit([mousename '  Template#' num2str(itemplate)], 'FontSize', 16, 'xoff', -0.12, 'yoff', .035);

end




