function sam_hb_exp_DB()

%==========================================================================
% Details: get hb details during pre, post and cond trials
%
% INPUTS:
%       
%
% OUTPUT:
%
% NOTES: needs HeartBeatInfo.mat (from Dat2Mat_##.m [calling
%           DetectHeartBeats_EmbReact_SB.m])
%       
%       see: 
%
%
%   Written by Samuel Laventure - 29-10-2019
%   Edited by Dima Bryzgalov - 25-11-2019
%   Edited by Dima Bryzgalov - 27-05-2019
%   Completely rewritten by Dima Bryzgalov - 29-05-2019
%  
%==========================================================================

%% Parameters
% Directory to save and name of the figure to save
dir_out = '/MOBS_workingON/Dima/Ongoing_results/Heart/';
sav = false;
ntrial = 4;  %assuming that all pre, post and cond have the same number of trials 
% threshold for minimum heartbeat (get rid of outliers)
ttime = 3.5;  %time in sec before and after stim to calculate rate
stimdelay=1; %time in sec after stim onset
prebuff=500; %before stim in s1e4
postbuff=0.3*1e4; %after stim
stimdur=1200;

% Get directories
% Dir = PathForExperimentsERC_Dima('UMazePAG');
Dir = PathForExperimentsERC_DimaMAC('UMazePAG');
Dir = RestrictPathForExperiment(Dir,'Group', 'ECG');

%% PreAllocate memory
% Original data structures
behav = cell(length(Dir.path),1);
TTL = cell(length(Dir.path),1);
HB = cell(length(Dir.path),1);

% Auxiliary arrays
sessions = cell(length(Dir.path),1);
FreezeEpoch = cell(length(Dir.path),1);
NonFreezeEpoch = cell(length(Dir.path),1);
tstim = cell(length(Dir.path),ntrial);
prestartend = nan(length(Dir.path),ntrial,2);
condstartend = nan(length(Dir.path),ntrial,2);
poststartend = nan(length(Dir.path),ntrial,2);
nbstim = zeros(length(Dir.path), ntrial);
StimSent = cell(length(Dir.path),1);

% Heart rate arrays
hbpre = cell(length(Dir.path),ntrial);
hbcond = cell(length(Dir.path),ntrial);
hbpost = cell(length(Dir.path),ntrial);
hbcond_freeze = cell(length(Dir.path),ntrial);
hbcond_nonfreeze = cell(length(Dir.path),ntrial);

% Resulting data
hbpre_mean = nan(length(Dir.path),ntrial);
hbcond_mean = nan(length(Dir.path),ntrial);
hbpost_mean = nan(length(Dir.path),ntrial);
hbcond_freeze_mean = nan(length(Dir.path),ntrial);
hbcond_nonfreeze_mean = nan(length(Dir.path),ntrial);
hbpre_std = nan(length(Dir.path),ntrial);
hbcond_std = nan(length(Dir.path),ntrial);
hbpost_std = nan(length(Dir.path),ntrial);
hbcond_freeze_std = nan(length(Dir.path),ntrial);
hbcond_nonfreeze_std = nan(length(Dir.path),ntrial);

ratepre_mean = nan(length(Dir.path), ntrial, ttime*1e4);
ratepost_mean = nan(length(Dir.path), ntrial, (ttime+stimdelay)*1e4);

%% Load the data
for imouse=1:length(Dir.path)
    behav{imouse} = load([Dir.path{imouse}{1} '/behavResources.mat'], 'behavResources','SessionEpoch', 'FreezeAccEpoch');
    TTL{imouse} = load([Dir.path{imouse}{1} 'behavResources.mat'], 'TTLInfo');
    HB{imouse} = load([Dir.path{imouse}{1} '/HeartBeatInfo.mat'], 'EKG');
end

%% Prepare auxiliary arrays
for imouse=1:length(Dir.path)
    sessions{imouse} = struct2cell(behav{imouse}.SessionEpoch);  %transform the object structure into cell to allow for navigation using id
    st = Start(sessions{imouse}{1});
    en = End(sessions{imouse}{end});
    TotalEpoch = intervalSet(st,en);
    
    FreezeEpoch{imouse} = behav{imouse}.FreezeAccEpoch;
    NonFreezeEpoch{imouse} = TotalEpoch - FreezeEpoch{imouse};
    
    StimSent{imouse} = TTL{imouse}.TTLInfo.StimEpoch;
    
    % Get sessions{imouse} id
    idpre = GetSessionIdx(behav{imouse}.behavResources, 'TestPre');
    idcond = GetSessionIdx(behav{imouse}.behavResources, 'Cond');
    idpost = GetSessionIdx(behav{imouse}.behavResources, 'TestPost');
    
    % get hb per trials
    for itrial=1:ntrial
        hbpre{imouse,itrial} = Restrict(Restrict(HB{imouse}.EKG.HBRate,sessions{imouse}{idpre(itrial),1}),HB{imouse}.EKG.GoodEpoch) ;
        hbpost{imouse,itrial} = Restrict(Restrict(HB{imouse}.EKG.HBRate,sessions{imouse}{idpost(itrial),1}),HB{imouse}.EKG.GoodEpoch);
        hbcond{imouse,itrial} = Restrict(HB{imouse}.EKG.HBRate,sessions{imouse}{idcond(itrial),1});
        hbcond_freeze{imouse,itrial} = Restrict(hbcond{imouse,itrial}, FreezeEpoch{imouse});
        hbcond_nonfreeze{imouse,itrial} = Restrict(hbcond{imouse,itrial},NonFreezeEpoch{imouse});
        
        %get start time
        prestartend(imouse,itrial,1) = behav{imouse}.behavResources(idpre(itrial)).PosMat(1,1);
        prestartend(imouse,itrial,2) = behav{imouse}.behavResources(idpre(itrial)).PosMat(end,1);
        condstartend(imouse,itrial,1) = behav{imouse}.behavResources(idcond(itrial)).PosMat(1,1);
        condstartend(imouse,itrial,2) = behav{imouse}.behavResources(idcond(itrial)).PosMat(end,1);
        poststartend(imouse,itrial,1) = behav{imouse}.behavResources(idpost(itrial)).PosMat(1,1);  
        poststartend(imouse,itrial,2) = behav{imouse}.behavResources(idpost(itrial)).PosMat(end,1);
        
        % Find stim timestamps for each session and allocate memory for trial by trial data
        tstim{imouse,itrial} = Start(and(StimSent{imouse}, sessions{imouse}{idcond(itrial),1}));
        nbstim(imouse,itrial) = length(tstim{imouse,itrial});
    end
end

%% PreAllocate lfp

tstim_unroll = tstim(:);
maxmun_stim = max(sum(reshape(arrayfun(@(x) length(tstim_unroll{x}), 1:numel(tstim_unroll)),...
    size(nbstim,1), size(nbstim,2)),2));
% PreAllocate
rateprelfp = nan(length(Dir.path), ntrial,...
    maxmun_stim, ttime*1E4);
ratepostlfp = nan(length(Dir.path), ntrial,...
    maxmun_stim, (ttime+stimdelay)*1E4);
FreezeInStim = cell(length(Dir.path), ntrial, maxmun_stim);
clear tstim_unroll

numstim_PreFreeze = 0;
numstim_PreNONFreeze = 0;
%% Calculate data
for imouse=1:length(Dir.path)
    for itrial=1:ntrial
        % calculate mean and std
        [hbpre_mean(imouse,itrial), hbpre_std(imouse,itrial)] = GetNanDescStats(Data(hbpre{imouse,itrial}));
        [hbcond_mean(imouse,itrial), hbcond_std(imouse,itrial)] = GetNanDescStats(Data(hbpre{imouse,itrial}));
        [hbpost_mean(imouse,itrial), hbpost_std(imouse,itrial)] = GetNanDescStats(Data(hbpre{imouse,itrial}));
        [hbcond_freeze_mean(imouse,itrial), hbcond_freeze_std(imouse,itrial)] = GetNanDescStats(Data(hbpre{imouse,itrial}));
        [hbcond_nonfreeze_mean(imouse,itrial), hbcond_nonfreeze_std(imouse,itrial)] = GetNanDescStats(Data(hbpre{imouse,itrial}));
              
        % stim react
        HBTimes = Data(HB{imouse}.EKG.HBTimes);
        ratepre = cell(length(tstim{imouse,itrial}),1);
        ratepost = cell(length(tstim{imouse,itrial}),1);
        
        isfreeze_beforestim = false(length(tstim{imouse,itrial}),1);
        for istim=1:length(tstim{imouse,itrial})
            FreezeInStim{imouse, itrial, istim} = and(and(FreezeEpoch{imouse},sessions{imouse}{idcond(itrial)}),...
                intervalSet(tstim{imouse,itrial}(istim)-ttime*1e4, tstim{imouse,itrial}(istim)+(ttime+stimdelay)*1e4));
            FreezeInStim{imouse, itrial, istim} = intervalSet(Start(FreezeInStim{imouse, itrial, istim})-tstim{imouse,itrial}(istim)+ttime*1e4,...
                End(FreezeInStim{imouse, itrial, istim})-tstim{imouse,itrial}(istim)+ttime*1e4);
            
            % Find out whether the mouse was freezing before stimulus
            UnCorrIStart = Start(FreezeInStim{imouse, itrial, istim}) - ttime*1e4;
            if sum(UnCorrIStart<0) > 0
                isfreeze_beforestim(istim) = true;
                numstim_PreFreeze = numstim_PreFreeze + 1;
            else
                numstim_PreNONFreeze = numstim_PreNONFreeze + 1;
            end
            
            % HB pre stim
            ratepre_id = find((HBTimes<tstim{imouse,itrial}(istim)) ...
                & (HBTimes>tstim{imouse,itrial}(istim)-ttime*1E4));
            ratepre{istim} = HBTimes(ratepre_id);
            hr_pre = 1./movmedian(diff(ratepre{istim}/1e4),[3 0]);
            delta_pre = diff(ratepre{istim}/1e4);
            
            % setting up rate on LFP points to allow for average between
            % stim
            fst = 1;
            for ihb=1:length(ratepre_id)
                if ihb>1
                    rateprelfp(imouse, itrial,istim,fst:single(ceil(ratepre{istim}(ihb)-ratepre{istim}(1)))) = hr_pre(ihb-1);
                    deltaprelfp(istim,fst:single(ceil(ratepre{istim}(ihb)-ratepre{istim}(1)))) = delta_pre(ihb-1);
                    fst = (ratepre{istim}(ihb)+1)-ratepre{istim}(1);
                end
            end
            
            % HB post stim
            ratepost_id = find((HBTimes>tstim{imouse,itrial}(istim)) ...
                & (HBTimes<tstim{imouse,itrial}(istim)+(stimdelay+ttime)*1E4));
            ratepost{istim} = HBTimes(ratepost_id);
            hr_post = 1./movmedian(diff(ratepost{istim}/1e4),[0 3]);
            delta_post = diff(ratepost{istim}/1e4);
            
            % setting up rate on LFP points to allow for average between
            % stim
            fst = (stimdelay+ttime)*1E4;
            tstart = tstim{imouse,itrial}(istim);
            for ihb=length(ratepost_id):-1:1
                if ihb>1
                    if ratepost{istim}(ihb)-tstim{imouse,itrial}(istim) > 1500
                        ratepostlfp(imouse, itrial, istim,floor(ratepost{istim}(ihb)-tstart):fst) = hr_post(ihb-1);
                        deltapostlfp(istim,floor(ratepost{istim}(ihb)-tstart):fst) = delta_post(ihb-1);
                        fst = floor(ratepost{istim}(ihb)-tstart)-1;
                    end
                end
            end
        end
        
        % calculate mean per trial
        if nbstim(imouse,itrial) >= 1
            ratepre_mean(imouse,itrial,1:ttime*1E4) = nanmean(squeeze(rateprelfp(imouse,itrial,:,:)));
            ratepost_mean(imouse,itrial,1:(ttime+stimdelay)*1E4) = nanmean(squeeze(ratepostlfp(imouse,itrial,:,:)));
            
            ratepre_PreFreeze_mean(imouse,itrial,1:ttime*1E4) = nanmean(squeeze(rateprelfp(imouse,itrial,isfreeze_beforestim,:)));
            ratepre_PreNONFreeze_mean(imouse,itrial,1:ttime*1E4) = nanmean(squeeze(rateprelfp(imouse,itrial,~isfreeze_beforestim,:)));
            ratepost_PreFreeze_mean(imouse,itrial,1:(ttime+stimdelay)*1E4) = nanmean(squeeze(ratepostlfp(imouse,itrial,isfreeze_beforestim,:)));
            ratepost_PreNONFreeze_mean(imouse,itrial,1:(ttime+stimdelay)*1E4) = nanmean(squeeze(ratepostlfp(imouse,itrial,~isfreeze_beforestim,:)));
        else
            ratepre_mean(imouse,itrial,1:ttime*1E4) = nan(1,ttime*1E4);
            ratepost_mean(imouse,itrial,1:(ttime+stimdelay)*1E4) = nan(1,(ttime+stimdelay)*1E4);
            
            ratepre_PreFreeze_mean(imouse,itrial,1:ttime*1E4) = nan(1,ttime*1E4);
            ratepre_PreNONFreeze_mean(imouse,itrial,1:ttime*1E4) = nan(1,ttime*1E4);
            ratepost_PreFreeze_mean(imouse,itrial,1:(ttime+stimdelay)*1E4) = nan(1,(ttime+stimdelay)*1E4);
            ratepost_PreNONFreeze_mean(imouse,itrial,1:(ttime+stimdelay)*1E4) = nan(1,(ttime+stimdelay)*1E4);
        end
    end
end
     

% Prepare data for figures
% Overall means
hbpre_meanall = squeeze(mean(hbpre_mean,2));
hbpost_meanall = squeeze(mean(hbpost_mean,2));
hbcond_meanall = squeeze(mean(hbcond_mean,2));
hbcond_freeze_meanall = squeeze(mean(hbcond_freeze_mean,2));
hbcond_nonfreeze_meanall = squeeze(mean(hbcond_nonfreeze_mean,2));

ratepre_all_mean = squeeze(nanmean(nanmean(ratepre_mean,2)));
ratepre_all_std = squeeze(nanstd(nanmean(ratepre_mean,2)));
ratepost_all_mean = squeeze(nanmean(nanmean(ratepost_mean,2)));
ratepost_all_std = squeeze(nanstd(nanmean(ratepost_mean,2)));

ratepre_PreFreeze = squeeze(nanmean(nanmean(ratepre_PreFreeze_mean,2)));
ratepre_PreNONFreeze = squeeze(nanmean(nanmean(ratepre_PreNONFreeze_mean,2)));
ratepost_PreFreeze = squeeze(nanmean(nanmean(ratepost_PreFreeze_mean,2)));
ratepost_PreNONFreeze = squeeze(nanmean(nanmean(ratepost_PreNONFreeze_mean,2)));
ratepre_PreFreeze_std = squeeze(nanstd(nanmean(ratepre_PreFreeze_mean,2)));
ratepre_PreNONFreeze_std = squeeze(nanstd(nanmean(ratepre_PreNONFreeze_mean,2)));
ratepost_PreFreeze_std = squeeze(nanstd(nanmean(ratepost_PreFreeze_mean,2)));
ratepost_PreNONFreeze_std = squeeze(nanstd(nanmean(ratepost_PreNONFreeze_mean,2)));

ratepre_trial_mean = nanmean(nanmean(ratepre_mean(:,:,ttime*1E4-19999:ttime*1E4),3),2);
ratepost_trial_mean = nanmean(nanmean(ratepost_mean(:,:,stimdur+1:stimdur+20000),3),2);

nbstim_all = sum(sum(nbstim)); 

%% 
%#####################################################################
%#
%#                        F I G U R E S
%#
%#####################################################################

%% Single stim mouse trial check-up

% DO ANIMATION
nbstim_mouse = sum(nbstim,2);
for imouse=1:length(Dir.path)
    f = figure('units', 'normalized', 'outerposition', [0 0 0.5 0.5], 'Name', Dir.name{imouse});
    title(Dir.name{imouse});
    tabs = arrayfun(@(itab) uitab('Title', ['Trial' num2str(itab)]), 1:nbstim_mouse(imouse));
    ax = arrayfun(@(tab) axes(tab, 'NextPlot', 'add'), tabs);
    
    axnum = 0;
    for itrial=1:ntrial
        for istim = 1:nbstim(imouse,itrial)
            axnum = axnum+1;
            rectangle(ax(axnum),'Position',[ttime*1E4-prebuff,0,prebuff+postbuff,20],'FaceColor',[.9 .9 .9],'EdgeColor','none',...
                'LineWidth',.01)           
            rectangle(ax(axnum),'Position',[ttime*1E4,0,stimdur,20],'FaceColor',[.6 .6 .6],'EdgeColor','none',...
                'LineWidth',.01)
            
            plot(ax(axnum),[5001:ttime*1E4-prebuff],squeeze(rateprelfp(imouse,itrial,istim,5001:end-prebuff)),...
                'LineWidth',2) % starts 500 ms after because median is culculated on trailing data (inverse for post)
            plot(ax(axnum),[ttime*1E4+postbuff:(ttime*2+stimdelay)*1e4-5000],...
                squeeze(squeeze(squeeze(ratepostlfp(imouse,itrial,istim,postbuff:end-5000)))),...
                'LineWidth',2)
            
            freezeauxdat = tsd([5001:(ttime*2+stimdelay)*1e4-5000],ones(length([5001:(ttime*2+stimdelay)*1e4-5000]),1));
            for ii = 1:length(Start(FreezeInStim{imouse,itrial,istim}))
                plot(ax(axnum),Range(Restrict(freezeauxdat,subset(FreezeInStim{imouse,itrial,istim},ii))),...
                    Data(Restrict(freezeauxdat,subset(FreezeInStim{imouse,itrial,istim},ii)))*13.8,'c','linewidth',5)
            end
            
            ylim(ax(axnum),[6 14.5])
            xlim(ax(axnum),[5000 (ttime*2+stimdelay)*1e4-5000])
            set(ax(axnum), 'Xtick', 5001:10000:((ttime*2+stimdelay)*1e4)-5000,...
                'Xticklabel', num2cell([-1*(ttime-.5):1:(stimdelay+ttime-.5)]), 'FontSize', 16, 'FontWeight', 'bold')
%             ax(axnum).Layer = 'top';
            
            xlabel(ax(axnum),'time (s)')
            ylabel(ax(axnum),'Hz')
            title(ax(axnum),[Dir.name{imouse} ' Tr.' num2str(itrial) ' Stim.' num2str(istim)], 'FontSize', 14)
            
        end
    end
end

%%  Overall separation PREFREEZE
supertit = 'Heart Rate Dynamics locked to stim';
f1 = figure('Color',[1 1 1], 'units', 'normalized', 'outerposition', [0 0 0.605 0.63],'Name', supertit, 'NumberTitle','off')

rectangle('Position',[ttime*1E4-prebuff,0,prebuff+postbuff,20],'FaceColor',[.9 .9 .9],'EdgeColor','none',...
    'LineWidth',.01)
hold on

rectangle('Position',[ttime*1E4,0,stimdur,20],'FaceColor',[.6 .6 .6],'EdgeColor','none',...
    'LineWidth',.01)
hold on

h_freeze = shadedErrorBar([5001:ttime*1E4-prebuff],ratepre_PreFreeze(5001:end-prebuff),ratepre_PreFreeze_std(5001:end-prebuff),...
    {'LineWidth', 2, 'Color', 'b'},1);
h_nonfreeze = shadedErrorBar([5001:ttime*1E4-prebuff],ratepre_PreNONFreeze(5001:end-prebuff),ratepre_PreNONFreeze_std(5001:end-prebuff),...
    {'LineWidth', 2, 'Color', 'r'},1);
hold on
shadedErrorBar([ttime*1E4+postbuff:(ttime*2+stimdelay)*1e4-5000],ratepost_PreFreeze(postbuff:end-5000),...
    ratepost_PreFreeze_std(postbuff:end-5000),{'LineWidth', 2,'Color', 'b'},1);
shadedErrorBar([ttime*1E4+postbuff:(ttime*2+stimdelay)*1e4-5000],ratepost_PreNONFreeze(postbuff:end-5000),...
    ratepost_PreNONFreeze_std(postbuff:end-5000),{'LineWidth', 2,'Color', 'r'},1);

ylim([9 14])
xlim([5000 (ttime*2+stimdelay)*1e4-5000])
set(gca, 'Xtick', 5001:5000:((ttime*2+stimdelay)*1e4),...
    'Xticklabel', num2cell([-1*(ttime-.5):.5:(stimdelay+ttime)]),...
    'FontSize', 16, 'FontWeight', 'bold')
ax = gca;
ax.Layer = 'top';
hold off
xlabel('time (s)', 'FontWeight', 'bold')
ylabel('Hz')
title('Heart Rate Dynamics locked to PAG stim')
legend([h_freeze.mainLine h_nonfreeze.mainLine],{['Stim when freezing (N=' num2str(numstim_PreFreeze) ')'],...
    ['Stim when moving (N=' num2str(numstim_PreNONFreeze) ')']}, 'Location', 'SouthEast')
if sav
    saveFigure(f1,'Mean_PAGHeart_FreezeSep',[dropbox dir_out]);
    saveas(f1,[dropbox dir_out 'Mean_PAGHeart_FreezeSep.fig']);
end

%% HR freezing/nonfrezing
Pl = {hbcond_nonfreeze_meanall, hbcond_freeze_meanall};
Cols = {[0.9 0.5 0.5], [0.2 0.2 0.9]};
xlabs = {'Moving','Freezing'};
tit = ['Average heart rate (N=' num2str(length(Dir.path)) ')'];
f3 = figure('units', 'normalized', 'outerposition', [0 0 0.4 0.55]);
MakeBoxPlot_DB(Pl,Cols,1:2,[],1);
set(gca,'XTick', [1:2], 'XTickLabel', xlabs, 'FontName', 'Helvetica',...
    'FontSize', 16, 'FontWeight', 'bold');
ylabel('Heart rate (Hz)', 'FontName', 'Helvetica', 'FontSize', 16, 'FontWeight', 'bold');
title(tit, 'FontName', 'Helvetica', 'FontSize', 16, 'FontWeight', 'bold');
ylim([9 14])
[p,h5,stats] = signrank(hbcond_nonfreeze_meanall, hbcond_freeze_meanall);
if p < 0.05
    sigstar_DB({{1,2}},p,0, 'StarSize',14);
end

%% HR with stim
% delta beat hotmaps
figure, imagesc(sort(deltaprelfp,1,'descend'))
    title('Pre-stim')
    ylabel('stims')  
    xlabel('time (s)')
    
figure, imagesc(sort(deltapostlfp,1,'ascend'))
    title('Post-stim')
    ylabel('stims')  
    xlabel('time (s)')


% per trial per mouse
supertit = 'Heart Rate Dynamics locked to stim';
figure('Color',[1 1 1], 'rend','painters','pos',[1 1 2100 400],'Name', supertit, 'NumberTitle','off', 'WindowState',...
    'maximized')
    for imouse=1:length(Dir.path)
        for itrial=1:ntrial
            subplot(length(Dir.path),ntrial,itrial+((imouse-1)*ntrial))
                rectangle('Position',[ttime*1E4-prebuff,0,prebuff+postbuff,20],'FaceColor',[.9 .9 .9],'EdgeColor','none',...
                'LineWidth',.01)
                hold on

                rectangle('Position',[ttime*1E4,0,stimdur,20],'FaceColor',[.6 .6 .6],'EdgeColor','none',...
                'LineWidth',.01)

                hold on
                plot([5001:ttime*1E4-prebuff],squeeze(squeeze(ratepre_mean(imouse,itrial,5001:end-prebuff))),...
                    'LineWidth',2) % starts 500 ms after because median is culculated on trailing data (inverse for post)
                hold on
                plot([ttime*1E4+postbuff:(ttime*2+stimdelay)*1e4-5000], squeeze(squeeze(ratepost_mean(imouse,itrial,postbuff:end-5000))),...
                    'LineWidth',2)
                                
                ylim([10 14.5])
                xlim([5000 (ttime*2+stimdelay)*1e4-5000])
                set(gca, 'Xtick', 5001:10000:((ttime*2+stimdelay)*1e4)-5000,...
                    'Xticklabel', num2cell([-1*(ttime-.5):1:(stimdelay+ttime-.5)]), 'FontSize', 16, 'FontWeight', 'bold')
                ax = gca;
                ax.Layer = 'top';

                hold off
                if imouse == length(Dir.path)
                    xlabel('time (s)')
                end
                ylabel('Hz')
                title([Dir.name{imouse} ' Tr.' num2str(itrial) ' ' num2str(nbstim(imouse,itrial)) ' stims'], 'FontSize', 14)
        end
    end
    % Save figure
    if sav
        print([dir_out 'heartrate_stim'], '-dpng', '-r300');
    end

% Overall
supertit = 'Heart Rate Dynamics locked to stim';
figure('Color',[1 1 1], 'rend','painters','pos',[1 1 1200 800],'Name', supertit, 'NumberTitle','off')

        rectangle('Position',[ttime*1E4-prebuff,0,prebuff+postbuff,20],'FaceColor',[.9 .9 .9],'EdgeColor','none',...
        'LineWidth',.01)
        hold on
        
        rectangle('Position',[ttime*1E4,0,stimdur,20],'FaceColor',[.6 .6 .6],'EdgeColor','none',...
        'LineWidth',.01)
        hold on

%         plot([5001:ttime*1E4-prebuff],ratepre_all_mean(5001:end-prebuff))
        shadedErrorBar([5001:ttime*1E4-prebuff],ratepre_all_mean(5001:end-prebuff),ratepre_all_std(5001:end-prebuff),...
        {'LineWidth', 2, 'Color', 'r'});
        hold on
        shadedErrorBar([ttime*1E4+postbuff:(ttime*2+stimdelay)*1e4-5000],ratepost_all_mean(postbuff:end-5000),...
            ratepost_all_std(postbuff:end-5000),{'LineWidth', 2,'Color', 'b'});
%         plot([ttime*1E4+postbuff:(ttime*2+stimdelay)*1e4-5000], ratepost_all_mean(postbuff:end-5000))

        ylim([10 14])
                xlim([5000 (ttime*2+stimdelay)*1e4-5000])
                set(gca, 'Xtick', 5001:5000:((ttime*2+stimdelay)*1e4),...
                    'Xticklabel', num2cell([-1*(ttime-.5):.5:(stimdelay+ttime)]),...
                    'FontSize', 16, 'FontWeight', 'bold')
        ax = gca;
        ax.Layer = 'top';
        t=annotation('textbox',[.8 .7 .2 .1],'String',['n=' num2str(nbstim_all) ' stim'],'FitBoxToText','on',...
            'FontWeight','bold');
        sz = t.FontSize;
        t.FontSize = 14;
        t.EdgeColor =  'none';
        hold off
        xlabel('time (s)', 'FontWeight', 'bold')
        ylabel('Hz')    
        title('Heart Rate Dynamics locked to PAG stim')
    % Save figure
    if sav
        print([dir_out 'heartrate_stim_all'], '-dpng', '-r300');
    end

% Overall - ALL data points
supertit = 'Heart Rate Dynamics locked to stim - no buffer';
figure('Color',[1 1 1], 'rend','painters','pos',[1 1 600 400],'Name', supertit, 'NumberTitle','off')
        rectangle('Position',[ttime*1E4,0,stimdur,20],'FaceColor',[.6 .6 .6],'EdgeColor','none',...
        'LineWidth',.01)
        hold on

        plot([5001:ttime*1E4],ratepre_all_mean(5001:end))
        hold on
        plot([ttime*1E4+1:(ttime*2+stimdelay)*1e4-5000], ratepost_all_mean(1:end-5000))
        hold on
        lx=(ttime+.075*2)*1e4+1200;
        plot([lx lx], [10 13]); 
        hold on
        lx=(ttime+.075)*1e4+1200;
        plot([lx lx], [10 13]); 
        ylim([10 13])
                xlim([5000 (ttime*2+stimdelay)*1e4-5000])
                set(gca, 'Xtick', 5001:5000:((ttime*2+stimdelay)*1e4)-5000,...
                    'Xticklabel', num2cell([-1*(ttime-.5):.5:(stimdelay+ttime-.5)]))
        ax = gca;
        ax.Layer = 'top';
        t=annotation('textbox',[.8 .7 .2 .1],'String',['n=' num2str(nbstim_all) ' stim'],'FitBoxToText','on');
        sz = t.FontSize;
        t.FontSize = 10;
        t.EdgeColor =  'none';
        hold off
        xlabel('time (s)')
        ylabel('Hz')    
        title('Heart Rate Dynamics locked to stim')
    % Save figure
    if sav
        print([dir_out 'heartrate_stim_all_nobuffer'], '-dpng', '-r300');
    end
    
    
supertit = 'Heart rate difference pre/post stim (+/-2s) ';
figure('Color',[1 1 1], 'rend','painters','pos',[1 1 400 400],'Name', supertit, 'NumberTitle','off')
    
        [p,h, her] = PlotErrorBarN_SL([ratepre_trial_mean ratepost_trial_mean],...
            'colorpoints',1,'barcolors', [.3 .3 .3], 'barwidth', 0.6, 'newfig', 0);
         ylim([10.5 13.5]);
        set(gca,'Xtick',[1:2],'XtickLabel',{'Pre-stim','Post-stim'});
        set(gca, 'FontSize', 12);
        h.FaceColor = 'flat';
        h.CData(1,:) = [.3 .3 .3];
        h.CData(2,:) = [.6 .6 .6];
        
        set(h, 'LineWidth', 2);
        set(her, 'LineWidth', 2);
        ylabel('Hz');
        title('Average HR difference pre/post stim', 'FontSize', 16);
        
    % Save figure
    if sav
        print([dir_out 'average_hr_diff_2sec'], '-dpng', '-r300');
    end

%% HB dynamics

ymax=max([Data(hbpre{imouse,itrial}); Data(hbcond{imouse,itrial}); Data(hbpost{imouse,itrial})]);
ymin=min([Data(hbpre{imouse,itrial}); Data(hbcond{imouse,itrial}); Data(hbpost{imouse,itrial})]);
supertit = 'Heart Rate Dynamics per trial - Pre-Tests';
figure('Color',[1 1 1], 'rend','painters','pos',[1 1 2100 400],'Name', supertit, 'NumberTitle','off')
    for imouse=1:length(Dir.path)
        for itrial=1:ntrial
            subplot(length(Dir.path),ntrial,itrial+((imouse-1)*ntrial))
                plot((Range(hbpre{imouse,itrial})./1E4)-prestartend(imouse,itrial,1),Data(hbpre{imouse,itrial}))
                ylim([ymin-ymin*.15 ymax+ymax*.15])
                xlim([0 prestartend(imouse,itrial,2)-prestartend(imouse,itrial,1)])
                set(gca, 'Xtick', 0:20:prestartend(imouse,itrial,2)-prestartend(imouse,itrial,1),...
                    'Xticklabel', num2cell([0:20:prestartend(imouse,itrial,2)-prestartend(imouse,itrial,1)])) 
                if itrial==1
                    xlabel('time')
                    ylabel('Hz')    
                    title([Dir.name{imouse}])
                end
        end
    end
    % Save figure
    if sav
        print([dir_out 'heartrate_dyna_pre'], '-dpng', '-r300');
    end

supertit = 'Heart Rate Dynamics per trial - Cond';
figure('Color',[1 1 1], 'rend','painters','pos',[1 1 2100 400],'Name', supertit, 'NumberTitle','off')
    for imouse=1:length(Dir.path)
        for itrial=1:ntrial
            subplot(length(Dir.path),ntrial,itrial+((imouse-1)*ntrial))
                plot((Range(hbcond{imouse,itrial})./1E4)-condstartend(imouse,itrial,1),Data(hbcond{imouse,itrial}))
                hold on
                ylim([ymin-ymin*.15 ymax+ymax*.15])
                %plotting stim markers
                %   fucking long line that could be summarized as this:
                %   plot(PosMat(PosMat:,4)==1,1)-time_of_start,PosMat(PosMat:,4)==1,1)*0+a
                %   bit more than the max on the y axis);
                plot(behav{imouse}.behavResources(idcond(itrial)).PosMat(behav{imouse}.behavResources(idcond(itrial)).PosMat(:,4)==1,1)-condstartend(imouse,itrial,1)...
                    ,behav{imouse}.behavResources(idcond(itrial)).PosMat(behav{imouse}.behavResources(idcond(itrial)).PosMat(:,4)==1,1)*0+max(ylim)*.95,...
                    'g*')
                xlim([0 condstartend(imouse,itrial,2)-condstartend(imouse,itrial,1)])
                set(gca, 'Xtick', 0:80:condstartend(imouse,itrial,2)-condstartend(imouse,itrial,1),...
                    'Xticklabel', num2cell([0:80:condstartend(imouse,itrial,2)-condstartend(imouse,itrial,1)])) 
                hold off
                if itrial==1
                    xlabel('time')
                    ylabel('Hz')    
                    title([Dir.name{imouse}])
                end
        end
    end
    % Save figure
    if sav
        print([dir_out 'heartrate_dyna_cond'], '-dpng', '-r300');
    end
    
supertit = 'Heart Rate Dynamics per trial - Post-Tests';
figure('Color',[1 1 1], 'rend','painters','pos',[1 1 2100 400],'Name', supertit, 'NumberTitle','off')
    for imouse=1:length(Dir.path)
        for itrial=1:ntrial
            subplot(length(Dir.path),ntrial,itrial+((imouse-1)*ntrial))
                plot((Range(hbpost{imouse,itrial})./1E4)-poststartend(imouse,itrial,1),Data(hbpost{imouse,itrial}))
                ylim([ymin-ymin*.15 ymax+ymax*.15])
                xlim([0 poststartend(imouse,itrial,2)-poststartend(imouse,itrial,1)])
                set(gca, 'Xtick', 0:20:poststartend(imouse,itrial,2)-poststartend(imouse,itrial,1),...
                    'Xticklabel', num2cell([0:20:poststartend(imouse,itrial,2)-poststartend(imouse,itrial,1)])) 
                if itrial==1
                    xlabel('time')
                    ylabel('Hz')    
                    title([Dir.name{imouse}])
                end
        end
    end
    % Save figure
    if sav
        print([dir_out 'heartrate_dyna_post'], '-dpng', '-r300');
    end
    
    %% HB dist
    supertit = 'Heart Rate Histograms per trial - Pre-Tests';
    figure('Color',[1 1 1], 'rend','painters','pos',[1 1 900 400],'Name', supertit, 'NumberTitle','off')
    for imouse=1:length(Dir.path)
        for itrial=1:ntrial
            subplot(length(Dir.path),ntrial,itrial+((imouse-1)*ntrial))
            if ~isempty(Data(hbpre{imouse,itrial}))
                hist(Data(hbpre{imouse,itrial}),length(Data(hbpre{imouse,itrial})))
                ylim([0 225])
                xlim([8 14.5])
                if itrial==1
                    xlabel('time')
                    ylabel('Hz')
                    title([Dir.name{imouse}])
                end
            end
        end
    end
    % Save figure
    if sav
        print([dir_out 'heartrate_hist_pre'], '-dpng', '-r300');
    end
    
    supertit = 'Heart Rate Histograms per trial - Cond';
    figure('Color',[1 1 1], 'rend','painters','pos',[1 1 900 400],'Name', supertit, 'NumberTitle','off')
    for imouse=1:length(Dir.path)
        for itrial=1:ntrial
            subplot(length(Dir.path),ntrial,itrial+((imouse-1)*ntrial))
            if ~isempty(Data(hbpre{imouse,itrial}))
                hist(Data(hbcond{imouse,itrial}),length(Data(hbcond{imouse,itrial})))
                ylim([0 225])
                xlim([8 14.5])
                if itrial==1
                    xlabel('Hz')
                    ylabel('Number')
                    title([Dir.name{imouse}])
                end
            end
        end
    end
    % Save figure
    if sav
        print([dir_out 'heartrate_hist_cond'], '-dpng', '-r300');
    end
    
    supertit = 'Heart Rate Histograms per trial - Post-Tests';
    figure('Color',[1 1 1], 'rend','painters','pos',[1 1 900 400],'Name', supertit, 'NumberTitle','off')
    for imouse=1:length(Dir.path)
        for itrial=1:ntrial
            subplot(length(Dir.path),ntrial,itrial+((imouse-1)*ntrial))
            if ~isempty(Data(hbpre{imouse,itrial}))
                hist(Data(hbpost{imouse,itrial}),length(Data(hbpost{imouse,itrial})))
                ylim([0 225])
                xlim([8 14.5])
                if itrial==1
                    xlabel('Hz')
                    ylabel('Number')
                    title([Dir.name{imouse}])
                end
            end
        end
    end
    % Save figure
    if sav
        print([dir_out 'heartrate_hist_post'], '-dpng', '-r300');
    end


%% averaged HB
maxy=14;
supertit = 'Heart Rate';
figure('Color',[1 1 1], 'rend','painters','pos',[1 1 900 400],'Name', supertit, 'NumberTitle','off')
    subplot(2,3,1:3)
        [p,h, her] = PlotErrorBarN_SL([hbpre_meanall hbcond_meanall hbpost_meanall],...
            'colorpoints',1,'barcolors', [.3 .3 .3], 'barwidth', 0.6, 'newfig', 0);
        ylim([11 maxy]);
        set(gca,'Xtick',[1:3],'XtickLabel',{'Pre', 'Cond', 'Post'});
        set(gca, 'FontSize', 12);
        h.FaceColor = 'flat';
        h.CData(1,:) = [.3 .3 .3];
        h.CData(2,:) = [.6 .6 .6];
        h.CData(3,:) = [1 1 1];
        set(h, 'LineWidth', 2);
        set(her, 'LineWidth', 2);
        ylabel('Beats/sec');
        title('Heart Rate across sessions{imouse}', 'FontSize', 16);

    subplot(2,3,4)
        [p,h, her] = PlotErrorBarN_SL(hbpre_mean,...
            'colorpoints',1,'barcolors', [.3 .3 .3], 'barwidth', 0.6, 'newfig', 0);
        ylim([11 maxy]);
        set(gca,'Xtick',[1:4]);
        set(gca, 'FontSize', 12);
        set(h, 'LineWidth', 2);
        set(her, 'LineWidth', 2);
        ylabel('Beats/sec');
        xlabel('Trials');
        title('Pre', 'FontSize', 16);   
            
    subplot(2,3,5)
        [p,h, her] = PlotErrorBarN_SL(hbcond_mean,...
            'colorpoints',1,'barcolors', [.6 .6 .6], 'barwidth', 0.6, 'newfig', 0);
        ylim([11 maxy]);
        set(gca,'Xtick',[1:4]);
        set(gca, 'FontSize', 12);
        set(h, 'LineWidth', 2);
        set(her, 'LineWidth', 2);
        ylabel('Beats/sec');
        xlabel('Trials');
        title('Cond', 'FontSize', 16); 
        
    subplot(2,3,6)
        [p,h, her] = PlotErrorBarN_SL(hbpost_mean,...
            'colorpoints',1,'barcolors', [1 1 1], 'barwidth', 0.6, 'newfig', 0);
        ylim([11 maxy]);
        set(gca,'Xtick',[1:4]);
        set(gca, 'FontSize', 12);
        set(h, 'LineWidth', 2);
        set(her, 'LineWidth', 2);
        ylabel('Beats/sec');
        xlabel('Trials');
        title('Post', 'FontSize', 16);         

    % Save figure
    if sav
        print([dir_out 'heartrate_sessions{imouse}'], '-dpng', '-r300');
    end

    

