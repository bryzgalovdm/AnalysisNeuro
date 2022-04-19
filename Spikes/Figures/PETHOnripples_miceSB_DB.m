%% Parameters
IsSave=0;
mice = [490 507 508 509 510 512 514];
wi = 0.5; % in sec
binsize = 0.01; % in sec
nbins = 101;
titles = {'PreSleep', 'Conditioning', 'PostSleep'};

%% Get data
% PreSleep
Dir_temp = PathForExperimentsEmbReact('SleepPreUMaze');
n=1;
for mm = 1:length(Dir_temp.ExpeInfo)
    if ismember(Dir_temp.ExpeInfo{mm}{1}.nmouse,mice)
        Dir_SleepPre.path{n} = Dir_temp.path{mm};
        Dir_SleepPre.ExpeInfo{n} = Dir_temp.ExpeInfo{mm};
        n=n+1;
    end
end
for idir = 1:length(Dir_SleepPre.path)
    s_pre{idir} = ConcatenateDataFromFolders_SB(Dir_SleepPre.path{idir},'Spikes');
    r_pre{idir} = ConcatenateDataFromFolders_SB(Dir_SleepPre.path{idir},'Ripples');
end

% Cond
Dir_temp = PathForExperimentsEmbReact('UMazeCond');
n=1;
for mm = 1:length(Dir_temp.ExpeInfo)
    if ismember(Dir_temp.ExpeInfo{mm}{1}.nmouse,mice)
        Dir_Cond.path{n} = Dir_temp.path{mm};
        Dir_Cond.ExpeInfo{n} = Dir_temp.ExpeInfo{mm};
        n=n+1;
    end
end
for idir = 1:length(Dir_Cond.path)
    s_cond{idir} = ConcatenateDataFromFolders_SB(Dir_Cond.path{idir},'Spikes');
    r_cond{idir} = ConcatenateDataFromFolders_SB(Dir_Cond.path{idir},'Ripples');
end

% PostSleep
Dir_temp = PathForExperimentsEmbReact('SleepPostUMaze');
n=1;
for mm = 1:length(Dir_temp.ExpeInfo)
    if ismember(Dir_temp.ExpeInfo{mm}{1}.nmouse,mice)
        Dir_SleepPost.path{n} = Dir_temp.path{mm};
        Dir_SleepPost.ExpeInfo{n} = Dir_temp.ExpeInfo{mm};
        n=n+1;
    end
end
for idir = 1:length(Dir_SleepPost.path)
    s_post{idir} = ConcatenateDataFromFolders_SB(Dir_SleepPost.path{idir},'Spikes');
    r_post{idir} = ConcatenateDataFromFolders_SB(Dir_SleepPost.path{idir},'Ripples');
end
%% Calculate PETH
%PreSleep
cnt=0;
for imouse = 1:length(Dir_SleepPre.path)
    for num = 1:length(s_pre{imouse})
        if ~isempty(s_pre{imouse})
            cnt=cnt+1;
            pre(cnt,:) = PETH_KJ(Data(s_pre{imouse}{num}), Range(r_pre{imouse}), binsize*1e4, nbins);
            pre_norip(cnt,:) = PETH_KJ(Data(s_pre{imouse}{num}), Range(r_pre{imouse})+0.55*1e4, binsize*1e4, nbins);
            
        end
    end
end
% Cond
cnt=0;
for imouse = 1:length(Dir_Cond.path)
    for num = 1:length(s_cond{imouse})
        if ~isempty(s_cond{imouse})
            cnt=cnt+1;
            cond(cnt,:) = PETH_KJ(Data(s_cond{imouse}{num}), Range(r_cond{imouse}), binsize*1e4, nbins);
            cond_norip(cnt,:) = PETH_KJ(Data(s_cond{imouse}{num}), Range(r_cond{imouse})+0.55*1e4, binsize*1e4, nbins);
            
        end
    end
end
% PostSleep
cnt=0;
for imouse = 1:length(Dir_SleepPost.path)
    for num = 1:length(s_post{imouse})
        if ~isempty(s_post{imouse})
            cnt=cnt+1;
            post(cnt,:) = PETH_KJ(Data(s_post{imouse}{num}), Range(r_post{imouse}), binsize*1e4, nbins);
            post_norip(cnt,:) = PETH_KJ(Data(s_post{imouse}{num}), Range(r_post{imouse})+0.55*1e4, binsize*1e4, nbins);
            
        end
    end
end

%% Prepare data to plot
% Sort conditioning
temp = zscore(cond')';
SustVal = nanmean(temp(:,40:65),2);
[~,ind1] = sort(SustVal);
cond_toplot = temp(ind1,:);
temp = zscore(cond_norip')';
cond_norip_toplot = temp(ind1,:);
% Sort PreSleep
temp = zscore(pre')';
pre_toplot = temp(ind1,:);
temp = zscore(pre_norip')';
pre_norip_toplot = temp(ind1,:);
% Sort PostSleep
temp = zscore(post')';
post_toplot = temp(ind1,:);
temp = zscore(post_norip')';
post_norip_toplot = temp(ind1,:);

% Organize
onRipples = {pre_toplot, cond_toplot, post_toplot};
noRipples = {pre_norip_toplot, cond_norip_toplot, post_norip_toplot};

%% Figure
f1 = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
ax = arrayfun(@(i) subplot(2,3,i, 'NextPlot', 'add', 'Box', 'off'), 1:6);
for isess = 1:length(titles)
    axes(ax(isess));
    DatNormZ = onRipples{isess};
    imagesc(-wi*1e3:binsize*1e3:wi*1e3,1:size(DatNormZ,1),DatNormZ);
    caxis([-0.3 0.3])
    axis ij
    xlim([-wi*1e3 wi*1e3])
    ylim([0 length(DatNormZ)])
    ylabel('Neuron #')
    xlabel('Time around a ripple (ms)')
    title(titles{isess})
    makepretty
    
    axes(ax(isess+length(titles)))
    DatNormZ = noRipples{isess};
    imagesc(-wi*1e3:binsize*1e3:wi*1e3,1:size(DatNormZ,1),DatNormZ);
    caxis([-0.3 0.3])
    axis ij
    xlim([-wi*1e3 wi*1e3])
    ylim([0 length(DatNormZ)])
    ylabel('Neuron #')
    xlabel('Time around a no-ripple timepoint (ms)')
    title(titles{isess})
    makepretty
end

%% Save figure
if IsSave
    foldertosave = ChooseFolderForFigures_DB('Spikes');
    saveas(f1, [foldertosave '/Ripples/PETHOnRipples3.fig']);
    saveFigure(f1, 'PETHOnRipples3', [foldertosave '/Ripples']);
end




