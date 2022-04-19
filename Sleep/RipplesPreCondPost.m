function RipplesPreCondPost(mice, expe, varargin)
%
% This functions plots barplots comparing ripples density during 
% pre-tests, learning and post-tests
%
%  INPUT
%       
%       mice            array with mice numbers that will harvested from ERC
%                       PathForExperiments
%       expe            type of ERC experiment ('PAG', 'MFB' or 'Novel')
%       IsSave          (optional) true or 1 if to save the figure (default=false)
% 
% 
%  OUTPUT
%
%       Figure
% 
% Coded by Dima Bryzgalov and Samuel Laventure, MOBS team, Paris, France
% 09/2021
% github.com/bryzgalovdm

%% Parameters
sav=0;
% Sleep time to restrict
SleepTimeToRestrict = 2*60*60*1e4; % 2 hours

%% Optional Arguments
for i=1:2:length(varargin)    
    switch(lower(varargin{i})) 
        case 'issave'
            sav = varargin{i+1};
            if sav ~= 1 && sav ~= 0
                error('Incorrect value for property ''IsSave'' (type ''help PreTestCharacteristics'' for details).');
            end
    end
end

%% Optional Arguments
for i=1:2:length(varargin)    
    switch(lower(varargin{i})) 
        case 'issave'
            sav = varargin{i+1};
            if sav ~= 1 && sav ~= 0
                error('Incorrect value for property ''IsSave'' (type ''help PreTestCharacteristics'' for details).');
            end
        case 'showpoints'
            SP = varargin{i+1};
            if SP ~= 1 && SP ~= 0
                error('Incorrect value for property ''ShowPoints'' (type ''help PreTestCharacteristics'' for details).');
            end
    end
end

%% Manage experiment
if strcmp(expe, 'PAG')
    fetchpaths = 'UMazePAG';
    FigName = 'RipplesPAG';
elseif strcmp(expe, 'MFB')
    fetchpaths = 'StimMFBWake';
    FigName = 'RipplesMFB';
elseif strcmp(expe, 'Novel')
    fetchpaths = 'Novel';
    FigName = 'RipplesNovel';
end

%% Allocate data arrays
% Get paths of each individual mouse
Dir = PathForExperimentsERC(fetchpaths);
Dir = RestrictPathForExperiment(Dir,'nMice',mice);
numsessions = CountNumSesionsERC(Dir);

% Data
Session = cell(length(numsessions), 1);
% Data
b = cell(length(numsessions), 1);
Sleep = cell(length(numsessions), 1);
Rip = cell(length(numsessions), 1);
%% Load the data

cnt=1;
for imouse=1:length(Dir.path)
    for isession = 1:length(Dir.path{imouse})
        try
            Rip{cnt} = load([Dir.path{imouse}{isession} 'SWR.mat'], 'ripples');
        catch
            Rip{cnt} = load([Dir.path{imouse}{isession} 'Ripples.mat'], 'ripples');
        end
        Session{cnt} = load([Dir.path{imouse}{isession} 'behavResources.mat'], 'behavResources', 'SessionEpoch', 'FreezeAccEpoch');
        try
            Sleep{cnt} = load([Dir.path{imouse}{isession} 'SleepScoring_OBGamma.mat'],'SWSEpoch','REMEpoch','Wake');
        catch
            Sleep{cnt} = load([Dir.path{imouse}{isession} 'SleepScoring_Accelero.mat'],'SWSEpoch','REMEpoch','Wake'); % REM and Sleep are not used
        end
        cnt=cnt+1;
    end
end

%% Preallocate arrays
id_Pre = cell(1,length(Session));
id_Cond = cell(1,length(Session));
id_Post = cell(1,length(Session));

CondFreezeEpoch = cell(1,length(Session));
PreNREMEpoch = cell(1,length(Session));
PostNREMEpoch = cell(1,length(Session));

PreRipples = cell(1,length(Session));
PostRipples = cell(1,length(Session));
CondRipples = cell(1,length(Session));
EarlyNREM_pre = cell(1,length(Session));
EarlyNREM_post = cell(1,length(Session));

PreR_density = nan(1,length(Session));
PostR_density = nan(1,length(Session));
CondR_density = nan(1,length(Session));
EarlyPreR_density = nan(1,length(Session));
EarlyPostR_density = nan(1,length(Session));

%% Find indices of PreTests and PostTest session in the structure

for i=1:length(Session)
    id_Pre{i} = FindSessionID_ERC(Session{i}.behavResources, 'TestPre');
    id_Cond{i} = FindSessionID_ERC(Session{i}.behavResources, 'Cond');
    id_Post{i} = FindSessionID_ERC(Session{i}.behavResources, 'TestPost');
end

%% Prepare intervalSets for ripples

for i=1:numsessions
    if strcmp(expe, 'PAG')
        [~, ~, CondEpoch{i}, ~, ~] = ReturnMnemozyneEpochs(Session{i}.SessionEpoch);
        CondFreezeEpoch{i} = and(CondEpoch{i}, Session{i}.FreezeAccEpoch);
    else
        [~, ~, CondEpoch{i}, ~, ~] = ReturnMnemozyneEpochs(Session{i}.SessionEpoch);
    end
    try
        PreNREMEpoch{i} = and(RestrictToTime(b{i}.SessionEpoch.PreSleep, SleepTimeToRestrict), Sleep{i}.SWSEpoch);
        PostNREMEpoch{i} = and(RestrictToTime(b{i}.SessionEpoch.PostSleep, SleepTimeToRestrict), Sleep{i}.SWSEpoch);
    catch
        PreNREMEpoch{i} = and(Session{i}.SessionEpoch.PreSleep, Sleep{i}.SWSEpoch);
        PostNREMEpoch{i} = and(Session{i}.SessionEpoch.PostSleep, Sleep{i}.SWSEpoch);
    end
end


%% Calculate ripples density

% Restrict sleeps to first 30 min
for i = 1:numsessions
    try
        EarlyNREM_pre{i} = RestrictToTime(PreNREMEpoch{i}, 30*60*1e4);
        EarlyNREM_post{i} = RestrictToTime(PostNREMEpoch{i}, 30*60*1e4);
    catch
        EarlyNREM_pre{i} = PreNREMEpoch{i};
        EarlyNREM_post{i} = PostNREMEpoch{i};
    end
end

% Extract ripples for presleep and postsleep (during detected sleep only)
for i = 1:length(Session)
    ripplesPeak=ts(Rip{i}.ripples(:,2)*1e4);
    PreRipples{i}=Restrict(ripplesPeak,PreNREMEpoch{i});
    PostRipples{i}=Restrict(ripplesPeak,PostNREMEpoch{i});
    if strcmp(expe, 'PAG')
        CondRipples{i}=Restrict(ripplesPeak,CondFreezeEpoch{i});
    else
        CondRipples{i}=Restrict(ripplesPeak,CondEpoch{i});
    end
    
    
    PreR_density(i) = length(Range(PreRipples{i})) / tot_length(PreNREMEpoch{i}, 's');
    PostR_density(i) = length(Range(PostRipples{i})) / tot_length(PreNREMEpoch{i}, 's');
    if strcmp(expe, 'PAG')
        CondR_density(i) = length(Range(CondRipples{i})) / tot_length(CondFreezeEpoch{i}, 's');
    else
        CondR_density(i) = length(Range(CondRipples{i})) / tot_length(CondEpoch{i}, 's');
    end
    
    EarlyPreR_density(i) = length(Restrict(ripplesPeak,EarlyNREM_pre{i})) / tot_length(EarlyNREM_pre{i}, 's');
    EarlyPostR_density(i) = length(Restrict(ripplesPeak,EarlyNREM_post{i})) / tot_length(EarlyNREM_post{i}, 's');
    
end


%% Plot
if strcmp(expe, 'PAG')
%     cols = {[.9856, .7372, .2537], [1 0 1]};
    cols = {[.9856, .7372, .2537], [0.9 0 0]};
elseif strcmp(expe, 'MFB')
%     cols = {[.9856, .7372, .2537], [.16 .95 1]};
    cols = {[.9856, .7372, .2537], [0 0.9 0]};
elseif strcmp(expe, 'Novel')
%     cols = {[.9856, .7372, .2537], [0 0 0]};
    cols = {[.9856, .7372, .2537], [1 1 1]};
end
data_toplot = {PreR_density, PostR_density};

f1 = figure('units', 'normalized', 'outerposition', [0 0 0.35 0.6]);
[b, p]=MakeBoxPlot_DB(data_toplot, cols, [1 2], {'PreSleep', 'PostSleep'}, 1);
for iplot = 1:length(b)
    b{iplot}.handles.upperWhiskers.Visible = 'off';
    b{iplot}.handles.lowerWhiskers.Visible = 'off';
    b{iplot}.handles.lowerWhiskerTips.Visible = 'off';
    b{iplot}.handles.upperWhiskerTips.Visible = 'off';
    if strcmp(expe, 'Novel') && iplot == 2
        b{iplot}.lineColor = 'k';
    else
        b{iplot}.lineColor = 'w';
    end
    b{iplot}.boxAlpha = .7;
end
for iplot = 1:length(p)
    set(p{iplot}{1}, 'MarkerSize', 40);
end
p = DoWilcoxonOnArray(data_toplot, {[1 2]});
if p <= 0.05
    sigstar_DB([1 2],p,0,'LineWigth',16,'StarSize',24);
end
ylabel('Ripples/s');
ylim([.37 1.6])
title('Full NREM ripples')
makepretty_DB
if sav
    foldertosave = ChooseFolderForFigures_DB('LFP');
    saveas(f1, [foldertosave filesep 'Rip_density_full_' expe '.fig']);
    saveFigure(f1, ['Rip_density_full_' expe], foldertosave);
end

% First 30 min
data_toplot = {EarlyPreR_density, EarlyPostR_density};

f2 = figure('units', 'normalized', 'outerposition', [0 0 0.35 0.6])
[b, p]=MakeBoxPlot_DB(data_toplot, cols, [1 2], {'PreSleep', 'PostSleep'}, 1);
for iplot = 1:length(b)
    b{iplot}.handles.upperWhiskers.Visible = 'off';
    b{iplot}.handles.lowerWhiskers.Visible = 'off';
    b{iplot}.handles.lowerWhiskerTips.Visible = 'off';
    b{iplot}.handles.upperWhiskerTips.Visible = 'off';
    if strcmp(expe, 'Novel') && iplot == 2
        b{iplot}.lineColor = 'k';
    else
        b{iplot}.lineColor = 'w';
    end
    b{iplot}.boxAlpha = .7;
end
for iplot = 1:length(p)
    set(p{iplot}{1}, 'MarkerSize', 40);
end
p = DoWilcoxonOnArray(data_toplot, {[1 2]});
if p <= 0.05
    sigstar_DB([1 2],p,0,'LineWigth',16,'StarSize',24);
end
ylabel('Ripples/s');
ylim([.37 1.6])
% title('30 min of NREM ripples')
makepretty_DB
if sav
    foldertosave = ChooseFolderForFigures_DB('LFP');
    saveas(f2, [foldertosave filesep 'Rip_density_' expe '.fig']);
    saveFigure(f2, ['Rip_density_' expe], foldertosave);
end

end