function SleepBeforeAfter(mice, expe, varargin)

%
% This functions plots barplots comparing sleep architecture before and after
% experient
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
SP = true;
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
        case 'showpoints'
            SP = varargin{i+1};
            if SP ~= 1 && SP ~= 0
                error('Incorrect value for property ''ShowPoints'' (type ''help PreTestCharacteristics'' for details).');
            end
    end
end

%% Manage experiment
FigName = 'SleepArchitecture';
if strcmp(expe, 'PAG')
    fetchpaths = 'UMazePAG';
elseif strcmp(expe, 'MFB')
    fetchpaths = 'StimMFBWake';
elseif strcmp(expe, 'Novel')
    fetchpaths = 'Novel';    
end

%% Allocate data arrays
% Get paths of each individual mouse
Dir = PathForExperimentsERC(fetchpaths);
Dir = RestrictPathForExperiment(Dir,'nMice',mice);
numsessions = CountNumSesionsERC(Dir);

% Data
b = cell(length(numsessions), 1);
sleepscored = cell(length(numsessions), 1);
%% Load the data

cnt=1;
for imouse=1:length(Dir.path)
    for isession = 1:length(Dir.path{imouse})
        b{cnt} = load([Dir.path{imouse}{isession} 'behavResources.mat'],'SessionEpoch');
        try
            sleepscored{cnt} = load([Dir.path{imouse}{isession} 'SleepScoring_OBGamma.mat'],'SWSEpoch','REMEpoch','Wake');
        catch
            sleepscored{cnt} = load([Dir.path{imouse}{isession} 'SleepScoring_Accelero.mat'],'SWSEpoch','REMEpoch','Wake'); % REM and Sleep are not used
        end
        cnt=cnt+1;
    end
end


%% Prepare arrays
NREMPre = cell(1, length(b));
REMPre = cell(1, length(b));
WakePre = cell(1, length(b));
NREMPost = cell(1, length(b));
REMPost = cell(1, length(b));
WakePost = cell(1, length(b));

PreSleep = nan(1, length(b));
PostSleep = nan(1, length(b));

PreWake = nan(1,length(b));
PreNREM = nan(1,length(b));
PreREM = nan(1,length(b));

PostWake = nan(1,length(b));
PostNREM = nan(1,length(b));
PostREM = nan(1,length(b));

%% Create epochs

for i=1:numsessions
    try % Restrict to the first two hours if applicable
        NREMPre{i} = and(RestrictToTime(b{i}.SessionEpoch.PreSleep, SleepTimeToRestrict), sleepscored{i}.SWSEpoch);
        REMPre{i} = and(RestrictToTime(b{i}.SessionEpoch.PreSleep, SleepTimeToRestrict), sleepscored{i}.REMEpoch);
        WakePre{i} = and(RestrictToTime(b{i}.SessionEpoch.PreSleep, SleepTimeToRestrict), sleepscored{i}.Wake);
        PreSleep(i) = SleepTimeToRestrict/1e4;
    catch
        NREMPre{i} = and(b{i}.SessionEpoch.PreSleep, sleepscored{i}.SWSEpoch);
        REMPre{i} = and(b{i}.SessionEpoch.PreSleep, sleepscored{i}.REMEpoch);
        WakePre{i} = and(b{i}.SessionEpoch.PreSleep, sleepscored{i}.Wake);
        PreSleep(i) = End(b{i}.SessionEpoch.PreSleep, 's') - Start(b{i}.SessionEpoch.PreSleep, 's');
    end
    try % Restrict to the first two hours if applicable
        NREMPost{i} = and(RestrictToTime(b{i}.SessionEpoch.PostSleep, SleepTimeToRestrict), sleepscored{i}.SWSEpoch);
        REMPost{i} = and(RestrictToTime(b{i}.SessionEpoch.PostSleep, SleepTimeToRestrict), sleepscored{i}.REMEpoch);
        WakePost{i} = and(RestrictToTime(b{i}.SessionEpoch.PostSleep, SleepTimeToRestrict), sleepscored{i}.Wake);
        PostSleep(i) = SleepTimeToRestrict/1e4;
    catch
        NREMPost{i} = and(b{i}.SessionEpoch.PostSleep, sleepscored{i}.SWSEpoch);
        REMPost{i} = and(b{i}.SessionEpoch.PostSleep, sleepscored{i}.REMEpoch);
        WakePost{i} = and(b{i}.SessionEpoch.PostSleep, sleepscored{i}.Wake);
        PostSleep(i) = End(b{i}.SessionEpoch.PostSleep, 's') - Start(b{i}.SessionEpoch.PostSleep, 's');
    end
end


%% Calculate

% Individually
for i=1:numsessions
    PreWake(i) = (sum(End(WakePre{i}, 's') - Start(WakePre{i}, 's'))) / PreSleep(i) * 100;
    PreNREM(i) = sum(End(NREMPre{i}, 's') - Start(NREMPre{i}, 's')) / PreSleep(i) * 100;
    PreREM(i) = sum(End(REMPre{i}, 's') - Start(REMPre{i}, 's')) / PreSleep(i) * 100;
    
    PostWake(i) = sum(End(WakePost{i}, 's') - Start(WakePost{i}, 's')) / PostSleep(i) * 100;
    PostNREM(i) = sum(End(NREMPost{i}, 's') - Start(NREMPost{i}, 's')) / PostSleep(i) * 100;
    PostREM(i) = sum(End(REMPost{i}, 's') - Start(REMPost{i}, 's')) / PostSleep(i) * 100;
end

%% F
f1 = figure('units','normalized', 'outerposition', [0 1 0.41 0.82]);
subplot(222)
[~,h] = PlotErrorBarN_DB([PreWake' PostWake'],...
    'barcolors', [0 0 0], 'barwidth', 0.5, 'newfig', 0, 'showpoints', SP);
h.FaceColor = 'flat';
h.CData(1,:) = [.9856, .7372, .2537];
if strcmp(expe, 'PAG')
%     h.CData(2,:) = [1 0 1];
    h.CData(2,:) = [.9 0 0];
elseif strcmp(expe, 'MFB')
%     h.CData(2,:) = [.16 .95 1];
    h.CData(2,:) = [0 .9 0];
elseif strcmp(expe, 'Novel')
%     h.CData(2,:) = [0 0 0];
    h.CData(2,:) = [1 1 1];
end
set(gca,'Xtick',[1:2],'XtickLabel',{});
ylabel('%')
ylim([0 80])
title('Wake');
makepretty_DB

subplot(223)
[~,h] = PlotErrorBarN_DB([PreNREM' PostNREM'], 'barcolors',...
    [0 0 0], 'barwidth', 0.5, 'newfig', 0, 'showpoints',SP);
h.FaceColor = 'flat';
h.CData(1,:) = [.9856, .7372, .2537];
if strcmp(expe, 'PAG')
%     h.CData(2,:) = [1 0 1];
    h.CData(2,:) = [.9 0 0];
elseif strcmp(expe, 'MFB')
%     h.CData(2,:) = [.16 .95 1];
    h.CData(2,:) = [0 .9 0];
elseif strcmp(expe, 'Novel')
%     h.CData(2,:) = [0 0 0];
    h.CData(2,:) = [1 1 1];
end
set(gca,'Xtick',[1:2],'XtickLabel',{});
ylabel('%')
ylim([0 80])
title('NREM');
makepretty_DB

subplot(224)
[~,h] = PlotErrorBarN_DB([PreREM' PostREM'],...
    'barcolors', [0 0 0], 'barwidth', 0.5, 'newfig', 0, 'showpoints',SP);
h.FaceColor = 'flat';
h.CData(1,:) = [.9856, .7372, .2537];
if strcmp(expe, 'PAG')
%     h.CData(2,:) = [1 0 1];
    h.CData(2,:) = [.9 0 0];
elseif strcmp(expe, 'MFB')
%     h.CData(2,:) = [.16 .95 1];
    h.CData(2,:) = [0 .9 0];
elseif strcmp(expe, 'Novel')
%     h.CData(2,:) = [0 0 0];
    h.CData(2,:) = [1 1 1];
end
ylim([0 12])
set(gca,'Xtick',[1:2],'XtickLabel',{});
ylabel('%')
title('REM');
makepretty_DB

if sav
    foldertosave = ChooseFolderForFigures_DB('Sleep');
    saveas(f1, [foldertosave filesep FigName '_' expe]);
    saveFigure(f1, [FigName '_' expe], foldertosave);
end

end