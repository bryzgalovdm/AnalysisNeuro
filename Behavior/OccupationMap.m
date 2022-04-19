function OccupMap = OccupationMap(mice, expe, epochname, varargin)
% 
% This function returns a sum of occupation maps for a certain type of
% recording sesion in a certain type of an experiment
% 
% INPUT
%   
%   mice                    mice to analyse
%   expe                    type of experiment ('PAG', or 'MFB', or 'Novel')
%   epochname               type of a session to get ('Hab', or 'PreTest', or 'UMaze',
%                           or 'Cond', or 'Task', or 'PostTest')
%   IsSave (optional)       Whether to save a figure or not
%                           (default = false)
% 
% OUTPUT
% 
%   OccupMap                resulting occupation map
% 
% EXAMPLE
% 
%   OccupMap = OccupationMap(mice_PAG, 'PAG', 'Hab')
%   OccupMap = OccupationMap(mice_PAG, 'PAG', 'Hab', 'IsSave', true)
% 
% SEE ALSO
%
%   ReturnMnemozyneEpochs, OccupationMap3Stages
% 
% By Dima Bryzgalov, MOBS team, Paris,
% 08-09/2021
% github.com/bryzgalovdm


%% Parameters
sav=0;
sizeMap = 50;
smo = 1.5;
% Maze borderss
mazeMap = [4 6; 4 59; 59 59; 59 6; 38 6; 38 42; 23 42; 23 6; 4 6];
ShockZoneMap = [4 6; 4 30; 23 30; 23 6; 4 6];

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

%% Manage experiment
FigName = 'Occup';
if strcmp(expe, 'PAG')
    fetchpaths = 'UMazePAG';
elseif strcmp(expe, 'MFB')
    fetchpaths = 'StimMFBWake';
elseif strcmp(expe, 'Novel')
    fetchpaths = 'Novel';    
end

%% Manage epoch names
hab = false;
pre = false;
umaze = false;
cond = false;
task = false;
post = false;

if strcmp(epochname, 'Hab')
    hab = true;
elseif strcmp(epochname, 'PreTest')
    pre = true;
elseif strcmp(epochname, 'UMaze')
    umaze = true;
elseif strcmp(epochname, 'Cond')
    cond = true;
elseif strcmp(epochname, 'Task')
    task = true;
elseif strcmp(epochname, 'PostTest')
    post = true;
end


%% Allocate data arrays
% Get paths of each individual mouse
Dir = PathForExperimentsERC(fetchpaths);
Dir = RestrictPathForExperiment(Dir,'nMice',mice);
numsessions = CountNumSesionsERC(Dir);

% Data
b = cell(length(numsessions), 1);

% Auxiliary
Epoch = cell(length(numsessions), 1);
OccupMap_temp = cell(length(numsessions), 1);

%% load data
cnt=1;
for imouse=1:length(Dir.path)
    for isession = 1:length(Dir.path{imouse})
        b{cnt} = load([Dir.path{imouse}{isession} 'behavResources.mat'],...
            'AlignedXtsd', 'AlignedYtsd', 'SessionEpoch');
        cnt=cnt+1;
    end
end

%% Get epoch
for isession = 1:numsessions
    if strcmp(expe, 'Novel')
        [HabEpoch, PreEpoch, UMazeEpoch, CondEpoch, TaskEpoch] = ReturnMnemozyneEpochs(b{isession}.SessionEpoch);
    else
        [HabEpoch, PreEpoch, UMazeEpoch, CondEpoch, TaskEpoch, PostEpoch] = ReturnMnemozyneEpochs(b{isession}.SessionEpoch);
    end
    
    if hab
        Epoch{isession} = HabEpoch;
    elseif pre
        Epoch{isession} = PreEpoch;
    elseif umaze
        Epoch{isession} = UMazeEpoch;
    elseif cond
        Epoch{isession} = CondEpoch;
    elseif task
        Epoch{isession} = TaskEpoch;
    elseif post
        Epoch{isession} = PostEpoch;
    end
end

%% Calculate occupamcy
for isession = 1:numsessions
    
    OccupMap_temp{isession} = hist2d(Data(Restrict(b{isession}.AlignedXtsd,Epoch{isession})),...
            Data(Restrict(b{isession}.AlignedYtsd,Epoch{isession})),sizeMap,sizeMap);
        OccupMap_temp{isession} = OccupMap_temp{isession}/sum(OccupMap_temp{isession}(:));
        largerMatrix = zeros(sizeMap+floor(sizeMap/4),sizeMap+floor(sizeMap/4));
        largerMatrix(1+floor(sizeMap/8):sizeMap+floor(sizeMap/8),1+floor(sizeMap/8):sizeMap+floor(sizeMap/8)) = (OccupMap_temp{isession})';
        OccupMap_temp{isession}=SmoothDec(largerMatrix,[smo,smo]);
end

OccupMap = OccupMap_temp{1};
if length(OccupMap_temp) > 1
    for isession = 2:length(OccupMap_temp)
        OccupMap = OccupMap + OccupMap_temp{isession};
    end
end


%% Figure
f = figure('units', 'normalized', 'outerposition', [0 0 0.36 0.65]);
imagesc(OccupMap);
axis xy
hold on
plot(mazeMap(:,1),mazeMap(:,2),'w','LineWidth',3)
if strcmp(expe, 'PAG')
    plot(ShockZoneMap(:,1),ShockZoneMap(:,2),'r','LineWidth',3);
    set(gca,'XTick', [], 'XTickLabel',{},'YTick', [], 'YTickLabel',{});
    caxis([0 .08]);
elseif strcmp(expe, 'MFB')
    plot(ShockZoneMap(:,1),ShockZoneMap(:,2),'g','LineWidth',3);
    set(gca,'XTick', [], 'XTickLabel',{},'YTick', [], 'YTickLabel',{}, 'XDir', 'reverse');
    caxis([0 .03]);
elseif strcmp(expe, 'Novel')
    set(gca,'XTick', [], 'XTickLabel',{},'YTick', [], 'YTickLabel',{});
    caxis([0 .03]);    
end
colormap hot
makepretty_DB

%% Save figure
if sav
    foldertosave = ChooseFolderForFigures_DB('Behavior');
    saveas(f, [foldertosave filesep FigName '_' expe '_' epochname]);
    saveFigure(f, [FigName '_' expe '_' epochname], foldertosave);
end

end