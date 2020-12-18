%% Parameters
IsSave = true;

ShockArea = [7 8; 25 30];
SafeArea = [38 8; 55 30];

%% Get data
mice = [797 798 828 861 882 905 906 911 912 977 994];
ripprob_pc = CalcParticipationProbab_ripples(mice, 'SepArea', {ShockArea, SafeArea});

%% Parse and flatten data arrays
% Allocate
% Shock
ripprob_shz.Pre = cell(length(ripprob_pc.Pre),1);
ripprob_shz.Cond = cell(length(ripprob_pc.Cond),1);
ripprob_shz.Post = cell(length(ripprob_pc.Post),1);
% Safe
ripprob_saz.Pre = cell(length(ripprob_pc.Pre),1);
ripprob_saz.Cond = cell(length(ripprob_pc.Cond),1);
ripprob_saz.Post = cell(length(ripprob_pc.Post),1);
% Others
ripprob_oc.Pre = cell(length(ripprob_pc.Pre),1);
ripprob_oc.Cond = cell(length(ripprob_pc.Cond),1);
ripprob_oc.Post = cell(length(ripprob_pc.Post),1);
% Parse
for imouse = 1:length(mice)
    % Shock
    ripprob_shz.Pre{imouse} = ripprob_pc.Pre{imouse}{1};
    ripprob_shz.Cond{imouse} = ripprob_pc.Cond{imouse}{1};
    ripprob_shz.Post{imouse} = ripprob_pc.Post{imouse}{1};
    % Safe
    ripprob_saz.Pre{imouse} = ripprob_pc.Pre{imouse}{2};
    ripprob_saz.Cond{imouse} = ripprob_pc.Cond{imouse}{2};
    ripprob_saz.Post{imouse} = ripprob_pc.Post{imouse}{2};
    % Others
    ripprob_oc.Pre{imouse} = ripprob_pc.Pre{imouse}{3};
    ripprob_oc.Cond{imouse} = ripprob_pc.Cond{imouse}{3};
    ripprob_oc.Post{imouse} = ripprob_pc.Post{imouse}{3};
end
% Flatten
% Shock and others
sz_pre = flatten_cellarray(ripprob_shz.Pre, 1);
sz_cond = flatten_cellarray(ripprob_shz.Cond, 1);
sz_post = flatten_cellarray(ripprob_shz.Post, 1);
% Safe and others
sa_pre = flatten_cellarray(ripprob_saz.Pre, 1);
sa_cond = flatten_cellarray(ripprob_saz.Cond, 1);
sa_post = flatten_cellarray(ripprob_saz.Post, 1);
%Others
oc_pre = flatten_cellarray(ripprob_oc.Pre, 1);
oc_cond = flatten_cellarray(ripprob_oc.Cond, 1);
oc_post = flatten_cellarray(ripprob_oc.Post, 1);

%% Plot participation probability during ripples
disp(['Number of shock zone place cells: ' num2str(size(sz_pre,1))]);
disp(['Number of safe zone cells: ' num2str(size(sa_pre,1))]);
disp(['Number of other neurons: ' num2str(size(oc_pre,1))]);

% Pre
% Number of stimulations in REM - separate figure for only SZ cells
Pl = {{sz_pre, sa_pre, oc_pre}, {sz_cond, sa_cond, oc_cond}, {sz_post, sa_post, oc_post}};
Cols = {[0.9 0.2 0.2], [0.2 0.2 0.9], [0.6 0.6 0.6]};
tits = {'PreSleep', 'Conditioning', 'PostSleep'};
f1 = figure('units', 'normalized', 'outerposition', [0 0 0.9 0.75]);
ax = arrayfun(@(i) subplot(1,3,i, 'NextPlot', 'add', 'Box', 'off'), 1:3);
yl =[0 0.3];
for ifig = 1:length(ax)
    axes(ax(ifig));
    MakeBoxPlot_DB(Pl{ifig},Cols,1:3, {'Shock', 'Safe', 'Other'}, 0);
    set(gca, 'FontSize', 16, 'FontWeight', 'bold');
    if ifig == 1
        ylabel('Participation probability');
    end
    title(tits{ifig}, 'FontSize', 16);
    ylim(yl);
end
% mtit('Probability to participate in ripples for place cells', 'FontSize', 18, 'xoff', 0, 'yoff', 0.05)

%% Save figure
if IsSave
    foldertosave = ChooseFolderForFigures_DB('Spikes');
    saveas(f1, [foldertosave '/ParticipProbRipples_boxplot.fig']);
    saveFigure(f1, 'ParticipProbRipples_boxplot', foldertosave);
end
