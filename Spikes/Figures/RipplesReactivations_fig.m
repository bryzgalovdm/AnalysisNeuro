%% Parameters
mazeMap = [8 8; 8 56; 55 56; 55 8; 39 8; 39 42; 24 42; 24 8; 8 8];
ShockZoneMap = [8 8; 8 30; 24 30; 24 8; 8 8];
SafeZoneMap = [55 8; 55 30; 39 30; 39 8; 55 8];

l = [4 62];
%%
mice = [797 798 828 861 882 905 906 911 912 977 994];
[map_react, allfields, coefs] = MapPCReactivatedRipples(mice);

%%
dat = {map_react.Pre, map_react.Cond, map_react.Post};
norm = {coefs.Pre*allfields, coefs.Post*allfields, coefs.Post*allfields};
tits = {'PreSleep', 'Conditioning', 'PostSleep'};
f1 = figure('units', 'normalized', 'outerposition', [0 0 0.9 0.4]);
ax = arrayfun(@(i) subplot(1,3,i, 'NextPlot', 'add', 'Box', 'off'), 1:3);
for ifig = 1:length(ax)
    axes(ax(ifig));
    
    dat{ifig}([1:8 56:62],:) = 0;
    dat{ifig}(:,[1:8 55:62]) = 0;
    dat{ifig}(1:42,24:39) = 0;

    imagesc(dat{ifig}./norm{ifig});
    %     imagesc(dat{ifig});
    plot(mazeMap(:,1),mazeMap(:,2),'w','LineWidth',4)
    plot(ShockZoneMap(:,1),ShockZoneMap(:,2),'r','LineWidth',4)
    title(tits{ifig}, 'FontSize', 16);
    set(gca, 'XTickLabel', {}, 'YTickLabel', {});
    xlim(l);
    ylim(l);
end
% mtit('Probability to participate in ripples for place cells', 'FontSize', 18, 'xoff', 0, 'yoff', 0.05)
