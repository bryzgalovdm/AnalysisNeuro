

function makepretty_DB()
% set some graphical attributes of the current axis

set(get(gca, 'XLabel'), 'FontSize', 22);
set(get(gca, 'YLabel'), 'FontSize', 22);
set(gca, 'FontSize', 20);
box off
set(gca,'Linewidth',2)
set(get(gca, 'Title'), 'FontSize', 20);

ch = get(gca, 'Children');

for c = 1:length(ch)
    thisChild = ch(c);
    if strcmp('line', get(thisChild, 'Type'))
        if strcmp('.', get(thisChild, 'Marker'))
            % if get(thisChild, 'MarkerSize')<15
            %             set(thisChild, 'MarkerSize', 15);
            % end
        end
        if strcmp('-', get(thisChild, 'LineStyle'))
            set(thisChild, 'LineWidth', 2);
        end
    elseif strcmp('stair', get(thisChild, 'Type'))
                    set(thisChild, 'LineWidth', 2);
    elseif strcmp('errorbar', get(thisChild, 'Type'))
                    set(thisChild, 'LineWidth', 2);

    end
end
