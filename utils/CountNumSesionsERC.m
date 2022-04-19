function [numsessions, micename] = CountNumSesionsERC(Dir, varargin)



%% Parameters
verbose = false;

%% Optional parameters
for i=1:2:length(varargin)
    switch(lower(varargin{i}))
        case 'verbose'
            verbose = varargin{i+1};
            if IsII ~= 1 && IsII ~= 0
                error('Incorrect value for property ''Verbose'' (type ''help CountNumSesionsERC'' for details).');
            end
    end
end

%% Count Num Sessions
numsessions = 0;
for imouse = 1:length(Dir.path)
    for isession = 1:length(Dir.path{imouse})
        numsessions = numsessions + 1;
        micename{isession} = Dir.name{imouse};
    end
end

if verbose
    print(['Analysing ' num2str(numsessions) ' sessions....']);
end
end







