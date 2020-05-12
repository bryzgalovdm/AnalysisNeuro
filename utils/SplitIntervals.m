function [I,Res] = SplitIntervals(is, len)

% [I,Res] = SplitIntervals(is, len) splits IntervalSet to interval sets of specific length
%
% This code makes your interval to be specific length (len) and give residual interval set (Res)
%
%
%
% USAGE
%
%    [I,Res] = SplitIntervals(is, len)
%
%    INPUTS
%
%    is         an interval set
%    len        desired total length of the resulting interval set in t
%
%    OUTPUTS
%
%    I      resulting splited interval sets
%    Res    residual interval set, less than len ([] if no resudue)
%
%
%    See also RestrictToTime, Restrict, intervalSet
%
% the current program is an addition to Francesco P. Battaglia's tsd plugin
% copyright (c) 2004 Francesco P. Battaglia, (c) 2019 Dmitri Bryzgalov
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% released under the GPL

% Find length of each interval
for i=1:length(Start(is))
    l=(End(is)-Start(is));
end

% Sanity check-up
if sum(l)<len
    error('Length of interval set is shorter than time you are trying to restrict it to');
end

% Housekeeping
count=0;
j=0;
numI = ceil(sum(l)/len);
I=cell(1,numI-1);

for i = 1:numI
    % Find the interval that exceeds the desired length
    if count == 0
        while count<len
            j=j+1;
            count = count + l(j);
        end
        last(i) = j;
    else
        if sum(l)-(i-1)*len > len
            count = l(last(i-1))-delta;
            j=last(i-1);
            while count < len
                j=j+1;
                count = count + l(j);
            end
            last(i) = j;
        else
            % Residual interval
            if i > 2
                if last(i-2) ~= last(i-1)
                    Res = or(intervalSet(Start(subset(is, last(i-1)))+delta, End(subset(is, last(i-1)))),...
                        subset(is, last(i-1)+1:length(l)));
                    break
                else
                    tempend = End(is);
                    tempstart = End(I{i-1});
                    Res = intervalSet(tempstart(end),tempend(end));
                    break
                end
            else
                Res = is-is_result;
                break
            end
        end
    end
    
    
    % Crop the interval that exceeds the desired length
    if count > len
        if i == 1
            almost = sum(l(1:last-1));
            delta = len - almost;
            
            is_almost = subset(is, 1:last-1);
            is_delta = intervalSet(Start(subset(is, last)), Start(subset(is, last))+delta);
            
            is_result = or(is_almost, is_delta);
        else
            if last(i-1) ~= last(i)
                is_beginning = intervalSet(Start(subset(is, last(i-1)))+delta, End(subset(is, last(i-1))));
                l_beg = End(is_beginning) - Start(is_beginning);
                almost = sum(l((last(i-1)+1):(last(i)-1)));
                delta = len - almost - l_beg;
                
                is_almost = subset(is, last(i-1)+1:last(i)-1);
                is_delta = intervalSet(Start(subset(is, last(i))), Start(subset(is, last(i)))+delta);
                
                is_result = or(is_beginning, or(is_almost, is_delta));
            else
                TTT = End(I{i-1});
                if length(TTT) == 1
                    is_result = intervalSet(TTT, TTT+len);
                else
                    is_result = intervalSet(TTT(end), TTT(end)+len);
                end
            end
        end
    elseif count == len
        if i == 1
            is_result = subset(is, 1:last);
        else
            if last(i-1) ~= last(i)
                is_beginning = intervalSet(Start(subset(is, last(i-1)))+delta, End(subset(is, last(i-1))));
                is_almost = subset(is, last(i-1):last(i));
                
                is_result = or(is_almost, is_beginning);
            else
                is_result = intervalSet(End(I{i-1}), End(I{i-1})+len);
            end
        end
    else
        error('The author is a fool and did not anticipate some important case. Let him know - dmitri.bryzgalov@espci.fr');
    end
    
    % Final assignment
    I{i} = is_result;
end

% Residual interval
if ~exist('Res')
    Res=[];
end

end


