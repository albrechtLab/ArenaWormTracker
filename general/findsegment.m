%
% Analyze vector for consecutive values.
%
% USAGE:
%   [pos,seglen,data2] = findsegment(data);
%   
%   data: vector of index values
%   pos: ending value for each group of consecutive values
%   seglen: length of each group of consecutive values
%
%   EXAMPLE:  [a,b] = findsegment(find([0 0 0 1 1 1 0 0 0 1]))
%
%                 a =
%                      6    10
% 
%                 b =
%                      3     1
%
%           i.e. first segment ends at position 6, 3 elements long.
     
%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 07-Feb-2008 11:31:55 
%---------------------------- 

function [pos,seglen,data2] = findsegment(data);

segnum = 1; count = 1; 
data2 = []; pos = []; seglen = [];

if length(data) > 0
    data(length(data)+1) = data(length(data)) + 2;

    prevframe = data(1);
    for i = 2:length(data)
        frame = data(i);
        if frame == prevframe+1;
            count = count+1;
        else
            pos(segnum) = prevframe;
            seglen(segnum) = count;
            data2 = [data2, count*ones(1,count)];
            segnum = segnum + 1;
            count = 1;
        end
        prevframe = frame;
    end
end