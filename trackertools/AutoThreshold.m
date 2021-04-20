%
% Automatically determine optimal threshold for video frame processing.
%
% USAGE:
%   [threshold, animalpix] = AutoThreshold(Movdiv, thrlist)
% 
%   Movdiv: background-divided movie frame
%   thrlist: patthern of thresholds to test [default = 0.10:-0.005:0.005]
%   threshold: optimal threshold value
%   animalpix: size of average animal object in pix^2

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 02-Feb-2010 08:07:07 
%---------------------------- 

function [threshold, animalpix] = AutoThreshold(Movdiv, thrlist)
if nargin < 2 thrlist = 0.10:-0.005:0.030; end

pix = [];
N = [];
t = [];

F2 = VideoFilter(Movdiv);

for th = 1:length(thrlist)

    thr = thrlist(th);

    pix = [pix, sum(sum(F2 >= thr))];   
    [L,NUM] = bwlabel(im2bw(F2,thr)); 
    N = [N, NUM];
    if NUM > 500 break; end;
end

avgobjsize = pix ./ N;
threshold = thrlist(find(avgobjsize == max(avgobjsize)));
threshold = threshold(end);
[L,NUM] = bwlabel(im2bw(F2,threshold));

PROPS = regionprops(L);
objsize = [PROPS.Area];

hrange = 0:5:500;
shist = hist(objsize,hrange); shist(find(hrange <= 20)) = 0; smhist = smooth(shist,3);
animalpix = hrange(find(smhist == max(smhist))); animalpix = animalpix(1);
