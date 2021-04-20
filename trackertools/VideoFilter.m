%
% Filters video by taking difference between narrow and broad 2-D filters.
%
% USAGE:
%   output = VideoFilter(input,narrow,broad)
%
%   input: input video frame
%   narrow: pixels for narrow filter (default = 3)
%   broad: pixels for broad filter (default = 20)
%   output: filtered video frame

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 28-May-2009 12:37:41 
%---------------------------- 

function output = VideoFilter(input,narrow,broad)

if ~exist('narrow') || isempty(narrow) narrow = 3; end
if ~exist('broad') || isempty(broad) broad = 20; end

F0 = input;
F1a = imfilter(F0,fspecial('average',narrow));
F1b = imfilter(F0,fspecial('average',broad));
F2 = F1a - F1b;

output = F2;
