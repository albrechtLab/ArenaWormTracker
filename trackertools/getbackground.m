%
% Calculate average background image for specified video range and interval.
%
% USAGE:
%   background = getbackground(RawMovieName,StartFrame,EndFrame,FrameInterval)
% 
%   RawMovieName: movie filename
%   filename: name to display on figure

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 26-Mar-2008 20:44:44 
%---------------------------- 
function background = getbackground(RawMovieName,StartFrame,EndFrame,FrameInterval)

MovieObj = VideoReader(RawMovieName);
m = MovieObj.Width;
n = MovieObj.Height;
cdatasum = zeros(n,m,'double');   %for 8-bit movies

FrameNum = MovieObj.NumberOfFrames;
if nargin < 4, FrameInterval = 20; end
if nargin < 3, EndFrame = FrameNum; end
if nargin < 2, StartFrame = 1; end

disp(['Background calculating from ',int2str(StartFrame),' to ',int2str(EndFrame),' in increments of ',int2str(FrameInterval)]);

progbars = 10;
for Frame = StartFrame:FrameInterval:EndFrame
    Mov = read(MovieObj, Frame);
    MovX64 = double(Mov(:,:,2))/255;
    cdatasum = cdatasum + MovX64;
    if mod(Frame * progbars,(EndFrame-StartFrame+1)) < progbars*FrameInterval fprintf(':'); end
end
fprintf('\n');

cdataaverage = cdatasum./round((EndFrame-StartFrame+1)/FrameInterval);
background = uint8(round(cdataaverage*255));
