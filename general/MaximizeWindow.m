%
% Maximize figure window size
%
% USAGE:
%   MaximizeWindow(figh)
% 
%   figh: figure handle

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 07-Mar-2010 10:54:16 
%---------------------------- 

function MaximizeWindow(figh)

screensize = get(0,'ScreenSize');

winmenu = 30; % pixels for windows bar at bottom
menubar = 74; % pixels for standard menu at top
figurepos = screensize + [0 winmenu 0 -winmenu-menubar];

if nargin < 1
    set(gcf,'Position',figurepos);
else
    for i = 1:numel(figh)
        figure(figh(i));
        set(gcf,'Position',figurepos);
    end
end
       