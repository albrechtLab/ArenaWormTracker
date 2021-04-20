%
% Script to sequentially run all analysis scripts:
%     ArenaTracker
%     SegmentTracks
%     Ethogram
%     WormDensity
%  

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 12-Apr-2011 17:13:58 
%---------------------------- 

% Track Video File
File = ArenaTracker;

% Segment & Link Tracks
SegFileList = SegmentTracks({File.TrackFile});

% Obtain Ethogram
Ethogram(SegFileList);

% Summarize Data by binning over space and time
WormDensity(SegFileList);

