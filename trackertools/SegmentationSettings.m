% 
% Define "Settings" structure containing all behavior segmentation
% settings. Open file to view definitions and [defaults]
%

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 01-Dec-2009 15:01:57 
%---------------------------- 

Settings = struct(  'StallDistance', 0.022, ...     % Max instantaneous speed for stall (mm/s) [0.022]
                    'FrameRate', 2, ...             % Video framerate (frames/s) [2]
                    'PixelSize', 36, ...            % Spatial calibration (pixels/mm) [36]
                    'FixJitter', true, ...          % obsolete
                    'SmoothWinSize', 3, ...         % Size of Window for smoothing track data (in frames) [3]
                    'MaxBodyPathAngleDev', 13, ...  % Max deviation (deg) of path vs body orientation in Fwd Run [13]
                    'MaxPathAngleDev', 11, ...      % Max deviation (deg) of path from 60deg for Fwd Run [11]
                    'MaxFwdRunEcc', 0.4, ...        % Max eccentricity for Fwd Run [0.4]
                    'MinFwdRunFr', 8, ...           % Min consecutive frames for 'real' Fwd Run [8]
                    'MinRevAngVel', 259, ...        % (was 210) Min angular velocity (deg/s) for reversal [259]
                    'MaxRevTime', 13, ...           % Max consecutive reversal time (s) [13]
                    'MaxShortTurnLen', 3, ...       % Max short turn (s) [3]
                    'MaxSmoothAngVel', 120, ...     % Maximum angular velocity for smooth, forward turns (deg/s) [120]
                    'MinPirEcc', 0.4, ...           % Min eccentricity for pirouette 'curl' state [0.4]
                    'MinOmegaEcc', 0.72, ...        % Min eccentricity for clear omega turn [0.72]
                    'Omega2PostMinEcc', 0.66, ...   % Min eccentricity for omega around 2 posts [0.66]
                    'Omega2PostTestTime', 3.5, ...  % Test time for angular velocity after omega (s) [3.5]
                    'Omega2PostMinAngVel', 18, ...  % Minimum angular velocity average after omega (deg/s) [18]
                    'MinOmegaInterval', 3, ...      % Minimum interval between omegas (s) [3]
                    'MinPirEccTime', 2, ...         % Min time of eccentric 'curl' states (s) [2]
                    'MinErraticAngAcc', 200, ...    % Min angular acceleration for erratic movement (or swim) (deg/s^2) [200]
                    'MinErrTime', 6, ...            % Min time of erratic movements (s) [6]
                    'ReversalHistory', 15, ...      % Time history (s) of reversals before omega or gradual curve [15]
                    'LongReversalTime', 4, ...      % MinTime for long reversal [4]
                    'XBorderPadding', 0, ...        % Distance (pix) to shrink the selected tracking boundary [0]
                    'XBorderEffectTime', 2, ...     % Time (s) after crossing X-border to exclude [2]
                    'CollisionRelSize', 1.5, ...    % relative size threshold for collision detection [1.5]
                    'CollEffectTime', 5, ...        % Time (s) to ignore data after collision detection [5]
                    'PriorRevTestTime', 15, ...     % Time (s) to test for reversals prior to omega/turn [15]
                    'MinPriorRevTime', 4, ...       % Minimum time (s) of reversals prior to omega within test time [4]
                    'BehDef', [ 1 1 1 1 7 2 6 7 4 8;  ... % Behavior code: 1-F, 2-LF, 3-R, 4-P, 5-OmR, 6-Omf, 7-?, 8-OOB
                                3 3 3 3 7 3 5 7 4 8;  ...
                                4 4 4 4 4 4 4 4 4 8 ]);
                