%% Tune the PID controller master gain to achieve a specific phase margin
% Step 7 of the 10 step Controller Design Process
% Return the controller dynamic parameter master gain
% Input is just an LTI system which has been tuned to it's maximum Phase
% Margin and the target phase margin.

function K = masterGainTune(sys,PM_targ)
% Tunes the PID controller master gain for a specific phase margin
%   Detailed explanation goes here
prop = struct('orig',1,'up',[],'down',[]);
dev = struct('orig',[],'up',[],'down',[]);

if(~test_margin(sys*prop.orig))
    % begin
    [Gm Pm] = margin(sys*prop.orig);
    dev.orig = abs(PM_targ-Pm);
else
    error('origin kref invalid');
end
gRain = 1;


while(gRain<15)
[prop.up Pm] = searchPM(prop.orig,sys,gRain,'up');
dev.up = abs(Pm-PM_targ);
[prop.down Pm] = searchPM(prop.orig,sys,gRain,'down');
dev.down = abs(Pm-PM_targ);

dev_arr = [dev.orig dev.up dev.down];
[M,dir] = min(dev_arr);

switch dir
    case 1
        % Origin reached, increasing grain.
        gRain = gRain+1;
    case 2
        prop.orig = prop.up;
        dev.orig = dev.up;
        gRain = gRain-1;
    case 3
        prop.orig = prop.down;
        dev.orig = dev.down;
        gRain = gRain-1;
end
%disp(['Proportion is: ',num2str(prop.orig), 'dev is: ',num2str(dev.orig)]);

% Reset structures
prop.up = [];
dev.up = [];
prop.down = [];
dev.down = [];

end
K = prop.orig;
return;
end

%% search for phase margin, return a phase margin after a little variation
function [prop_out PM] = searchPM(prop0,sys,gRain,dir)
c =1;
while(1)
    switch dir
        case 'up'
            prop = prop0+c*10^(-gRain);
        case 'down'
            prop = prop0-c*10^(-gRain);
    end

    if(~test_margin(prop*sys))
        [Gm Pm] = margin(prop*sys);
        PM = Pm;
        prop_out = prop;
        break;
    else
        c = c+1;
    end
end

end

%% Test if margin() works
% internal flag. test_margin checks that the margin() function works. For
% some esoteric reason, margin() can fail randomly. It is essential to
% catch the error or else it will halt the program. smh.
% INPUT: sys (an open loop LTI system)
% OUTPUT: Boolean (0 or 1) if the margin() function can activate
%       1 means it can't
%       0 means it can
function C = test_margin(sys)
C = 0;
try
    bigboi = allmargin(sys);
catch E
    % Margin is unable to be initiated, this is a debugging piece of
    % code output
    %disp('Cannot initiate margin()');
    C = 1;
    % continue
end
end