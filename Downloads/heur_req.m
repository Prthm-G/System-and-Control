%% check if step info requirements are met.
% inputs: stepinfo response structure
%       param = 'overshoot','risetime', etc
%       overunder = 'over' or 'under'
%       val = value you are checking (double)
% OUTPUTS: boolean true or false. 
% Troubleshooting. Check to make sure that the name of param is correct.
% idk I haven't checked it long enough.
% Eg.   S = stepinfo(sys);
%       fit_req = heur_req(S,'overshoot','under',20);
% For ELEC 341, 2022
% Date: December 14, 2022
function H = heur_req(info_struct,param,overunder,val)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
H = 0;
s = tf('s');
names = fieldnames(stepinfo(1/s));
order = contains(names,param,'IgnoreCase',true);
index = find(order==1);
if(index)
end

% If over is defined, use the over test. (what about undefined?)
if(strcmpi('over',overunder))
    if(info_struct.(names{index}) > val)
        H = 1;
        return;
    else
        H = 0; % fails condition
        return;
    end
elseif (strcmpi('under',overunder))
    if(info_struct.(names{index}) < val)
        H = 1;
        return;
    else
        H = 0;
        return;
    end
else
    error('Invalid Condition Call');
end
end