%% HostName=gethostname()
% works on Windows

function HostName=gethostname()

% [~,HostName]=system('hostname');
[~,HostName]=unix('hostname');
HostName=strtrim(HostName);