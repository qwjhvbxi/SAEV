% works on Windows

function HostName=getHost()

% [~,HostName]=system('hostname');
[~,HostName]=unix('hostname');
HostName=strtrim(HostName);