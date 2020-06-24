%% sensitivity analysis

addpath plots functions utilities

[P,R]=generateplotline3('NYC2016',[],'Operations.maxwait',[Inf],'m',13000,'TransportLayer.ts',10:16,'TransportLayer.tr',6:18);



[~,Rw]=generateplotline3('NYC2016','peakwait','Operations.maxwait',[Inf],'m',13000,'TransportLayer.ts',10:16,'TransportLayer.tr',6:18);

% best values: ts = 12-13; tr = 7-10 (24 to 26 minutes; 14 to 20 minutes)
% chosen values: ts=12, tr=10
