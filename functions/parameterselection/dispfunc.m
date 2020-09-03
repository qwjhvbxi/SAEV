
function [state,options,optchanged] = dispfunc(options,state,flag)
    optchanged=0;
    [~,II]=min(state.Score);
    state.Population(II,:)
end