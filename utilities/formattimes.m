%% a=FORMATTIMES(seconds)
% Return string with formatted times in the form m:ss

function a=formattimes(seconds)
    a=sprintf('%d:%0.2d',floor(seconds/60),round(rem(seconds,60)));
end