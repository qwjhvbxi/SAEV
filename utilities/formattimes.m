%% a=FORMATTIMES(seconds)
% Return string with formatted times in the form m:ss or h:mm:ss

function a=formattimes(seconds)
    if seconds<3600
        a=sprintf('%d:%0.2d',floor(seconds/60),round(rem(seconds,60)));
    else
        a=sprintf('%dh%0.2d:%0.2d',floor(seconds/3600),round(rem(seconds/60,60)),round(rem(seconds,60)));
    end
end