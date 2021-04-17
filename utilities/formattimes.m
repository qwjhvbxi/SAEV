function a=formattimes(seconds)
    a=sprintf('%d:%0.2d',floor(seconds/60),round(rem(seconds,60)));
end