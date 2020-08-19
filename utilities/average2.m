%% [av] = AVERAGE2(t,n)
% averages matrix t every n values
% 

function av=average2(t,n)

if (rem(n,1)==0)
    
    if (n>0)
        
        s=size(t);
        
        if (rem(s(1),n)==0)
            
            av=squeeze(mean(reshape(t,n,s(1)/n,s(2))));

            if s(2)==1 
                
                av=av';
                
            end
            
        else
            
            error('The number of rows of matrix t has to be a multiple of n!');
            
        end
        
    else
        
        error('n must be positive!');
        
    end
else
    
    error('n must be an integer!');
    
end

end