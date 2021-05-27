%% timeLeft=DISPLAYETA(iteration,totaliterations,timeelapsed)
% Print progress and estimated time to complete simulation.
% Also returns the estimated time left.

function timeLeft=displayeta(iteration,totaliterations,timeelapsed)

LastIterations=100;

if numel(timeelapsed)==1
    % global average
    timespent=timeelapsed;
    relatedtime=iteration;
else
    % consider only last iterations
    timespent=timeelapsed(max(1,iteration))-timeelapsed(max(1,iteration-LastIterations+1));
    relatedtime=min(LastIterations,iteration);
end

timeLeft=timespent/relatedtime*(totaliterations-iteration);

fprintf('\n\n');
fprintf('Progress: %0.1f%%\n',iteration/totaliterations*100);
fprintf('Time elapsed: %s \n',formattimes(timeelapsed));
fprintf('Time left: %s \n',formattimes(timeLeft));
progressbar(iteration,totaliterations,20);
fprintf('\n\n');

end