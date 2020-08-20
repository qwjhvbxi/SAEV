function TimeLeft=displayprogress(iteration,totaliterations,timeelapsed)

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

DisplayResolution=20;
TimeLeft=timespent/relatedtime*(totaliterations-iteration);

fprintf('\n\n');
fprintf('Progress: %0.1f%%\n',iteration/totaliterations*100);
displaytimes(TimeLeft,'ETA');
fprintf(char('*'*ones(1,floor(iteration/totaliterations*DisplayResolution)  )));
fprintf(char('-'*ones(1,DisplayResolution-floor(iteration/totaliterations*DisplayResolution) )));
fprintf('\n\n');

end