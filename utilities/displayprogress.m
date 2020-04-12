function TimeLeft=displayprogress(iteration,totaliterations,timeelapsed)

if numel(timeelapsed)==1
    % global average
    timespent=timeelapsed;
    relatedtime=iteration;
else
    % consider only last 50 iterations
    timespent=timeelapsed(max(1,iteration))-timeelapsed(max(1,iteration-50+1));
    relatedtime=min(50,iteration);
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