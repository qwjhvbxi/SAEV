function progressbar(iteration,totaliterations,dispresolution)

fprintf(char('*'*ones(1,floor(iteration/totaliterations*dispresolution)  )));
fprintf(char('-'*ones(1,dispresolution-floor(iteration/totaliterations*dispresolution) )));

end