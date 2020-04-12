function reply = hasParallelToolbox(S)
% Get availability of the Parallel Processing Toolbox:
%   R = hasParallelToolbox()
% Disable the flag for an existing PPT manually:
%   hasParallelToolbox(false)
% Re-enable the actual status:
%   hasParallelToolbox(true)   % Considers a missing toolbox also
%   or
%   clear('hasParallelToolbox')
persistent hasPPT usePPT
if isempty(hasPPT)
  try
    % One of these 2 tests is sufficient already:
    hasPPT = ~isempty(getCurrentTask()) && (matlabpool('size') > 0);
  catch ME
    if ~strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
      rethrow(ME);
    end
    hasPPT = false;
  end
  usePPT = hasPPT;
end
if nargin > 0
   usePPT = (any(S)) && hasPPT;
end
reply = usePPT;
end