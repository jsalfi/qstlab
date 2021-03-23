function[]=main(matrix_values, directory)
% This function runs multiple files, taking as input a save directory and a
% matrix of values that represent the voltages you want to apply on each of
% the five gates.

% note that the matrix format is: 
%
% [ PlungerL PlungerR BarrierL BarrierR BarrierCenter]
%
% The code will execute, in sequence:
% 1. KL Solver.m, which in turn will call DataProcessing.m and GetXYZParams.m to extract the potental profile from a COMSOL file. 
% 2. KL_Inspect_Mod which calculates important device parameters using the reuslts of KL_Solver.m 
% 
% Then, you must manually run the MCI code in Python, inputting your save directory so that the results can be extracted in the code.
% You can then run Occupation_Number.py.
% 
% In terms of file structure, the results will be saved as follows:
%    ----> yyyymmddHHMMSS (Master)
%       [----> PythonResults
%

% Input the folder that you will use to store your results
cd D:\Lab\Salfi\KLVb_3D_Matlab\Results\Reduced_Gates
%
% if we don't specify a directory, create one to write in
if nargin < 2
    directory = datestr(now,'yyyymmddHHMMSS');
end

% Run the KL solver with inputted matrix values array, and working directory. 
iters = KL_Solver_Mod(matrix_values, directory);

% Run the KL Inspect function by inputting the number of iterations (number
% of files to process) and the directory
KL_inspect_Mod(iters,directory)



end


%  xt = get(gca, 'XTick');
% set(gca, 'XTickLabel', fliplr(xt))

