function[]=main(matrix_values, directory)
cd D:\Lab\Salfi\KLVb_3D_Matlab\Results
% note that the matrix fomrmat is: 
%
% [ PlungerL PlungerR BarrierL BarrierR BarrierCenter]
%
%
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

% Occupation_Number(iters,directory)

end

