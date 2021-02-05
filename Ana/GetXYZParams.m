%% Version
function[set_x,set_y,set_z] = GetXYZParams(filename_str)

if nargin < 1 
    filename_str="D:\Lab\Salfi\Electrostatics_3D_Comsol\COMSOL_Potential.txt";
end

unsorted=readmatrix(filename_str);
[a,b]=size(unsorted);    % find size of imported matrix


set_x = unique(unsorted(1:a,1));
set_y = unique(unsorted(1:a,2));
set_z = unique(unsorted(1:a,3));

end