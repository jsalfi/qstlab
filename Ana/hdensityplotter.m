%% Version
load('D:\Lab\Salfi\KLVb_3D_Matlab\Results\20201029221052\hdensity.mat');
len= 24;


for i=1:len
    id=strcat('nyz_',num2str(i));
    a=hdensity.(id);
 [maxI,Id] = max(a(:));
 [maxRow,maxCol] = ind2sub(size(a), Id);
 slice_y.(id)=a(maxRow,1:end);
 slice_z.(id)=a(1:end,maxCol);

end

