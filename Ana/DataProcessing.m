function[sum] = DataProcessing(val1 ,val2 ,val3, val4, val5,dirstr)
unsorted=readmatrix("D:\Lab\Salfi\Electrostatics_3D_Comsol\QDElectrostaticsNew(terminalmod)");
[a,b]=size(unsorted);    % find size of imported matrix
numgates= b-3;           % usually 5 gates

set_x = unique(unsorted(1:a,1));
set_y = unique(unsorted(1:a,2));
set_z = unique(unsorted(1:a,3));
delta_x=set_x(2)-set_x(1);
delta_z=set_z(2)-set_z(1);

xlen=length(set_x);
ylen=length(set_y);
zlen=length(set_z);


% create an empty array 
gatearr= zeros(xlen,ylen,zlen);


%%
for n=1:numgates
    rowindex=1;
    gateindex=n+3;
    gateid=strcat('gate',num2str(n));
    for z=1:zlen
        for y=1:ylen
            for x=1:xlen
                gatearr(y,x,z)=unsorted(rowindex,gateindex)*-1.602e-19;
                rowindex=rowindex+1;

            end

        end

    end 
    gates.(gateid)=gatearr;

end
%%
% figure
% for i=1:1:numgates
%     gateid=strcat('gate',num2str(i));
%     subplot(2,3,i)
%     surf(-1*gates.(gateid)(:,:,zlen))
%     view(2)
%     shading interp
%     axis tight
% end
% hold on
% sgtitle("Gates 1 to 5, and their sum")

%%
sum=zeros(xlen,ylen,zlen);
matrix= [val1 val2 val3 val4 val5];
for n=1:1:numgates
    gateid=strcat('gate',num2str(n));
    gates.(gateid)=gates.(gateid)*matrix(n);
    sum=sum+gates.(gateid);
    grad=sum./1.602e-19;
    [E_fieldx, E_fieldy, E_fieldz]=gradient(grad,delta_x,delta_x,delta_z);
end

% norm=sqrt(E_fieldx.^2+E_fieldy.^2+E_fieldz.^2);

% subplot(2,3,6)
% figure;
% surf(sum(:,:,zlen));
% title("Potential Distribution")
% xlabel("points in x")
% ylabel("points in y")
% zlabel("Potential Energy(J)")
% % view(2)
% shading interp
% axis tight
% % 
matrix_string=strjoin(string(matrix));
str=sprintf('%s/E_field [%s].mat',dirstr,matrix_string);
save(str,'E_fieldx', 'E_fieldy','E_fieldz');
title(matrix_string)
end




  




