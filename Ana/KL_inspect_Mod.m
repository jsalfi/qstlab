% INPUTS: number of iterations of the code required (based on size of voltage matrix), save directory.
% OUTPUTS: Inspect_results.mat file with all calculated quantities

function []= KL_inspect_Mod(iters,directory)
N=iters;
qs=zeros(1,N);

% calls function that extracts COMSOL grid dimensions
[x,y,z] = GetXYZParams();
[X, Y ,Z]=meshgrid(x,y,z);


cd 'D:\Lab\Salfi\KLVb_3D_Matlab\Results\Reduced_Gates'
for d=1:1:N
    
    %% This section loads the KL solver results as VE [d] and the electric field files as E_field [d]
    fn = sprintf('%s/VE%d',directory,d);   
    load (fn);

    str=sprintf('%s/E_field [%s].mat',directory, strjoin(string(voltages)));
    load(str);
    middle_values(d)=voltages(:,1);
    
    %% This section finds the wavefunction + zeeman splitting from the KL Solver imported files
    qs(d)=E_sorted(2)-E_sorted(1);                 % find the Zeeman splitting
    psi=squeeze(reshape(V_sorted,L,L,Lz,4,K));     % turn the wavefunction array into something properly indexed
    nzy1=squeeze(sum(abs(psi(L/2,:,:,:,1)).^2,4)); % find the hole PDF
    nzyid=strcat('nyz_',num2str(d));
    hdensity.(nzyid)=nzy1;
    
    % This section plots and saves the hole pdf. I suggest you don't run
    % this code. 
%     figure;
%     imagesc((1:1:Lz)*a_eff*Delta,(1:1:L)*a_eff*Delta,nzy1);
%     xlabel('z(nm)');
%     ylabel('y(nm)');
%     colorbar;
%     axis image;
%     saveas(gcf,['image',num2str(d),'.pdf']);
%     close 


%% This section calculates the dipole operator pij_x,y,z

pij_mj_x=zeros(K,K,K);    % a 4x4x4 matrix to hold the ij and mj values for the x dipole operator
pij_mj_y=zeros(K,K,K);    % a 4x4x4 matrix to hold the ij and mj values for the y dipole operator
pij_mj_z=zeros(K,K,K);    % a 4x4x4 matrix to hold the ij and mj values for the z dipole operator

% Create the pij matrix (an intermediate form that depends on mj still)
for i=1:4
    for j=1:4
        for mj=1:4
            pij_intermediate_x=psi(:,:,:,mj,i).*X.*conj(psi(:,:,:,mj,j))*q_e;
            pij_intermediate_y=psi(:,:,:,mj,i).*Y.*conj(psi(:,:,:,mj,j))*q_e;
            pij_intermediate_z=psi(:,:,:,mj,i).*Z.*conj(psi(:,:,:,mj,j))*q_e;

            pij_mj_x(i,j,mj)= sum(pij_intermediate_x,[1,2,3]);
            pij_mj_y(i,j,mj)= sum(pij_intermediate_y,[1,2,3]);
            pij_mj_z(i,j,mj)= sum(pij_intermediate_z,[1,2,3]); 

        end
    end
end

% Now sum over the 3rd dimension, which means summing over mj
pij_x=sum(pij_mj_x,3);
pij_y=sum(pij_mj_y,3);
pij_z=sum(pij_mj_z,3);

pijx_final(d)=pij_x(1,2);
pijy_final(d)=pij_y(1,2);
pijz_final(d)=pij_z(2,2)-pij_z(1,1);

clear fn
clear str

%% This section calculates the average electric field value for Ex, Ey, Ez

for mj=1:4
    psi_int_x=conj(psi(:,:,:,mj,1)).*E_fieldx.*psi(:,:,:,mj,1);
    psi_new_x(mj)=sum(psi_int_x,[1,2,3]);
    
    psi_int_y=conj(psi(:,:,:,mj,1)).*E_fieldy.*psi(:,:,:,mj,1);
    psi_new_y(mj)=sum(psi_int_y,[1,2,3]);
    
    psi_int_z=conj(psi(:,:,:,mj,1)).*E_fieldz.*psi(:,:,:,mj,1);
    psi_new_z(mj)=sum(psi_int_z,[1,2,3]);
end
E_avg_x(d)=sum(psi_new_x);
E_avg_y(d)=sum(psi_new_y);
E_avg_z(d)=sum(psi_new_z);


%% This section calculates the QD radius 

r_hat= X+Y+Z;

for mj=1:4
    psi_R=conj(psi(:,:,:,mj,1)).*r_hat.*psi(:,:,:,mj,1);
    psi_new_R(mj)=sum(psi_R,[1,2,3]);
end
exp_R(d)=sum(psi_new_R);

for mj=1:4
    psi_R2=conj(psi(:,:,:,mj,1)).*(r_hat-exp_R(d)).*(r_hat-exp_R(d)).*psi(:,:,:,mj,1);
    psi_new_R2(mj)=sum(psi_R2,[1,2,3]);
end
SD_r(d)=abs(sqrt(sum(psi_new_R2)));


%% This section saves the wavefunction for later use.
dir = sprintf('%s/psi_%d.mat',directory,d); 
save(dir,'psi')

end


%% This section calculates the EDSR Rabi time from the pij matrix elements
hbar  = 1.05457148e-34;
h=hbar*2*pi;
time=h./(2*10^5*abs(pijx_final));

%% This section saves all calcualted quantities to a file, for later use
str=sprintf('%s/Inspect_Results.mat',dirstr);
save(str,'qs','SD_r','E_avg_x','E_avg_y','E_avg_z','pijx_final','pijy_final','pijz_final','time','middle_values')



end



