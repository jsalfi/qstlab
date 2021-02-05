%% Version
N=41;
qs=zeros(1,N);

for d=1:1:N
    dirstr='D:\Lab\Salfi\KLVb_3D_Matlab\Results\20201029213713\';
    fn=[dirstr,'\VE',num2str(d),'.mat'];
    load(fn);
    
    qs(d)=E_sorted(2)-E_sorted(1);
    
    psi=squeeze(reshape(V_sorted,L,L,Lz,4,K));
    nzy1=squeeze(sum(abs(psi(L/2,:,:,:,1)).^2,4));
% 
%     figure;
%     imagesc((1:1:Lz)*a_eff*Delta,(1:1:L)*a_eff*Delta,nzy1);
%     xlabel('z(nm)');
%     ylabel('x or y(nm)');
%     colorbar;
%     axis image;
%     saveas(gcf,['image',num2str(d),'.pdf']);
%     close 
       %% This section calculates the dipole operator pij
psi=squeeze(reshape(V_sorted,L,L,Lz,4,K));   % turn the wavefunction array into something properly indexed
pij_mj_x=zeros(K,K,K);    % a 4x4x4 matrix to hold the ij and mj values for the x dipole operator
pij_mj_y=zeros(K,K,K);    % a 4x4x4 matrix to hold the ij and mj values for the y dipole operator
pij_mj_z=zeros(K,K,K);    % a 4x4x4 matrix to hold the ij and mj values for the z dipole operator

x=(1:1:L).*(a_eff*1e-9);
y=(1:1:L).*(a_eff*1e-9);
z=(1:1:Lz).*(a_eff*1e-9);
[X, Y ,Z]=meshgrid(x,y,z);

% Create the pij matrix (an intermediate form that depends on mj still)
for i=1:4
    for j=1:4
        for mj=1:4
            pij_intermediate_x=psi(:,:,:,mj,i).*conj(psi(:,:,:,mj,j)).*X*q_e;
            pij_intermediate_y=psi(:,:,:,mj,i).*conj(psi(:,:,:,mj,j)).*Y*q_e;
            pij_intermediate_z=psi(:,:,:,mj,i).*conj(psi(:,:,:,mj,j)).*Z*q_e;
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

end
% 
hbar  = 1.05457148e-34;
h=hbar*2*pi;
time=h./(2*10^5*abs(pijx_final));


%%

figure;
plot(Ez,qs*1000);
xlabel('Field (V/m)');
ylabel('Qubit splitting (\mu eV)');