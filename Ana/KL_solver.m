clear all;
close all;

% How does this code work? 
% There are several sections
%
% "Fixed parameters" -> just constants, nothing to change here
%
% "Derived scales" 
%    a_eff = the default spacing between grid points in real space solution of
%            wavefunction, derived from other parameters (don't change it)
%    Delta = a way to make spacing between points smaller 
%            (you can change it to make make finer spacing, see next line)
%    a_eff/Delta = the actual spacing between grid points in real space
%                  solution of wavefunction
%
% "Problem Parameters"
%    L = # points in x and y, spaced a_eff/Delta apart (you can change)
%    Lz = # points in z, spaced a_eff/Delta apart (you can change)
%    B0 = magnetic field in Tesla (you can change)
%    Ez = electric field in V/m (you can change)
%    a_qd = radius of QD which defines the potential (you can change)
%    vU = the electrostatic potential, right now it is harmonic + vertical
%         field, but we want to change this to a realistic potential
%
% "Material Parameters"
%    e_xx = strain in material (x dir)
%    e_yy = strain in material (y dir)
%    gamma_i = Material dependent Luttinger parameters
%    g = Material dependent Lande g-factor
%    a = hydrostatic deformation potential in eV
%    b = deformation potential in eV
%    d = shear deformation potential in eV
%
% "Calculated quantities"
%    V_sorted = normalized wavefunction for K states. 
%       What are the indices? 
%       psi=squeeze(reshape(V_sorted,L,L,Lz,4,K))
%       index 1: # grid points, x
%       index 2: # grid points, y
%       index 3: # grid points, z
%       index 4: KL spinor index for m_J, 1(-3/2), 2(-1/2), 3(1/2), 4(3/2) 
%       index 5: state #
%       
%    E_sorted = energies of first K states
%       What are the induces?
%       index 1: state #


% Fixed parameters
% ========================================================================
hbar  = 1.05457148e-34;
m_e   = 9.10938188e-31;
eps_0 = 8.854187817620e-12;
q_e   = 1.60217646e-19;
jeV   = 6.24150974e18;
mu_B  = 9.274009994e-24; 

% Derived scales
% ========================================================================
a_0 = (4*pi*eps_0*hbar^2)/(m_e*q_e^2)*1e9;            %Bohr radius in nm
E_0 = -m_e/(2*hbar^2)*(q_e^2/(4*pi*eps_0))^2*jeV*1e3; %energy in meV

eps_r = 16.2;  %dimensionless

a_eff = a_0*eps_r/2 ;  %a_eff in nm
E_eff = 4*E_0/eps_r^2 ;%E_eff in meV

% Distance units scaling
Delta=2;  %originally 2

Escale = Delta*Delta;
%%
% Problem parameters
% ========================================================================

L  = 70;               % # grid points along x and y
Lz = 18;               % # grid points along z
B0   = 0.1 ;            % in Tesla
Ez   = (0:0.1:4)*1e7;  % in V/m
a_qd = 8;              % in nm (for the potential)

exx  = -0.0063; %strain
ezz  = 0.0044;  %strain

% Luttinger parameters
% gamma_1 = 4.285;
% gamma_2 = 0.339;
% gamma_3 = 1.446;
gamma_1=13.15;
gamma_2=4.4;
gamma_3=5.69;

g=20.36;

% Deformation potential
a     = 2.0;   %eV
b     = -2.16; %eV
d     = -6.06; %eV

% Calculated properties
m_p=m_e/(gamma_1+gamma_2);
omega_p=hbar/(m_p*(a_qd*1e-9)^2);
hz=B0*g*mu_B/(-q_e*E_eff/(Escale*1e3));
dlh=q_e*b*(2*exx-2*ezz)/(-q_e*E_eff/(Escale*1e3));

% Beginning of code
% =================
LL = L*L;
N = LL*Lz;

%Center of mass
x0 = round(L/2)+0.5;
y0 = x0;

Vapp = @(x,y,z,Ez)((q_e*Ez*(z*a_eff*1e-9)+...
                 1/2*m_e*omega_p^2*((x-x0).^2+(y-y0).^2)*(a_eff*1e-9).^2));

             
             %%

%How many states?
K = 4;

s3 = sqrt(3);

iP = zeros(7*N-6*LL,1);      jP = iP;vP = iP;kP = 0;
iQ = zeros(6*N-6*LL,1);      jQ = iQ;vQ = iQ;kQ = 0;
iR = zeros(8*N-12*LL+4*L,1); jR = iR;vR = iR;kR = 0;
iS = zeros(8*N-16*LL+8*L,1); jS = iS;vS = iS;kS = 0;
iU = zeros(1,N);             jU = iU;vU = iU;kU = 0;
iRe = zeros(1,N);            jRe = iRe;vRe = iRe;kRe = 0;

% Kinetic energy operator
tic;
for i = 1:N
    x = 1+floor(mod(i-1,L));
    y = 1+floor(mod(i-1,LL)/L);
    z = 1+floor((i-1)/LL);
    
    kP = kP +1;
    iP(kP) = i; jP(kP) = i; vP(kP) = 6*gamma_1;
    
    if (x ~= 1)
        kP = kP + 1;
        iP(kP) = i; jP(kP) = i-1; vP(kP) = -gamma_1;
        
        kQ = kQ + 1;
        iQ(kQ) = i; jQ(kQ) = i-1; vQ(kQ) = -gamma_2;
        
        kR = kR + 1;
        iR(kR) = i; jR(kR) = i-1; vR(kR) = gamma_2*s3;
        
        if (y ~= 1)
            kR = kR + 1;
            iR(kR) = i; jR(kR) = i-1-L; vR(kR) = -gamma_3*1i*s3*0.5;
        end
        if (y ~= L)
            kR = kR + 1;
            iR(kR) = i; jR(kR) = i-1+L; vR(kR) = gamma_3*1i*s3*0.5;
        end
        if (z ~= 1)
            kS = kS + 1;
            iS(kS) = i; jS(kS) = i-1-LL; vS(kS) = -gamma_3*s3*0.5;
        end
        if (z ~= Lz)
            kS = kS + 1;
            iS(kS) = i; jS(kS) = i-1+LL; vS(kS) = gamma_3*s3*0.5;
        end
    end
    if (x ~= L)     
        kP = kP + 1;
        iP(kP) = i; jP(kP) = i+1; vP(kP) = -gamma_1;
        
        kQ = kQ + 1;
        iQ(kQ) = i; jQ(kQ) = i+1; vQ(kQ) = -gamma_2;
        
        kR = kR + 1;
        iR(kR) = i; jR(kR) = i+1; vR(kR) = gamma_2*s3;
        
        if (y ~= 1)
            kR = kR + 1;
            iR(kR) = i; jR(kR) = i+1-L; vR(kR) = gamma_3*1i*s3*0.5;
        end
        if (y ~= L)
            kR = kR + 1;
            iR(kR) = i; jR(kR) = i+1+L; vR(kR) = -gamma_3*1i*s3*0.5;
        end
        if (z ~= 1)
            kS = kS + 1;
            iS(kS) = i; jS(kS) = i+1-LL; vS(kS) = gamma_3*s3*0.5;
        end
        if (z ~= Lz)
            kS = kS + 1;
            iS(kS) = i; jS(kS) = i+1+LL; vS(kS) = -gamma_3*s3*0.5;
        end
    end
    if (y ~= 1)
        kP = kP + 1;
        iP(kP) = i; jP(kP) = i-L; vP(kP) = -gamma_1;
        
        kQ = kQ + 1;
        iQ(kQ) = i; jQ(kQ) = i-L; vQ(kQ) = -gamma_2;
        
        kR = kR + 1;
        iR(kR) = i; jR(kR) = i-L; vR(kR) = -gamma_2*s3;
        
        if (z ~= 1)
            kS = kS + 1;
            iS(kS) = i; jS(kS) = i-L-LL; vS(kS) = gamma_3*1i*s3*0.5;
        end
        if (z ~= Lz)
            kS = kS + 1;
            iS(kS) = i; jS(kS) = i-L+LL; vS(kS) = -gamma_3*1i*s3*0.5;
        end
    end
    if (y ~= L)
        kP = kP + 1;
        iP(kP) = i; jP(kP) = i+L; vP(kP) = -gamma_1;
        
        kQ = kQ + 1;
        iQ(kQ) = i; jQ(kQ) = i+L; vQ(kQ) = -gamma_2;
        
        kR = kR + 1;
        iR(kR) = i; jR(kR) = i+L; vR(kR) = -gamma_2*s3;
        
        if (z ~= 1)
            kS = kS + 1;
            iS(kS) = i; jS(kS) = i+L-LL; vS(kS) = -gamma_3*1i*s3*0.5;
        end
        if (z ~= Lz)
            kS = kS + 1;
            iS(kS) = i; jS(kS) = i+L+LL; vS(kS) = gamma_3*1i*s3*0.5;
        end
    end
    if (z ~= 1)
        kP = kP + 1;
        iP(kP) = i; jP(kP) = i-LL; vP(kP)  = -gamma_1;
        
        kQ = kQ + 1;
        iQ(kQ) = i; jQ(kQ) = i-LL; vQ(kQ) = 2*gamma_2;
    end
    if (z ~= Lz)
        kP = kP + 1;
        iP(kP) = i; jP(kP) = i+LL; vP(kP)  = -gamma_1;
        
        kQ = kQ + 1;
        iQ(kQ) = i; jQ(kQ) = i+LL; vQ(kQ) = 2*gamma_2;
    end
end

disp('***********CPU Timing Results in seconds***********');
disp(['P,Q,R,S,U construction:                      ' num2str(toc)]);

V0 = rand(4*N,1);

E = zeros(K,length(Ez));

dirstr=datestr(now,'yyyymmddHHMMSS');
pijx_final=zeros(1,length(Ez));
pijy_final=zeros(1,length(Ez));
pijz_final=zeros(1,length(Ez));
% Potential energy
for d = 1:length(Ez)
    kU = 0;
    for i = 1:N
        x = 1+floor(mod(i-1,L));
        y = 1+floor(mod(i-1,LL)/L);
        z = 1+floor((i-1)/LL);
        
        % How to get energy in eV /Escale * E_eff*1000
        kU = kU + 1;
        iU(kU) = i; jU(kU) = i;
        vU(kU) = Vapp(x,y,z,Ez(d))./(-q_e*E_eff/(Escale*1e3));

        kRe = kRe + 1;
        iRe(kRe) = i; jRe(kRe) = i;
        if z <= 2
            vRe(kRe) = 0*sqrt(-1);
        else
            vRe(kRe) = 0;
        end
    end
    
    P  = sparse(iP,jP,vP);
    Q  = sparse(iQ,jQ,vQ);
    Rk = sparse(iR,jR,vR);
    Re = sparse(iRe,jRe,vRe);
    S  = sparse(iS,jS,vS);
    U  = sparse(iU,jU,vU);

    % On diagonal terms
    
    HHp = P + Q + U + (+3/2*hz-dlh/2)*speye(N,N);
    LHp = P - Q + U + (+1/2*hz+dlh/2)*speye(N,N);
    LHm = P - Q + U + (-1/2*hz+dlh/2)*speye(N,N);
    HHm = P + Q + U + (-3/2*hz-dlh/2)*speye(N,N);
    
    R = Rk + Re;

    [iHHp,jHHp,vHHp] = find(HHp);
    [iLHp,jLHp,vLHp] = find(LHp);
    [iLHm,jLHm,vLHm] = find(LHm);
    [iHHm,jHHm,vHHm] = find(HHm);
    [iR,jR,vR] = find(R);
    
    %    #1,1      #1,2      #1,3      #2,1      #2,2      #2,4      #3,1      #3,3      #3,4      #4,2      #4,3      #4,4   
    iH = [iHHp;    iS;       iR;       N+jS;     N+iLHp;   N+iR;     2*N+jR;   2*N+iLHm; 2*N+iS;   3*N+jR;   3*N+jS;   3*N+iHHm];
    jH = [jHHp;    N+jS;     2*N+jR;   iS;       N+jLHp;   3*N+jR;   iR;       2*N+jLHm; 3*N+jS;   N+iR;     2*N+iS;   3*N+jHHm];
    vH = [vHHp;    -vS;      vR;      -conj(vS); vLHp;     vR;       conj(vR); vLHm;     vS;       conj(vR); conj(vS); vHHm];

    H = sparse(iH,jH,vH);

    clear P;
    clear Q;
    clear R;
    clear S;
    clear U;
    
    opts.p = 60;
    opts.issym = 0;
    opts.isreal = 0;
    opts.disp = 0;
    opts.maxit = 10000;
    %opts.v0 = V0;
    
    tic;
    [V,D] = eigs(H,K,'sr',opts);
    [E_n,n] = sort(diag(D));
    % Calculated quantities
    % ========================================================================
    V_sorted = V(:,n);
    E_sorted = -E_n/Escale*E_eff; %answer in meV
    clear V;
    
    disp(['Found eigenvalues of H:                      ' num2str(toc)]);
    
    if (exist(dirstr,'dir'))
    else
        mkdir(dirstr);
    end
    
    save(sprintf('%s/VE%d',dirstr,d),...
        'D','Delta','E_0','E_eff','E_n','Ez','E_sorted','Escale',...
        'K','L','LL','Lz','N','V_sorted',...
        'a_qd','a_0','a_eff','dirstr','eps_0','eps_r',...
        'gamma_1','gamma_2','gamma_3','hbar','jeV','kP','kQ','kR',...
        'kRe','kS','kU','m_e','n','opts','q_e',...
        'x','x0','y','y0','z',...
        'H',...
        '-v7.3');

  

end

%% Plotting the potential 

x=linspace(0,70,1000);
y=linspace(0,70,1000);
[X,Y]= meshgrid(x,y);
z=1;
% plot the potential for Ez=0
V_app=((q_e*Ez(end)*(z*a_eff*1e-9)+1/2*m_e*omega_p^2*((X-x0).^2+(Y-y0).^2)*(a_eff*1e-9).^2));
mesh(X,Y,V_app)

