function [iters]= KL_Solver_Mod (matrix, directory)
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
%         (imported 3d profile from COMSOL) 
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
%       What are the indices?
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

a_eff = a_0*eps_r/2  ; %a_eff in nm
E_eff = 4*E_0/eps_r^2 ;%E_eff in meV

% Distance units scaling
Delta=2;
Escale = Delta*Delta;
%%
% Problem parameters
% ========================================================================
    dirstr = directory;
    if (exist(dirstr,'dir'))
    else
        mkdir(dirstr);
    end


V_imported=DataProcessing(1,0,0,0,0,directory);   % import the potential distribution using the function created. Indices indicate desired gate voltages.
[L,L,Lz]=size(V_imported);              % points along x,y , and z given by the size of imported 3x3 matrix
x=1:1:L;
y=1:1:L;
z=1:1:Lz;
% L  = 70;               % # grid points along x and y
% Lz = 18;               % # grid points along z: will be returned when the DataProcessing code is called
Bx   = 0 ;               % x-component of magnetic field in Tesla
By   = 0.1 ;             % y-component of magnetic field in Tesla
Bz   = 0 ;               % z-component of magnetic field in Tesla
% Ez   = (0:1:4)*1e7;    % in V/m: can be scrapped once we are importing our own potential
a_qd = 8;                % diameter in nm of the quantum dot(for the harmonic potential)

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
% omega_p=hbar/(m_p*(a_qd*1e-9)^2);
hx=Bx*g*mu_B/(-q_e*E_eff/(Escale*1e3));
hy=By*g*mu_B/(-q_e*E_eff/(Escale*1e3));
hz=Bz*g*mu_B/(-q_e*E_eff/(Escale*1e3));
dlh=q_e*b*(2*exx-2*ezz)/(-q_e*E_eff/(Escale*1e3));
%%
% Beginning of code
% =================
LL = L*L;
N = LL*Lz;

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


% dirstr=datestr(now,'yyyymmddHHMMSS');
%%
% Potential energy
[rows,~]=size(matrix);
iters=rows;
for d = 1:iters
    [V_imported]=DataProcessing(matrix(d,1),matrix(d,2),matrix(d,3),matrix(d,4),matrix(d,5),directory);   % import the potential distribution using the function created. Indices indicate desired gate voltages
    kU = 0;
    disp("Iteration number:", num2str(d));
    for i = 1:N
        x = 1+floor(mod(i-1,L));
        y = 1+floor(mod(i-1,LL)/L);
        z = 1+floor((i-1)/LL);
        
        % How to get energy in eV /Escale * E_eff*1000
        kU = kU + 1;
        iU(kU) = i; jU(kU) = i;
        vU(kU) = V_imported(x,y,z)./(-q_e*E_eff/(Escale*1e3));

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
    
	% in-plane Zeeman energy
	hp0 = hx-1i*hy;
	hp1 = sqrt(3)/2*hp0;
	[ihp,jhp,vhp] = find(hp0);
	
	[iShp,jShp,vShp] = find(S+hp0);
	
    %    #1,1      #1,2      #1,3      #2,1           #2,2      #2,3      #2,4      #3,1      #3,2       #3,3      #3,4      #4,2      #4,3          #4,4   
    iH = [iHHp;    iShp;     iR;       N+jShp;        N+iLHp;   N+ihp;    N+iR;     2*N+jR;   2*N+jhp;   2*N+iLHm; 2*N+iShp; 3*N+jR;   3*N+jShp;     3*N+iHHm];
    jH = [jHHp;    N+jShp;   2*N+jR;   iShp;          N+jLHp;   2*N+jhp;  3*N+jR;   iR;       N+ihp;     2*N+jLHm; 3*N+jShp; N+iR;     2*N+iShp;     3*N+jHHm];
    vH = [vHHp;    -vS+hp1;  vR;       conj(-vS+hp1); vLHp;     hp0;      vR;       conj(vR); conj(hp0); vLHm;     vS+hp1;   conj(vR); conj(vS+hp1); vHHm];

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
    

    
    voltages= [matrix(d,1) matrix(d,2) matrix(d,3) matrix(d,4) matrix(d,5)];
    save(sprintf('%s/VE%d',dirstr,d),...
        'D','Delta','E_0','E_eff','E_n','E_sorted','Escale',...
        'K','L','LL','Lz','N','V_sorted',...
        'a_qd','a_0','a_eff','dirstr','eps_0','eps_r',...
        'gamma_1','gamma_2','gamma_3','hbar','jeV','kP','kQ','kR',...
        'kRe','kS','kU','m_e','n','opts','q_e',...
        'x','y','z',...
        'H','voltages',...
        '-v7.3');
    
     
end
end
