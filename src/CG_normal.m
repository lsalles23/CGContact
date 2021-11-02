function [P,rz,phix,phiy,deltaz,uz] = CG_normal(profile,material,loading,options)
%CG_normal solves normal contact problem to find the pressure distribution 
%
% [P,rz,phix,phiy,deltaz,uz] =  CG_normal(profile,material,loading,option)
%
% Input data description
% profile
%        .x,.y                  : coordinates of all grid points (m)
%        .h                     : initial separation between the two surfaces (m)
% material
%         E1,E2 nu1,nu2         : Young's modulus and Poisson ratio of the material
%         H                     : hardness of the softer material - represents the
%                                 maximum allowable pressure 
% loading
%         N, Cx, Cy             : normal load, moments around x and y axis (N - N.m)
%         deltaz, phix, phiy    : rigid body motions : translation along z axis
%                                 rotations around x and y axis
% options
%        .niter                 : maximum number of iterations
%        .errlim_p              : tolerance to judge the convergence of the pressure
%        .info                  : diplay or not iteration error (=1 or =0)

%
% Ouput data description
% P                     : pressure distribution (Pa) or (Mpa)
% rz                    : residuals of the conjugate gradient - represent the 
%                         penetrations in the contact area (m)
% uz                    : elastic deflections in the vertical direction
%
% see also CG_tangential

%Copyright 2019 Imperial College London 
%authors Loic Salles Jason Armand

%profile of the contacting bodies
x = profile.x; y = profile.y; h = profile.h;
%material properties
E1 = material.E1; E2 = material.E2; nu1 = material.nu1; nu2 = material.nu2;

if isfield(material,'H')
    H = material.H; %hardness
else
    H = Inf; %pure elastic contact
end

%loading
% normal force
N = loading.N;
% Moment around x and y
Cx = loading.Cx;
Cy = loading.Cy;

deltaz = loading.deltaz;
phix = loading.phix;
phiy = loading.phiy;

%options
niter = options.niter;
errlim = options.errlim_p;
info = options.info;

%% Computation of the influence coefficients
%TODO for spline coupling
%build other compliance matrix with FFT fA
Nx=length(x); dx=x(2)-x(1); Ny=length(y); dy=y(2)-y(1);
A=conn(Nx,dx,Ny,dy,E1,nu1,E2,nu2);
fA=fft2(A);
size(fA);
clear A;

%% Pressure initialization (which verifies the equilibrium equations)
if (N~=0) || (Cx~=0) || (Cy~=0) %imposed forcing
   P=ones(Nx,Ny);
   [P,X0] = Equilibrium(P,ones(Nx,Ny),loading,profile,material);
else %imposed displacement
   P=ones(Nx,Ny);
end
%% Setting of the initial values of the CG variables
Gold=1; 
khi=0;              % variable used to reset the conjugate gradient
pk=zeros(Nx,Ny);    % descent directions
err=1;
errP=1;
it=0;

%initalise working array
PP=zeros(2*Nx,2*Ny);
PK=zeros(2*Nx,2*Ny);
rz=zeros(Nx,Ny);
Pold=zeros(Nx,Ny);
    
if info == 1
   disp('Iteration   error');
end
%main conjugate gradient loop
while ((errP>errlim)&&(it<niter))
    it = it+1;
    % active area
    %s=((P>0)&(P<H))|((P<=0)&(h<=0));
    s=(P>0)&(P<H);
     % inactive area: not in contact or limit of plasticity 
    sn=((P<=0)|(P>=H));
    PP(1:Nx,1:Ny)=P;
    % displacement due to contact pressures
    dd=real(ifft2(fA.*fft2(PP,2*Nx,2*Ny)));
    uz=dd(1:Nx,1:Ny); 
    clear dd;
    %TODO for spline coupling
    %   add displacement due to other compliance torque, bending and foundation
    %   uztot = uzcont+uzbend+uzfound
    % =====================================================================
    if (N~=0) || (Cx~=0) || (Cy~=0)
        % Estimation of the rigid body motions deltaz, phiy and phix (Eq. 3.47):
        X1 = RigidBodyMotion(uz,s,profile);
        deltaz=X1(1,1);
        phix=X1(2,1);
        phiy=X1(3,1);
    end
    % =====================================================================
    % Compute the residuals (Eq. 3.49)
    rz=zeros(Nx,Ny); 
    for i=1:Nx
        for j=1:Ny
            if (s(i,j)>0)
                rz(i,j)=h(i,j)+uz(i,j)-deltaz-phix*y(j)+phiy*x(i);
            end
        end
    end
    % =====================================================================
    % G norm of the residuals in the contact area
    G=sum(rz(s).*rz(s));
    % =====================================================================
    % Computation of the descent directions pk (Eq. 3.51)
    pk(s)=-rz(s)+khi*G/Gold*pk(s);
    pk(sn)=0;
    Gold=G;
    % =====================================================================
    % Computation of qk: deflections corresponding to the descent dir. (Eq. 3.53)
    PK(1:Nx,1:Ny)=pk;
    dd=real(ifft2(fA.*fft2(PK,2*Nx,2*Ny)));
    qk=dd(1:Nx,1:Ny); clear dd
    %TODO for spline coupling
    %   add displacement due to other compliance torque, bending and foundation
    %   qktot = qkcont+qkbend+qkfound
    % =====================================================================
    % Computation of the coefs alphak=dp: step to be made in the descent
    % directon (Eq 3.54)
    dp=sum(rz(s).*pk(s))/sum(qk(s).*pk(s));
    P(s)=P(s)-dp*pk(s); % Eq. 3.56
    % =====================================================================
    % Projection
    P=min(H,max(0,P));
   
    sol=find((P==0)&(rz<0)|(P==H)&(rz>0)); %points outside the active area where status changed
    
    if isempty(sol)
        khi=1;
    else
        P(sol)=P(sol)-dp*rz(sol); % separation negative -> the point must be transfered to the potential contact area
        % -dp*rz(sol) is positive for sure (Eq 3.58) !There is a typo in the thesis h(i,j) -> rz(i,j)
        khi=0; % re-initialization of the conjugate gradient
    end
    s=((P>0)&(P<H)); % update of the active contact area
    % =====================================================================
    % Ajustment of the pressures according to the static equilibrium (Eq. 3.60)
    %TODO spline coupling
    %modify equation of equilibrium for Torque equation
    if (N~=0) || (Cx~=0) || (Cy~=0)
        [P,X2] = Equilibrium(P,s,loading,profile,material);
    end
    % =====================================================================
    % Convergence check
    if it==1
        err =1;
        errP=1;
    end
    if it>1
        errP=abs((max(max(P))-max(max(Pold))))/max(max(Pold));
    end
    % =====================================================================
    Pold=P;
    % Error computation ? Eq. 3.64
    err=sqrt(Gold*dx*dy);
    if info == 1
       disp(num2str([it err errP],'%10d   %10.16e   %10.16e'));
    end
end

end

function A=conn(Nx,dx,Ny,dy,E1,nu1,E2,nu2)

%function computing the influence coefficients
%between pressure and displacement
% A:matrix of the coefficients between P and Uz

Ee = 0.5*((1-nu1^2)/E1 + (1-nu2^2)/E2);
Ee = 1./Ee;
C = 2./pi/Ee;
px2 = dx/2; py2=dy/2;
XL=dx*Nx; YL=dy*Ny;
xb=0:dx:XL-px2; yb=0:dy:YL-py2;
[xxb,yyb]=ndgrid(xb,yb);

xxm=xxb-px2; xxp=xxb+px2;yym=yyb-py2;yyp=yyb+py2;
%Influence coefficient matrix (Eq. 3.13)
A=FNF(xxm,yym)+FNF(xxp,yyp) - FNF(xxm,yyp)-FNF(xxp,yym);A=C*A;
%Extension of the domain for FFT (Fig. 3.4)
A(Nx+1,:)=A(Nx,:);
A(Nx+2:2*Nx,:)=A(Nx:-1:2,:);
A(:,Ny+1)=A(:,Ny);
A(:,Ny+2:2*Ny)=A(:,Ny:-1:2);

end

function reco=FNF(x,y)
%kernel (Green's function) based on Love's theory (Eq. 3.16)
r=sqrt(x.*x+y.*y); reco=y.*log(x+r)+x.*log(y+r);

end
