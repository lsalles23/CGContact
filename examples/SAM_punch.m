% Main program
clear all
addpath ../src/
%##########################################################################
%% Creation of the grid
Nx=256;                 % Number of points in the x direction
Ny=256;                 % Number of points in the y direction
coef1=1.5*1e-3;
coef=1.5*1e-3;
x=linspace(-coef1,coef1,Nx);    
y=linspace(-coef,coef,Ny);    
[Mx,My]=ndgrid(x,y);
%##########################################################################
%% Separation between the 2 contacting surfaces
% rounded punch
Ra=1000e-3;                   
Rb=1000e-3;                   
c=1e-3;
h=zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        if My(i,j)>c
            %h(i,j)=(sqrt(Mx(i,j)^2+My(i,j)^2)-c)^2/2/R;
            h(i,j)=(My(i,j)-c)^2/2/Ra;
        elseif My(i,j)<-c
            h(i,j)=(My(i,j)+c)^2/2/Rb;
        else
            h(i,j)=0.0;
        end
    end
end
profile.h =h;
profile.x=x; profile.y=y;
%##########################################################################
%% Applied rigid body displacements :
loading.deltaz=0;            % Translation in the z direction
loading.phix=0;               % Rotation around x
loading.phiy=0;               % Rotation around y
%##########################################################################
%% Applied loads 
loading.N=1.00;              % Normal load
loading.Cx=0.00;             % Moment around x
loading.Cy=0.0e-2;             % Moment around y
%##########################################################################
%% Material properties
material.E1=2.1e11;              % Young modulus
material.nu1=0.3;                % Poisson's ratio
material.E2=2.1e11;              
material.nu2=0.3;
%##########################################################################
%% Control data
options.niter = 100;
options.errlim_p=1e-15;
options.info = 1;
%% Normal contact solution 
tic;
[P,rk,phix,phiy,deltaz,u]=CG_normal(profile,material,loading,options);
toc;
% Display
plot(x,P(Nx/2,:));
xlabel('r');
ylabel('Pressure p')
