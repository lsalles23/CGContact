% Main program
clear all
addpath ../src
%##########################################################################

%% Creation of the grid
Nx=64;                 % Number of points in the x direction
Ny=64;                 % Number of points in the y direction
x=linspace(-0.1e-3,0.1e-3,Nx);    
y=linspace(-5e-3,5e-3,Ny);    
[hx,hy]=ndgrid(x,y);
%##########################################################################
%% Separation between the 2 contacting surfaces
%rough=zeros(Nx,Ny);
R=6e-3;                   % Equivalent radius : 1/R = 1/R1 + 1/R2 
h=zeros(Nx,Ny);
profile.x = x; profile.y = y;
profile.h=(hx.^2)/2/R;
%##########################################################################
%% Applied rigid body displacements :
loading.deltax=0.0;             % Translation in the x direction
loading.deltay=0.0;             % Translation in the y direction
loading.deltaz=0.0;            % Translation in the z direction
loading.phix=0.0;               % Rotation around x
loading.phiy=0.0;               % Rotation around y
loading.phiz=0.0;               % Torsional angle
%##########################################################################
%% Applied loads 
loading.N = 185;            % Normal load
loading.Cx = 0.4;             % Moment around x
loading.Cy = 0.0;             % Moment around y
%##########################################################################
%% Material properties
material.E1=2.0e11;              % Young modulus
material.nu1=0.3;                % Poisson's ratio
material.E2=2.0e11;              
material.nu2=0.3;
material.mu=0.5;                 % Friction coefficient
%##########################################################################
%% Control data
options.niter = 100;
options.errlim_p = 1e-8;
options.info = 0;
%% Normal contact solution 
[P,er,phix,phiy,deltaz,u]=CG_normal(profile,material,loading,options);
