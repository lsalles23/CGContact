function [P,X2] = Equilibrium(P,s,loading,profile,material)
% Equilibrium update pressures distribution to respect static equilibrium 
%
% [P,X2] = Equilibrium(P,s,loading,profile,material)
%
%
% Input data description
% P                          : array of pressure
% s                          : boolean matrix of active area 
% loading
%         N, Cx, Cy          : normal load, moments around x and y axis (N - N.m)
%         deltaz, phix, phiy : rigid body motions : translation along z axis
%                              rotations around x and y axis
% profile                    : profile structure
%        .x,.y               : coordinates of all grid points (m)
%        .h                  : initial separation between the two surfaces (m)
%
% material
%         E1,E2 nu1,nu2      : Young's modulus and Poisson ratio of the material
%         H                  : hardness of the softer material - represents the
%                              maximum allowable pressure 
% Ouput data description
% P                          : updated pressure 
% X2                         : coefficient of the transformation for updating
%
% see also RigidBodyMotion

%Copyright 2019 Imperial College London 
%authors Loic Salles Jason Armand

if isfield(material,'H')
    H = material.H; %hardness
else
    H = Inf; %pure elastic contact
end
x = profile.x; y = profile.y;
Nx=length(x); dx=x(2)-x(1); Ny=length(y); dy=y(2)-y(1);
sumpx=0.0;
sumpy=0.0;
sumpx2=0.0;
sumpy2=0.0;
sumpxy=0.0;
sump=0;
    
for i=1:Nx
   for j=1:Ny
       if (s(i,j)>0)||(P(i,j)==H)
           sump=sump+P(i,j);
           sumpx=sumpx+P(i,j)*x(i);
           sumpy=sumpy+P(i,j)*y(j);
           sumpx2=sumpx2+P(i,j)*x(i)^2;
           sumpy2=sumpy2+P(i,j)*y(j)^2;
           sumpxy=sumpxy+P(i,j)*x(i)*y(j);
        end
    end
end
    
M2=zeros(3,3);

M2(1,1)=sump;
M2(1,2)=sumpx;
M2(1,3)=sumpy;
M2(2,1)=sumpx;
M2(2,2)=sumpx2;
M2(2,3)=sumpxy;
M2(3,1)=sumpy;
M2(3,2)=sumpxy;
M2(3,3)=sumpy2;

N2=zeros(3,1);

N2(1,1) = loading.N/dx/dy;
N2(2,1) = loading.Cy/dx/dy;
N2(3,1) =-loading.Cx/dx/dy;
%solve equilibrium equation (Eq. 3.62)
X2=linsolve(M2,N2);
% adjust the pressure (Eq. 3.60)
for i=1:Nx
   for j=1:Ny
       if (s(i,j)>0)
           P(i,j)=P(i,j)*(X2(1,1)+X2(2,1)*x(i)+X2(3,1)*y(j));
       end
   end
end
