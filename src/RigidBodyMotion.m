function X1 = RigidBodyMotion(uz,s,profile)
% RigidBodyMotion computes the rigid motion 
%
% X1 = RigidBodyMotion(uz,s,profile)
%
%
% Input data description
% uz                    : normal displacement at the interface 
% s                     : boolean matrix of active area 
% profile               : profile structure
%        .x,.y          : coordinates of all grid points (m)
%        .h             : initial separation between the two surfaces (m)
%
% Ouput data description
% X1                    : rigid body motion [deltaz,phix,phiz]
%
% see also Equilibrium

%Copyright 2019 Imperial College London 
%authors Loic Salles Jason Armand

x = profile.x; y = profile.y;
Nx=length(x); dx=x(2)-x(1); Ny=length(y); dy=y(2)-y(1);
h = profile.h;

    sum1=0.0;
    sumx=0.0;
    sumy=0.0;
    sumxy=0.0;
    sumx2=0.0;
    sumy2=0.0;
    sumN1=0.0;
    sumN2=0.0;
    sumN3=0.0;
    
    for i=1:Nx
        for j=1:Ny
            if s(i,j)>0
                sum1=sum1+1.0;
                sumx=sumx+x(i);
                sumy=sumy+y(j);
                sumxy=sumxy+x(i)*y(j);
                sumx2=sumx2+x(i)^2;
                sumy2=sumy2+y(j)^2;
                sumN1=sumN1+uz(i,j)+h(i,j);
                sumN2=sumN2+(uz(i,j)+h(i,j))*y(j);
                sumN3=sumN3+(uz(i,j)+h(i,j))*x(i);
            end
        end
    end
    
    M1=zeros(3,3);
    
    M1(1,1)=sum1;
    M1(1,2)=sumy;
    M1(1,3)=-sumx;
    M1(2,1)=sumy;
    M1(2,2)=sumy2;
    M1(2,3)=-sumxy;
    M1(3,1)=sumx;
    M1(3,2)=sumxy;
    M1(3,3)=-sumx2;
    
    N1=zeros(3,1);
    
    N1(1,1)=sumN1;
    N1(2,1)=sumN2;
    N1(3,1)=sumN3;

    X1=linsolve(M1,N1); %(Eq. 3.47)
