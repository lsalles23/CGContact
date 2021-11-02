function [Qx,Qy,rkx,rky,delta0x,delta0y,delta0z,st,sl,s] = CG_tangential(x,y,P,Qoldx,Qoldy,Tx,Ty,E1,nu1,E2,nu2,mu,torque,coupling)

Nx=length(x); dx=x(2)-x(1);
Ny=length(y); dy=y(2)-y(1);

lbdij=zeros(Ny,Nx); % lambda = norm of the slips

[Ax,Ay,Axy]=conn(Ny,dy,Nx,dx,E1,nu1,E2,nu2);
fAx=fft2(Ax);
fAy=fft2(Ay);
fAxy=fft2(Axy);

clear Ax;
clear Ay;
clear Axy;

unity=ones(Nx,Ny);
[X,Y]=meshgrid(x,y);
%#################################################################################
% Area of contact : The entire contact area is supposed to be initially in stick
s=P>0;
st = s;             % stick region
sl = false(Ny,Nx);  % slip region
%#################################################################################
% Initialization of the shears that must verify the equilibrium equations

suminix=0.0;
suminiy=0.0;
suminixy=0.0;

for i=1:Nx
    for j=1:Ny
        if (st(i,j)==1)
            suminix=suminix+x(i);
            suminiy=suminiy+y(j);
            suminixy=suminixy+x(i)^2+y(j)^2;
        end
    end
end

M=zeros(3,3);
M(1,1)=sum(unity(st));
M(1,3)=-suminiy;
M(2,2)=sum(unity(st));
M(2,3)=suminix;
M(3,1)=-suminiy;
M(3,2)=suminix;
M(3,3)=suminixy;

TM=zeros(3,1);
TM(1,1)=Tx/dx/dy;
TM(2,1)=Ty/dx/dy;
TM(3,1)=torque/dx/dy;

Qx=zeros(Ny,Nx);
Qy=zeros(Ny,Nx);


AAA=linsolve(M,TM);

Qx=AAA(1,1)*st-AAA(3,1)*Y'.*st+Qoldx.*st;

Qy=AAA(2,1)*st+AAA(3,1)*X'.*st+Qoldy.*st;
%#################################################################################
errlim=1e-15; Gold=1; khi=0;

pkx=zeros(Ny,Nx); 
pky=zeros(Ny,Nx); 

err=1;it=0;

dQerrx=1000;
dQerry=1000;

QQx=zeros(2*Ny,2*Nx);
QQy=zeros(2*Ny,2*Nx);

PKx=zeros(2*Ny,2*Nx);
PKy=zeros(2*Ny,2*Nx);

while (((err>errlim)||abs(dQerrx>1)||abs(dQerry>1))&&(it<1000))
    it = it+1;
    %#############################################################################
    % Computation of the elastic deformations using the DC-FFT
    QQx(1:Ny,1:Nx)=(Qx-Qoldx);
    QQy(1:Ny,1:Nx)=(Qy-Qoldy);
    
    ddx=real(ifft2(fAx.*fft2(QQx,2*Ny,2*Nx)));
    ddy=real(ifft2(fAy.*fft2(QQy,2*Ny,2*Nx)));
    
    ddxduetoy=coupling*real(ifft2(fAxy.*fft2(QQy,2*Ny,2*Nx))); % contribution of Qy on ux
    ddyduetox=coupling*real(ifft2(fAxy.*fft2(QQx,2*Ny,2*Nx))); % contribution of Qx on uy
    
    ux=ddx(1:Ny,1:Nx)+ddxduetoy(1:Ny,1:Nx);
    uy=ddy(1:Ny,1:Nx)+ddyduetox(1:Ny,1:Nx);
    
    clear ddx;
    clear ddy;
    clear ddxduetoy;
    clear ddyduetox;
    %#############################################################################
    % Approximation of deltax, deltay and deltaz (rotation, needed when torsion
    % is applied)
    
    A=zeros(3,3);
    
    sumix=0.0;
    sumiy=0.0;
    sumixy=0.0;
    sumux=0.0;
    sumuy=0.0;
    sumt=0.0;
    
    for i=1:Nx
        for j=1:Ny
            if (st(i,j)==1)
                sumix=sumix+x(i);
                sumiy=sumiy+y(j);
                sumixy=sumixy+x(i)^2+y(j)^2;
                sumux=sumux+ux(i,j);
                sumuy=sumuy+uy(i,j);
                sumt=sumt-y(j)*ux(i,j)+x(i)*uy(i,j);
            end
        end
    end
    
    A(1,1)=sum(unity(st));
    A(3,1)=-sumiy;
    A(2,2)=sum(unity(st));
    A(3,2)=sumix;
    A(1,3)=sumiy;
    A(2,3)=-sumix;
    A(3,3)=-sumixy;
    
    B=zeros(3,1);
    
    B(1,1)=sumux;
    B(2,1)=sumuy;
    B(3,1)=sumt;

    %if (torque==0.0)&&(Ty==0.0)
       % keyboard
    %    delta0x=sumux/A(1,1);
    %    delta0y=0.0;
    %    delta0z=0.0;
    %end
    
    %if torque==0.0
    %    DELTA=linsolve(A(1:2,1:2),B(1:2,1));
    %    delta0x=DELTA(1,1); % rigid displacement in the x direction
    %    delta0y=DELTA(2,1); % rigid displacement in the y direction
    %    delta0z=0.0;        % rigid rotation around the z direction
    %end
    %    else
    DELTA=linsolve(A,B);
    delta0x=DELTA(1,1);
    delta0y=DELTA(2,1);
    delta0z=DELTA(3,1);
    %    end
    %#############################################################################
    % Calculations of the residuals rkx and rky
    
    rkx=ux-delta0x-delta0z*Y';
    rky=uy-delta0y+delta0z*X';
    
    clear ux;
    clear uy;
    
    %#############################################################################
    % Calculation of lbdij
    lbdij = lbdij.*not(sl)-sqrt((rkx.*sl).^2+(rky.*sl).^2).*sign((rkx.*sl).*(Qx.*sl)+(rky.*sl).*(Qy.*sl));
    %#############################################################################
    % Updtate of the stick region : if lambda(i,j) < 0 and (i,j) € sl, (i,j) is
    % added to st, and removed from sl
    stnew = (lbdij.*sl)<0; % logical
    
    if stnew~=0
        khi=0;
        st=st|stnew;
        sl=sl.*not(stnew);
        sl=logical(sl);
    end
    %#############################################################################
    % G norm2 of rk in the contact
    G=sum(sum((rkx.*st).*(rkx.*st)+(rky.*st).*(rky.*st)));
    %#############################################################################
    % Computation of the descent directions pk : in the slip zone, the descent directions are equal
    % to zero : the conjugate gradient is used in the stick region only
    pkx=pkx.*not(st)+rkx.*st+khi*G/Gold*(pkx.*st);%stick
    pky=pky.*not(st)+rky.*st+khi*G/Gold*(pky.*st);%stick
    pkx=pkx.*not(sl);%slip
    pky=pky.*not(sl);%slip
    
    Gold=G;
    %#############################################################################
    % Computation of qk (convolution of the influence coefficients matrix 
    % with the descent directions)
    PKx(1:Ny,1:Nx)=pkx;
    PKy(1:Ny,1:Nx)=pky;
    
    ddpx=real(ifft2(fAx.*fft2(PKx,2*Ny,2*Nx)));
    ddpxy=coupling*real(ifft2(fAxy.*fft2(PKy,2*Ny,2*Nx)));
    ddpy=real(ifft2(fAy.*fft2(PKy,2*Ny,2*Nx)));
    ddpyx=coupling*real(ifft2(fAxy.*fft2(PKx,2*Ny,2*Nx)));
    
    qkx=ddpx(1:Ny,1:Nx)+ddpxy(1:Ny,1:Nx);
    qky=ddpy(1:Ny,1:Nx)+ddpyx(1:Ny,1:Nx);
    
    clear ddpx;
    clear ddpxy;
    clear ddpy;
    clear ddpyx;
    %#############################################################################
    % Computation of the coeff's alpha=dp (length of the step to be made in
    % the computed descend direction)
    
    den1=sum(sum((qkx.*st).*(pkx.*st)));
    den2=sum(sum((qky.*st).*(pky.*st)));
    num1=sum(sum((rkx.*st).*(pkx.*st)));
    num2=sum(sum((rky.*st).*(pky.*st)));
    
    alpha=(num1+num2)/(den1+den2);
    
    Qx=Qx.*not(st)+Qx.*st-alpha*(pkx.*st);
    Qy=Qy.*not(st)+Qy.*st-alpha*(pky.*st);
    
    % Projection
    slnew=sqrt(Qx.^2.*st+Qy.^2.*st)>mu*P.*st; % Coulomb friction
    sl=sl|slnew;
    st=st.*not(slnew);
    st=logical(st);
    
    if (sqrt(sum(sum((Qx(sl)).^2+(Qy(sl)).^2))))~=0
        Qx(sl)=mu*P(sl).*Qx(sl)./sqrt(Qx(sl).^2+Qy(sl).^2);
        Qy(sl)=mu*P(sl).*Qy(sl)./sqrt(Qx(sl).^2+Qy(sl).^2);
    end
    %#############################################################################
    counterslnew=0;
    
    for i=1:Nx
        for j=1:Ny
            if(slnew(i,j)==1)
                counterslnew=counterslnew+1;
            end
        end
    end
    
    if (counterslnew>0)
        khi=0;
    else
        khi=1;
    end
    %#############################################################################
    % Calculation of the errors : 
    dQerrx = Tx-sum(sum(Qx(s)))*dx*dy;
    dQerry = Ty-sum(sum(Qy(s)))*dx*dy;
    err=sqrt(Gold*dx*dy);
    er(it)=err;
    %#############################################################################
    % Enforcement of the equilibrium equations : 
    
    sumx=0;
    sumy=0;
    sumxy=0;
    
    for i=1:Nx
        for j=1:Ny
            if (s(i,j)==1)
                sumx=sumx+x(i);
                sumy=sumy+y(j);
                sumxy=sumxy+x(i)^2+y(j)^2;
            end
        end
    end
    
    C=zeros(3,3);
    C(1,1)=sum(unity(s));
    C(2,2)=sum(unity(s));
    C(3,1)=-sumy;
    C(3,2)=sumx;
    C(1,3)=-sumy;
    C(2,3)=sumx;
    C(3,3)=sumxy;
    
    D=zeros(3,1);
    
    sumqx=0.0;
    sumqy=0.0;
    sumqxy=0.0;
    
    for i=1:Nx
        for j=1:Ny
            if (s(i,j)==1)
                sumqx=sumqx+Qx(i,j);
                sumqy=sumqy+Qy(i,j);
                sumqxy=sumqxy-y(j)*Qx(i,j)+x(i)*Qy(i,j);
            end
        end
    end
    
    D(1,1)=Tx/dx/dy-sumqx;
    D(2,1)=Ty/dx/dy-sumqy;
    D(3,1)=torque/dx/dy-sumqxy;
    
    ABC=linsolve(C,D);
    
    a=ABC(1,1);
    b=ABC(2,1);
    c=ABC(3,1);

    %keyboard
    
    Qx=Qx+a*st-c*Y'.*st;
    Qy=Qy+b*st+c*X'.*st;
    
    %#############################################################################
    disp(num2str([it err dQerrx dQerry],'%10d %10.2g %10.2g'));
    %#############################################################################
end

end

function [Ax,Ay,Axy]=conn(Nx,dx,Ny,dy,E1,nu1,E2,nu2)

%subroutine computing the influence coefficients d
%between pressure and displacement 
% A:matrix of the coefficients between P and u

C1 = -((1.0-nu1^2)/pi/ E1 + (1.0-nu2^2)/pi/ E2); 

C2 = -((1+nu1)/pi/E1+(1+nu2)/pi/E2);

C3 = -nu1*(1+nu1)/pi/E1;

C4 = -nu2*(1+nu2)/pi/E2;

C5 = (1+nu1)/pi/E1;

C6 = (1+nu2)/pi/E2;

px2 = dx/2; py2=dy/2;
XL=dx*Nx; YL=dy*Ny;
xb=0:dx:XL-px2; yb=0:dy:YL-py2;
[xxb,yyb]=ndgrid(xb,yb);

xxm=xxb-px2; xxp=xxb+px2;yym=yyb-py2;yyp=yyb+py2;

%Ax = FNF1(xxm,yym,C1,C2) + FNF1(xxp,yyp,C1,C2) - FNF1(xxm,yyp,C1,C2) - FNF1(xxp,yym,C1,C2);
Ax = FNF1(xxm,yym,C5,C6,nu1,nu2) + FNF1(xxp,yyp,C5,C6,nu1,nu2) - FNF1(xxm,yyp,C5,C6,nu1,nu2) - FNF1(xxp,yym,C5,C6,nu1,nu2);
Ay = FNF2(xxm,yym,C5,C6,nu1,nu2) + FNF2(xxp,yyp,C5,C6,nu1,nu2) - FNF2(xxm,yyp,C5,C6,nu1,nu2) - FNF2(xxp,yym,C5,C6,nu1,nu2);
Axy = FNF3(xxm,yym,C3,C4) + FNF3(xxp,yyp,C3,C4) - FNF3(xxm,yyp,C3,C4) - FNF3(xxp,yym,C3,C4);

Ax(Nx+1,:)=Ax(Nx,:);
Ax(Nx+2:2*Nx,:)=Ax(Nx:-1:2,:);
Ax(:,Ny+1)=Ax(:,Ny);
Ax(:,Ny+2:2*Ny)=Ax(:,Ny:-1:2);

Ay(Nx+1,:)=Ay(Nx,:);
Ay(Nx+2:2*Nx,:)=Ay(Nx:-1:2,:);
Ay(:,Ny+1)=Ay(:,Ny);
Ay(:,Ny+2:2*Ny)=Ay(:,Ny:-1:2);

Axy(Nx+1,:)=Axy(Nx,:);
Axy(Nx+2:2*Nx,:)=Axy(Nx:-1:2,:);
Axy(:,Ny+1)=Axy(:,Ny);
Axy(:,Ny+2:2*Ny)=Axy(:,Ny:-1:2);

end

function recox = FNF1(x,y,C1,C2,nu1,nu2) 
r=sqrt(x.*x+y.*y); 
recox=-C1*(y.*log(r+x)-y+(1-nu1)*x.*log(y+r))-C2*(y.*log(r+x)-y+(1-nu2)*x.*log(y+r));
end

function recoy = FNF2(x,y,C1,C2,nu1,nu2) 
r=sqrt(x.*x+y.*y); 
recoy=-C1*(x.*log(r+y)-x+(1-nu1)*y.*log(x+r))-C2*(x.*log(r+y)-x+(1-nu2)*y.*log(x+r));
end

function recoxy = FNF3(x,y,C1,C2) 
r=sqrt(x.*x+y.*y); 
recoxy=(C1+C2)*r;
end