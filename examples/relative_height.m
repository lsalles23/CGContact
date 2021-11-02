% clear variables; close all;
% Test=5;
% Pre=1; % 0=WORN SURFACES.  1=UNWORN SURFACES

function [h, x, y]=relative_height(x,y,A,B,Nx,Ny,Pre)

%%%%%%%%%%%%%%%%%
%Reduce size to Nx*Ny
xx=1:size(A,1);
yy=1:size(A,2);
x2 = interp1(0:size(xx,2)-1,xx,linspace(0,size(xx,2)-1,Nx));
y2 = interp1(0:size(yy,2)-1,yy,linspace(0,size(yy,2)-1,Ny));
[X2,Y2] = meshgrid(x2,y2);
%A2 = griddata(X(:),Y(:),A(:),X2,Y2,'nearest');
A2 = interpn(A,X2,Y2);
A2=imrotate(A2,90);
A2=flip(A2,1);
B2 = interpn(B,X2,Y2);
B2=imrotate(B2,90);
B2=flip(B2,1);
x2 = interp1(0:size(x,2)-1,x,linspace(0,size(x,2)-1,Nx));
y2 = interp1(0:size(y,2)-1,y,linspace(0,size(y,2)-1,Ny));
figure
subplot(1,2,1)
contourf(y,x,A)
title('original')
subplot(1,2,2)
contourf(y2,x2,A2)
title('Nx*Ny')
figure
subplot(1,2,1)
contourf(y,x,B)
title('original')
subplot(1,2,2)
contourf(y2,x2,B2)
title('Nx*Ny')
% PIXmm=PIXmm/10;%Pixel/mm

A=A2;
B=B2;
y=y2;
x=x2;

if Pre==1 %otherwise do not do it, because I performe cross-correlation of worn interfaces beforehand
    A=flip(A,1);
    %A=flip(A,2);
    %B=flip(B,1); %flip vertically
    %B=flip(B,2);%flip horizontally

    %Ruota di 90 grade i preworn (perche' li ho scannerizzati tutti nella
    %stessa direzione 
    if Pre==1
        B=imrotate(B,90);
    end
    B=-B;
    
    %Offset to put them on zero
    A=A-abs(max(A(:)));
    B=B+abs(min(B(:)));
end

%Centre x and y axes on zero (otherwise the code does not work)
x=x-mean(x(:));
y=y-mean(y(:));


figure;
subplot(1,2,1);
surf(y,x,A);
%axis equal;
xlabel('mm');ylabel('mm');zlabel('\mum');
title('Specimen A');
subplot(1,2,2);
contourf(y,x,A);
axis equal;
xlabel('mm');ylabel('mm');
title('Specimen A');
set(gcf,'Position',[100 200 1800 700])

figure;
subplot(1,2,1);
%mesh(Y,X,A)
surf(y,x,B)
xlabel('mm');ylabel('mm');zlabel('\mum');
title('Specimen B');
subplot(1,2,2);
contourf(y,x,B);
axis equal;
xlabel('mm');ylabel('mm');
title('Specimen B');
set(gcf,'Position',[100 200 1800 700])

x=x/1000;%convert in meters
y=y/1000;
h=A-B;
h=-h/1e6;%convert in meters
h=h-min(h(:));

figure;
subplot(1,2,1);
surf(y*1e3,x*1e3,h*1e6);
xlabel('mm');ylabel('mm');zlabel('Contact interface [\mum]');
title('Relative Height');
subplot(1,2,2);
contourf(y*1e3,x*1e3,h*1e6);
axis equal;
xlabel('mm');ylabel('mm');
title('Relative Height');
set(gcf,'Position',[100 200 1800 700])

%MeanB=take_mean(B,'B');

%lim=caxis;
%figure(a);
%caxis(lim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B_old1=B;
% A=A-mean(A(:));
% B=B-MeanB;
% B_old2=B;
% D=sum(A(:).*A(:));% For normalization %stdA*stdB*size(A,1)*size(A,2);
% 
% 
% %Rotate Surfaces
% deg=linspace(-20,20,1000); %degrees of rotation
% peak=zeros(size(deg));
% ij=zeros(size(deg));
% ji=zeros(size(deg));
% 
% for i=1:length(deg)
%     
%     B = imrotate(B_old2,deg(i));
%     %SurfRot=rotate(SurfA,direction,deg(i)) 
%     %figure;
%     %SurfA=surf(YA,XA,A);
% 
% crr = xcorr2(B,A);
% [peak(i),ind] = max(crr(:));
% %[peaks,inds] = findpeaks(crr(:),1,'MinPeakDistance',length(crr(:))/20);
% %[peaks,new]=sort(peaks,'descend');
% %inds=inds(new);
% [ij(i),ji(i)] = ind2sub(size(crr),ind);
% 
% end
% 
% figure;
% plot(deg,peak/D);
% xlabel('Degree of Rotation');
% ylabel('Cross-correlation');
% hold on;
% p = polyfit(deg,peak/D,4);
% Y=p(1)*deg.^4+p(2)*deg.^3+p(3)*deg.^2+p(4)*deg.^1+p(5);
% [m,index]=max(Y);
% plot(deg,Y);
% plot(deg(index),Y(index),'or');
% plot(deg(index),peak(index)/D,'or');
% 
% %Plot for the best Angle
% B = imrotate(B_old1,deg(index));
% Bmean=B(ij(index)+1-size(A,1):ij(index),ji(index)+1-size(A,2):ji(index));
% B=B-mean(Bmean(:));
% 
% crr = xcorr2(B,A);
% [peaks,inds] = findpeaks(crr(:),1,'MinPeakDistance',length(crr(:))/20);
% [peaks,new]=sort(peaks,'descend');
% inds=inds(new);
% [ij1,ji1] = ind2sub(size(crr),inds(1));
% D=sum(A(:).*A(:));%stdA*stdB*size(A,1)*size(A,2);
% 
% 
% %Plot Cross-Correlation
% figure;
% plot(crr(:)/D)
% title('Cross-Correlation')
% hold on
% plot(inds,peaks/D,'or')
% title(['Rotation Angle=',num2str(deg(index))]);
% 
% 
% % figure;
% % X=1:size(B,1);
% % Y=1:size(B,2);
% % %mesh(Y,X,A)
% % contourf(Y,X,B)
% % title('Specimen B');
% % hold on;
% % plot(ji,ij,'or')
% 
% %PLOT A rotated
% figure;
% X=1:size(A,1);
% Y=1:size(A,2);
% contourf(Y,X,A);
% title('Specimen A');
% % figure;
% % boxplot(A(:));
% % title('Specimen A - Boxplot');
% 
% %PLOT THE PORTION OF B that overlaps with A
% figure;
% X=1:size(A,1);
% Y=1:size(A,2);
% %mesh(Y,X,A)
% Bred=B(ij1+1-size(A,1):ij1,ji1+1-size(A,2):ji1);
% contourf(Y,X,Bred);
% title('Specimen B - rotated and cross correlated');
% % figure;
% % boxplot(Bred(:));
% % title('Specimen B - Boxplot');
% 
% figure; boxplot([A(:),Bred(:)],'PlotStyle','compact')
% pbaspect([1 10 1])
% names = {'A'; 'B'};
% set(gca,'xtick',[1:2],'xticklabel',names)
% title('Boxplots of final surfaces');
% 
% %Plot the differences
% figure;
% X=1:size(A,1);
% Y=1:size(A,2);
% %mesh(Y,X,A)
% diff=A-Bred;
% ERR=sum(dot(diff,diff))/length(diff);
% Av_ERR=sqrt(ERR)/15; %average micrometer height difference (from pixel to micrometer you divide by roughly 15)
% contourf(Y,X,diff);
% title('Difference');