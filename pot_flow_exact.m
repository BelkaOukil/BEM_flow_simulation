
clear all
clc

U_inf = 100;
a = 1 ;  
c =-a*2;
b =a*2;
n =a*50; 
% [x,y]=meshgrid([c:(b-c)/n:b],[c:(b-c)/n:b]');
X=linspace(-10,10,100);
Y=linspace(-10,10,100);
[x,y] = meshgrid(X,Y);

for i=1:length(x)
   for k=1:length(x)
      if sqrt(x(k,i).^2+y(k,i).^2)<a
         x(k,i)=1/0;
         y(k,i)=1/0;
      end
   end
end

r=sqrt(x.^2+y.^2);
theta=atan2(x,y);

psy = U_inf.*sin(theta).*r.*(1-(a^2./(r.^2)));
% phi = U_inf.*cos(theta).*r.*(1+(a^2./(r.^2)));

% n=100;
% r1=ones(1,n+1)*a;
% t=[0:2*pi/n:2*pi];

figure(1)
C_psy = contour(x,y,psy,25);
% clabel(C_psy,0)
hold on
% C_phi = contour(x,y,phi,50);
% clabel(C_phi,0)
viscircles([0 0],a,'color','k');
% polar(t,r1,'-k')
axis square
% title('potentiel de vitesse')
grid on