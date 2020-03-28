clear
close all
clc
more off

% PROBLEMA 1
% numero de Reynolds para esfera lisa com
% diametro D
% velocidade (modulo) V
% viscosidade visco
Re=@(rho,D,V,visco) rho*D*V/visco;

% velocidade em funcao de Re, rho, D e visco
Velocidade=@(Re,rho,D,visco) Re*visco/(rho*D);

% O Coeficiente de Arrasto (CA) para bola lisa, para cada numero de Re eh medido experimentalmente no tunel do vento.
% Contudo, Morrison escreveu uma funcao que calcula o CA aproximadamente para cada Re (1e-1<=Re<=1e6).
% Referencia: F A. Morrison, Data Correlation for Drag Coefficient for Sphere, http://www.chem.mtu.edu/~fmorriso/ DataCorrelationForSphereDrag2013.pdf
CAMorrison=@(re) 24./re+2.6*(re/5)./(1+(re/5).^(1.52))+0.411*(re/263000).^(-7.94)./(1+(re/263000).^(-8))+re.^0.8/461000;
%nReynolds=1e-1:1e6;
%loglog(nReynolds,CAMorrison(nReynolds))

% bola de tenis, ar
rho=1.224; % kg/m3, densidade do ar
D=6.5e-2; % 6.5cm, diametro bola de tenis
r=D/2; % raio da bola
A=pi*r*r; % area da bola
massa=0.058; % kg, massa da bola
visco=1.83e-5; %kg/(ms), viscosidade do ar
%visco=1e-3;


CM=1; % coeficiente de Magnus
g=9.81*[0 0 -1]'; % m/s2, aceleracao da gravidade
v0=[26 0 5]'; % m/s, velocidade inicial
rey=Re(rho,D,norm(v0),visco);
if rey<1e-1 || rey>1e6
  fprintf('A formula de Morrison pode nao ser aplicada...\n')
end
p0 = [0 0 1.5]'; % m, posicao inicial
omega = 500*[1 1 1]'; % Hz, velocidade de rotacao
%omega=20*omega;

% forca de arrasto
FArrasto=@(CA,rho,A,nv,vuni) 0.5*CA*rho*A*nv*nv*(-vuni);

% forca de sustentacao
FSustentacao=@(CM,rho,A,r,omega,vel) 0.5*CM*rho*A*r*cross(omega,vel);

% metodo de Euler para calcular a velocidade e posicao a partir da aceleracao
tini=0;
tfin=3.165;
n=70;
h=(tfin-tini)/n;
vel=v0;
pos=zeros(3,n);
pos(:,1)=p0;
tempo=0;

for i=1:n-1
  nv=norm(vel);
  if nv<1e-6
    CA=0;
    vuni=vel*0;
  else
    CA=CAMorrison(Re(rho,D,nv,visco));
    vuni=vel/nv;
  end
  FA=FArrasto(CA,rho,A,nv,vuni);
  FS=FSustentacao(CM,rho,A,r,omega,vel);
  FG=massa*g;
  vel=vel+(h/massa)*(FA+FS+FG);
  pos(:,i+1)=pos(:,i)+h*vel;
  if pos(3,i+1)<0
    break
  end
  tempo=tempo+h;
end

fprintf('Tempo de voo: %f s\n',tempo)
plot3(pos(1,1:i),pos(2,1:i),pos(3,1:i),'or')
hold

% quadra
rectangle ("Position", [0, -5.48, 23.77, 10.97]);
rectangle ("Position", [0, -4.11, 23.77, 8.23]);
rectangle ("Position", [5.48, -4.11, 12.81, 8.23]);
line ( [5.48 18.30],[0 0]);
line ([11.88 11.88],[-5.48 5.48]);

% rede
line([11.88 11.88],[5.48 5.48],[0 0.9])
line([11.88 11.88],[-5.48 -5.48],[0 0.9])
line([11.88 11.88],[5.48 -5.48],[0.9 0.9])

grid
view([0 0])
axis equal

% drawCircle(21,0,raio);
% angle = linspace(0, 2*pi,360);
% xCirc = cos(angle);
% yCirc = sin(angle);
% plot(xCirc,yCirc)
% axis('equal') 

