function [PN]=Reynolds_Tex(H,b2,Omega)
% Solves for the normal force and the apparent viscosity calculated
% using a 2nd order fluid constitutive model. 
% The pressure and flow fields with a Newtonian fluid are solved
% using the Reynolds equation in cylindrical coordinates
% throught a pseudo-spectral Galerkin method.
%
% The pressure with the 2nd order fluid is obtained using
% the Three-Dimensional Flow Theorem of Giesekus, which
% assumes that b2=b11, where b2 and b11 are 2nd order fluid
% coefficient that relate to the 1st and 2nd normal stress differences.
% The velocity field is the same as that obtained with a Newtonian fluid.
%
% The input parameters for the function are:
% H (a matrix of values for the gap height at the given meshpoints)
%   (for mesh independece, H should be at least 25x25)
% 
% The output parameters are:
% tstar (the shear stress/shear stress with no texture)
% Fn (the normal force on the flat surface)
% mu (the effective friction coefficient)

%format long
%Define variables used to generate the mesh
R1=0.01e-3; R2=20e-3; h0=0.25e-3 + 19e-6; Rc=14.25e-3; Rt=3e-3;
phi=2*pi/10; Ntex=2*pi/phi;
A=size(H); N=A(1)-1;
%Parameter for 2nd order fluid: b2=0 is Newtonian

%Define flow and fluid variables
b1=1.4;

%Create the 1-D matrices for the psuedo-spectral method
Kh=zeros(N+1); Mh=zeros(N+1); Ch=zeros(N+1); Dh=zeros(N+1);
z=zeros(N+1,1); w=zeros(N+1,1);
[Kh,Mh,Ch,Dh,z,w]=semhat(N);
q=zeros(1,N); q(N)=1; Q=[q;eye(N)];
I=eye(N+1); Res=I(2:N,:); Pro=Res';
Nuem=I(1:N,:)';
r=(R2-R1)/2*(z)+(R2+R1)/2; theta=phi/2*(z);
[Rmat,Theta]=ndgrid(r,theta);
Rdiff=R2-R1;

R1d=diag(r); R=kron(I,R1d);
Rinv1d=diag(1./r); Rinv=kron(I,Rinv1d);

dhdz=H*Dh';
Hmat=H;

dh=zeros((N+1)^2,1);
h=zeros((N+1)^2,1);
for j=1:N+1,
    dh((j-1)*(N+1)+1:j*(N+1))=dhdz(:,j);
    h((j-1)*(N+1)+1:j*(N+1))=H(:,j);
end

eta_mat=b1;

%Create the 2-D matrices used for solving the Reynolds equation
H=diag(h.^3);
Kb=-phi/Rdiff*(kron(I,Dh))'*(kron(Mh,Mh)*R*H)*kron(I,Dh)...
   -Rdiff/phi*(kron(Dh,I))'*(kron(Mh,Mh)*H*Rinv)*kron(Dh,I);
fb=3*eta_mat*Rdiff*Omega*kron(Mh,Mh)*R*dh;

%K=(kron(Q,Nuem))'*Kb*kron(Q,Nuem); K=sparse(K);
%f=kron(Q,Nuem)'*fb; f=sparse(f);
K=(kron(Q,I))'*Kb*kron(Q,I);
f=kron(Q,I)'*fb;
%K=(kron(Pro,Pro))'*Kb*kron(Pro,Pro); K=sparse(K);
%f=kron(Pro,Pro)'*fb; f=sparse(f);

Con=condest(K);

%tstar=0; Fn=0; mu=0;

%K=full(K); K=distributed(K);
%f=full(f); f=distributed(f);

%Solve for the Newtonian pressure

%
%if (Con>=10e-15)
%p=K\f;
p = pinv(K)*f;
%else 
 %   p = pinv(K)*f;
%end

P=(kron(Q,I))*p;

%P=gather(p);

%poolobj = gcp('nocreate');
%delete(poolobj);


%Create the Newtonian Pressure matrix
%{
Pn=zeros(N+1,N+1);
for j=1:N+1,
    for i=1:N+1,
        Pn(i,j)=P((i)+(j-1)*(N));
    end
end
PN=Pn;
%}

%
Pn=zeros(N+1,N+1);
for j=1:N+1,
    for i=1:N+1,
        Pn(i,j)=P((i)+(j-1)*(N+1));
    end
end
%Pn=I*Pn*Pro';

PN=Pn;
Pr0=PN(N+1,:)';
Pref=-Pr0(floor(N/2+1));
PN=PN+Pref*ones(N+1);

Pr0=PN(N+1,:)';
Pref=-1/2*(Pr0'*w);
PN=PN+Pref*ones(N+1);

pn= reshape(PN,(N+1)^2,1);
dd = [kron(I,Dh);kron(Dh,I)];
dp = dd*pn;
%}

%{
Pn=zeros(N+1,N);
for j=1:N,
    for i=1:N+1,
        Pn(i,j)=p((i)+(j-1)*(N+1));
    end
end
Pn=I*Pn*Q';

PN=Pn;
Pr0=-Pn(N+1,floor(N/2+1));
PN=PN+Pr0;
%}

%Calculate the Non-Newtonian Pressure matrix
P=PN+b2./b1.*Omega.*(2/phi*PN*Dh')+b2./2*((1./(2*b1).*(2/Rdiff*Dh*PN).*Hmat)^2+...
    (1./(2*b1.*Rmat).*(2/phi*PN*Dh').*Hmat).^2+(Omega*Rmat./Hmat)^2+...
    Omega/b1*(2/phi*PN*Dh'));
%Calculate tauzz
tauzz=b2*((1./(2*b1).*(2/Rdiff*Dh*PN).*Hmat)^2+...
    (1./(2*b1.*Rmat).*(2/phi*PN*Dh').*Hmat).^2+(Omega*Rmat./Hmat)^2+...
    Omega/b1*(2/phi*PN*Dh'));
%Compute Normal Force
Fn=Ntex*(Rdiff*phi/4*w'*(R1d'*(P-tauzz))*w);
%Compute torque and viscosity
tau=-(-1./(2*Rmat).*(2/phi*PN*Dh').*Hmat-b1.*Rmat.*Omega./(Hmat)-...
    b2.*Omega./(2*b1.*Rmat).*(2/phi*(2/phi*PN*Dh')*Dh').*Hmat+...
    b2.*Rmat.*Omega^2./(Hmat.^2).*(2/phi*Hmat*Dh'));
M=Ntex*(Rdiff*phi/4*w'*(R1d'.^2*tau)*w);
%eta_a=2*h0/(pi*R2^4)*M/Omega;
tstar=(2/(pi*R2^3)*M)/(b1*(Omega*R2)/h0);
%Compute effective friction coefficient
mu=(M/R2)/Fn;
%}

%Calculate velocity
%{
dPdr=2/Rdiff*(Dh*P); dPdth=2/phi*(P*Dh');
i=1; j=1;
while r(i)<Rc,
    i=i+1;
end
rcent=i;

while theta(j)<0,
    j=j+1;
end
theta0=j;

dprRc=dPdr(i,:); dprTc=dPdr(:,j);
dptRc=dPdth(i,:); dptTc=dPdth(:,j);
prc=p(i,:);

DPDRRC=ones(N+1,1)*dprRc; DPDRTC=dprTc*ones(1,N+1);
DPDTRC=ones(N+1,1)*dptRc; DPDTTC=dptTc*ones(1,N+1);
Prc=ones(N+1,1)*prc;

zRc=zeros(N+1); zTc=zeros(N+1);
for k=1:N+1,
    zRc(:,k)=Hmat(i,k)/2*(1+z);
    zTc(k,:)=Hmat(k,j)/2*(1+z);
end

hRc=Hmat(i,:); hthetac=Hmat(:,j); Prc=P(i,:); Ptc=P(:,j);
HRC=zeros(N+1); HTC=zeros(N+1); PRC=zeros(N+1); PTC=zeros(N+1);
for k=1:N+1,
    HRC(:,k)=hRc(k)*ones(N+1,1);
    PRC(:,k)=Prc(k)*ones(N+1,1);
    HTC(k,:)=hthetac(k)*ones(1,N+1);
    PTC(k,:)=Ptc(k)*ones(1,N+1);
end

%
vtRc=1/(2*b1*r(i)).*(DPDTRC).*(zRc.^2-zRc.*HRC)+r(i)*Omega*(HRC-zRc)./HRC;
%{
vzrc=-(1/r(i)*(2/Rdiff)*Dh*(Rmat.*DPDRRC)+1/r(i)*(2/phi)*(1./Rmat.*DPDTRC)*Dh').*(zRc.^3-HRC.^3)./(6*b1)...
     +(1/r(i)*(2/Rdiff)*Dh*(Rmat.*HRC.*DPDRRC)+1/r(i)*(2/phi)*(1./Rmat.*HRC.*DPDTRC)*Dh').*(zRc.^2-HRC.^3)./(4*b1);%...
     +Omega.*(zRc.^2-HRC.^2)./2.*(2/phi)*(1./HRC*Dh');
%}
vzrc=zeros(N+1);
vrtc=1/(2*b1).*(DPDRTC).*(zTc.^2-zTc.*HTC);
%{
vztc=-1./(2*b1*Rmat).*(zTc.^3/3-zTc.^2.*HTC/2).*(Dh*PTC+(PTC*Dh')*Dh')+...
    1./(2*b1*Rmat).*(PTC*Dh').*zTc.^2/2.*(HTC*Dh')+...
    1./(2*b1*Rmat).*(Dh*PTC).*zTc.^2/2.*(Dh*HTC)-...
    1./(2*b1).*(Dh*(Dh*PTC)).*(zTc.^3/3-zTc.^2.*HTC/2)-...
    Omega.*zTc.^2./(2*HTC).*(HTC*Dh');
%}
vztc=zeros(N+1);
%{
Color(:,:,1)=zeros(N+1);
Color(:,:,2)=zeros(N+1);
Color(:,:,3)=zeros(N+1);

map1=zeros(64,3); map2=zeros(64,3); map3=zeros(64,3);
map4=zeros(64,3); map5=zeros(64,3); map6=zeros(264,3);
for j=1:64,
    map1(j,3)=-0.8/63*(j-1)+1;
    map1(j,2)=0.5/63*(j-1);
    map1(j,1)=1/63*(j-1);
    
    map2(j,3)=1/64*(64-j);
    map2(j,2)=1/64*(64-j);
    map2(j,1)=1/64*(j);
    
    map3(j,3)=-0.064516129032258+0.065524193548387*j-0.001008064516129*j^2;
    map3(j,2)=0.5*(-0.064516129032258+0.065524193548387*j-0.001008064516129*j^2);
    if map3(j,3)>1,
        map3(j,3)=1;
    end
    map3(j,1)=(-0.256016385048644e-3+0.012000768049156e-3*j+0.244015616999488e-3*j^2);
    if map3(j,1)>1,
        map3(j,1)=1;
    end
    
    map4(j,3)=1.019559651817716-0.019666858678955*j+0.000107206861239*j^2;
    if map4(j,3)>1,
        map4(j,3)=1;
    end
    map4(j,2)=-0.024321556579621+0.024577572964670*j-0.000256016385049*j^2;
    map4(j,1)=-0.016385048643113+0.016393049155146*j-0.000008000512033*j^2;
    if map4(j,1)>1,
        map4(j,1)=1;
    end
    %
    if map4(j,3)<0,
        map4(j,3)=0;
    end
    if map4(j,2)<0,
        map4(j,2)=0;
    end
    if map4(j,1)<0,
        map4(j,1)=0;
    end
    %
    
    map5(j,3)=1.013108038914491-0.013114439324117*j+0.000006400409626*j^2;
    if map5(j,3)>1,
        map5(j,3)=1;
    end
    map5(j,2)=-0.030773169482847+0.031129992319508*j-0.000356822836662*j^2;
    map5(j,1)=-0.022836661546339+0.022945468509985*j-0.000108806963646*j^2;
    if map5(j,1)>1,
        map5(j,1)=1;
    end
    %
    if map5(j,3)<0,
        map5(j,3)=0;
    end
    if map5(j,2)<0,
        map5(j,2)=0;
    end
    if map5(j,1)<0,
        map5(j,1)=0;
    end
end
%}

vt0=Rmat*Omega; vr0=0*Rmat;
vx0=vr0.*cos(Theta)-vt0.*sin(Theta); vy0=vr0.*sin(Theta)+vt0.*cos(Theta);
magv0=sqrt(vt0.^2+vr0.^2);
vth0d2=1./(2*b1*Rmat).*(2/phi*P*Dh').*((h0/2)^2-(h0/2)*Hmat)+Rmat*Omega.*(Hmat-(h0/2))./Hmat;
vrh0d2=1./(2*b1).*(2/Rdiff*Dh*P).*((h0/2)^2-(h0/2)*Hmat);
vxh0d2=vrh0d2.*cos(Theta)-vth0d2.*sin(Theta);
vyh0d2=vrh0d2.*sin(Theta)+vth0d2.*cos(Theta);
magvh0d2=sqrt(vth0d2.^2+vrh0d2.^2);
vth0=1./(2*b1*Rmat).*(2/phi*P*Dh').*((h0)^2-(h0)*Hmat)+Rmat*Omega.*(Hmat-(h0))./Hmat;
vrh0=1./(2*b1).*(2/Rdiff*Dh*P).*((h0)^2-(h0)*Hmat);
vxh0=vrh0.*cos(Theta)-vth0.*sin(Theta);
vyh0=vrh0.*sin(Theta)+vth0.*cos(Theta);
magvh0=sqrt(vxh0.^2+vyh0.^2);

for j=1:264,
    map6(j,3)=1.003065045133951-0.003065133088643*j+0.000000087954693*j^2;
    if map6(j,3)>1,
        map6(j,3)=1;
    end
    map6(j,2)=-0.007259164659101+0.007279460204459*j-0.000020295545359*j^2;
    map6(j,1)=-0.005358023974690+0.005363916939106*j-0.000005892964416*j^2;
    if map6(j,1)>1,
        map6(j,1)=1;
    end
    %
    if map6(j,3)<0,
        map6(j,3)=0;
    end
    if map6(j,2)<0,
        map6(j,2)=0;
    end
    if map6(j,1)<0,
        map6(j,1)=0;
    end
end
%}
%{
for i=1:N+1,
    for j=1:N+1,
        if P(i,j)>100,
            Color(i,j,1)=1;
        elseif P(i,j)<-100
            Color(i,j,3)=1;
        else
            Color(i,j,3)=1;
            Color(i,j,1)=1;
        end
    end
end
%}
%
%Figures
%{
f1=figure();
surf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,P)
colormap(map6)
view([80,30])
xlabel('--X-- [mm]')
ylabel('--Y-- [mm]')
zlabel('--P-- [Pa]')
fname=['\Users\Jonathon\Dropbox\Grad Research\2015-03-03-Reynolds_No_Cavitation\Cylindrical\Pressure_Profiles\Surf_P_b' num2str(beta*180/pi) '_degrees.jpg'];
saveas(f1,fname,'jpeg')


f2=figure();
contourf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,P)
colormap(map6)
colorbar
xlabel('--X-- [mm]')
ylabel('--Y-- [mm]')
zlabel('--P-- [Pa]')
fname=['\Users\Jonathon\Dropbox\Grad Research\2015-03-03-Reynolds_No_Cavitation\Cylindrical\Pressure_Profiles\Contour_P_b' num2str(beta*180/pi) '_degrees.jpg'];
saveas(f2,fname,'jpeg')


f3=figure();
surf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,-Hmat*1e3)
colormap(map6)
view([80,30])
xlabel('--X-- [mm]')
ylabel('--Y-- [mm]')
zlabel('-H(r,\theta) [mm]')
axis([0,20,-10,10,-10,10])
fname=['\Users\Jonathon\Dropbox\Grad Research\2015-03-03-Reynolds_No_Cavitation\Cylindrical\Texture_Profiles\H_b' num2str(beta*180/pi) '_degrees.jpg'];
saveas(f3,fname,'jpeg')


f3=figure();
mesh(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,-Hmat*1e3)
%colormap(map6)
hold on
mesh(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,zeros(N+1))
colormap([0,0,0])
text(15,6,0,'\rightarrow r\Omega','FontSize',18)
text(15,0,-(h1*1e3)-1,'Fixed','FontSize',18);
view([90,0])
xlabel('--X-- [mm]')
ylabel('--Y-- [mm]')
zlabel('-H(r,\theta) [mm]')
axis([0,20,-10,10,-10,10])
fname=['\Users\Jonathon\Dropbox\Grad Research\2015-03-03-Reynolds_No_Cavitation\Cylindrical\Texture_Profiles\H_b' num2str(beta*180/pi) '_degrees+moving_surf.jpg'];
saveas(f3,fname,'jpeg')
%}
%{
f1=figure('Position', [100, 50, 650, 1250]);
subplot(3,2,1)
surf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,-Hmat*1e3,'EdgeColor',[0.5,0.5,0.5])
hold on
plot3((2*r)*1e3,zeros(N+1,1),0*ones(N+1,1),'-k')
text(R1*1e3-0.5,0-0.2,0.25,'\uparrow \alpha','FontSize',14)
text(30-1.5,0-0.2,0.25,'\uparrow \alpha*','FontSize',14)
plot3(Rc*cos(theta)*1e3,Rc*sin(theta)*1e3,0*ones(N+1,1),'--k')
text(Rc*cos(-phi/2)*1e3,Rc*sin(-phi/2)*1e3-1.75,0.5,'\beta \uparrow','FontSize',14)
text(Rc*cos(phi/2)*1e3,Rc*sin(phi/2)*1e3-.75,0.5,'\uparrow \beta*','FontSize',14)
colormap(map6)
view([80,30])
text(0,-9,4,'(A)','Fontsize',14);
xlabel('X [mm]')
ylabel('Y [mm]')
zlabel('-H(r,\theta) [mm]')
axis([0,30,-10,10,-3,3])

%{
subplot(3,2,2)
mesh(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,-Hmat*1e3)
colormap([0,0,0])
hold on
surf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,zeros(N+1))
colormap([0,0,0])
text(15,6,0,'\rightarrow r\Omega','FontSize',14)
text(15,0,-(h1*1e3)-0.5,'Fixed','FontSize',14);
text(0,-9,3.5,'(B)','Fontsize',14);
view([90,0])
xlabel('X [mm]')
ylabel('Y [mm]')
zlabel('-H(r,\theta) [mm]')
axis([0,20,-10,10,-3,3])
%}
subplot(3,2,2)
%{
%surf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,-Hmat*1e3-1.5,'EdgeColor',[0.5,0.5,0.5])
colormap(map6)
view([80,30])
hold on
contourf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,magvh0d2);
%contourf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,magvh0);
%contourf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,magv0);
c=colorbar;
hold on
%quiver(Rmat(1:4:end,1:4:end).*cos(Theta(1:4:end,1:4:end))*1e3,Rmat(1:4:end,1:4:end).*sin(Theta(1:4:end,1:4:end))*1e3,vx0(1:4:end,1:4:end),vy0(1:4:end,1:4:end),1,'k')
quiver(Rmat(1:4:end,1:4:end).*cos(Theta(1:4:end,1:4:end))*1e3,Rmat(1:4:end,1:4:end).*sin(Theta(1:4:end,1:4:end))*1e3,vxh0d2(1:4:end,1:4:end),vyh0d2(1:4:end,1:4:end),1,'k')
%quiver(Rmat(1:4:end,1:4:end).*cos(Theta(1:4:end,1:4:end))*1e3,Rmat(1:4:end,1:4:end).*sin(Theta(1:4:end,1:4:end))*1e3,vxh0(1:4:end,1:4:end),vyh0(1:4:end,1:4:end),1,'k')
text(0,-9,4,'(B)','Fontsize',14);
xlabel('X [mm]')
ylabel('Y [mm]')
ylabel(c,'m/s')
zlabel('-H(r,\theta) [mm]')
axis([0,20,-10,10,-3,3])
%}
contourf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,magvh0d2)
hold on
quiver(Rmat(1:4:end,1:4:end).*cos(Theta(1:4:end,1:4:end))*1e3,Rmat(1:4:end,1:4:end).*sin(Theta(1:4:end,1:4:end))*1e3,vxh0d2(1:4:end,1:4:end),vyh0d2(1:4:end,1:4:end),1,'k')
colormap(map6)
caxis([0,20e-3*Omega])
axis([0,20,-10,10])
ax=gca;
ax.YAxisLocation='right';
text(1.5,-8.5,'(B)','Fontsize',14);
camroll(-90)
colorbar
xlabel('X [mm]')
ylabel('Y [mm]')
zlabel('P [Pa]')


subplot(3,2,3)
quiver(Rmat*1e3,-zTc*1e3,vrtc,vztc,1,'-k')
hold on
plot(r*1e3,zeros(1,length(theta)),'-r')
hold on
plot(r*1e3,-zTc(:,N+1)*1e3,'-r')
plot(zeros(1,10),linspace(-0.75,0.125,10),'--k')
plot(20*ones(1,10),linspace(-0.75,0.125,10),'--k')
text(1.5,0.55,'(C)','Fontsize',14);
text(-1.5,-0.75,'\alpha','Fontsize',14);
text(21,-0.75,'\alpha*','Fontsize',14);
xlabel('r [mm]')
ylabel('z [mm]')
axis([-2,25,-2.5,0.25])

subplot(3,2,4)
quiver(Rc.*(Theta(1:4:end,1:4:end))*1e3,-zRc(1:4:end,1:4:end)*1e3,vtRc(1:4:end,1:4:end),vzrc(1:4:end,1:4:end),1,'-k')
hold on
plot(Rc.*(theta)*1e3,zeros(1,length(theta)),'-r')
hold on
plot(Rc.*(theta)*1e3,-zRc(N+1,:)*1e3,'-r')
plot(-Rc*phi/2*ones(1,10)*1e3,linspace(-0.75,0.125,10),'--k')
plot(Rc*phi/2*ones(1,10)*1e3,linspace(-0.75,0.125,10),'--k')
text(-4.25,0.55,'(D)','Fontsize',14);
text(-5.25,-0.75,'\beta','Fontsize',14);
text(4.75,-0.75,'\beta*','Fontsize',14);
xlabel('R_c\theta [mm]')
ylabel('z [mm]')
axis([-6,6,-2.5,0.25])

subplot(3,2,5)
surf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,P,'EdgeColor',[0.5,0.5,0.5])
caxis([-4000,4000])
text(0,-9,8000,'(E)','Fontsize',14);
colormap(map6)
view([80,30])
xlabel('X [mm]')
ylabel('Y [mm]')
zlabel('P [Pa]')
axis([0,20,-10,10,-5000,5000])

subplot(3,2,6)
[C,h]=contourf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,P*1e-3);
caxis([-4,4])
v=[0,-2,2];
colormap(map6)
axis([0,20,-10,10])
clabel(C,h,v);
ax=gca;
ax.YAxisLocation='right';
text(-2.5,-8.5,'(F)','Fontsize',14);
camroll(-90)
c=colorbar;
xlabel('X [mm]')
ylabel('Y [mm]')
ylabel(c,'kPa')
zlabel('P [Pa]')
%}
%{
f2=figure('Position', [100, 50, 800, 300]);
subplot(1,3,1)
contourf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,magv0)
hold on
quiver(Rmat(1:4:end,1:4:end).*cos(Theta(1:4:end,1:4:end))*1e3,Rmat(1:4:end,1:4:end).*sin(Theta(1:4:end,1:4:end))*1e3,vx0(1:4:end,1:4:end),vy0(1:4:end,1:4:end),1,'k')
colormap(map6)
caxis([0,20e-3*Omega])
axis([0,20,-10,10])
ax=gca;
ax.YAxisLocation='right';
text(1.5,-8.5,'(A)','Fontsize',14);
camroll(-90)
%colorbar
xlabel('X [mm]')
ylabel('Y [mm]')
zlabel('P [Pa]')

subplot(1,3,2)
contourf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,magvh0d2)
hold on
quiver(Rmat(1:4:end,1:4:end).*cos(Theta(1:4:end,1:4:end))*1e3,Rmat(1:4:end,1:4:end).*sin(Theta(1:4:end,1:4:end))*1e3,vxh0d2(1:4:end,1:4:end),vyh0d2(1:4:end,1:4:end),1,'k')
colormap(map6)
caxis([0,20e-3*Omega])
axis([0,20,-10,10])
ax=gca;
ax.YAxisLocation='right';
text(1.5,-8.5,'(B)','Fontsize',14);
camroll(-90)
%colorbar
xlabel('X [mm]')
ylabel('Y [mm]')
zlabel('P [Pa]')

subplot(1,3,3)
contourf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,magvh0)
hold on
quiver(Rmat(1:4:end,1:4:end).*cos(Theta(1:4:end,1:4:end))*1e3,Rmat(1:4:end,1:4:end).*sin(Theta(1:4:end,1:4:end))*1e3,vxh0(1:4:end,1:4:end),vyh0(1:4:end,1:4:end),1,'k')
colormap(map6)
caxis([0,20e-3*Omega])
axis([0,20,-10,10])
ax=gca;
ax.YAxisLocation='right';
text(1.5,-8.5,'(C)','Fontsize',14);
camroll(-90)
hp3 = get(subplot(1,3,3),'Position');
colorbar('Position', [hp3(1)+hp3(3)+0.01  hp3(2)  0.03  hp3(2)+hp3(3)*3.3])
xlabel('X [mm]')
ylabel('Y [mm]')
zlabel('P [Pa]')

%}
%{
f1=figure();
quiver(Rc.*(Theta(1:4:end,1:4:end))*1e3,-zRc(1:4:end,1:4:end)*1e3,vtRc(1:4:end,1:4:end),vzrc(1:4:end,1:4:end),1,'-k')
hold on
plot(Rc.*(theta)*1e3,zeros(1,length(theta)),'-r')
hold on
plot(Rc.*(theta)*1e3,-zRc(N+1,:)*1e3,'-r')
xlabel('R_c\theta [mm]')
ylabel('z [mm]')
axis([-Rc*phi/2*1e3,Rc*phi/2*1e3,-2.5,0.25])

f2=figure();
contourf(Rc.*Theta*1e3,-zRc*1e3,PRC)
caxis([-4000,4000])
colormap(map6)
xlabel('R_c\theta [mm]')
ylabel('z [mm]')
axis([-Rc*phi/2*1e3,Rc*phi/2*1e3,-2.5,0.25])
%}


%fname=['\Users\Jonathon\Dropbox\Grad Research\2015-03-03-Reynolds_No_Cavitation\Cylindrical\Pressure_Profiles\H_b' num2str(beta*180/pi) '_degrees_total_90deg_rotation.jpg'];
%saveas(f1,fname,'jpeg')
%}
%print(f1,'-dpdf',fname);
%{
figure()
mesh(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,Hmat*1e3)
colormap(map6)
xlabel('--X-- [mm]')
ylabel('--Y-- [mm]')
zlabel('--h(r,\theta)-- [mm]')
axis([0,20,-10,10,0,20])
%}

%{
f1=figure('Position', [50, 50, 700, 450]);
quiver(Rc.*(Theta+phi/2)*1e3,-zRc*1e3,vtRc,vzrc,1.5,'-k')
hold on
plot(Rc.*(theta+phi/2)*1e3,zeros(1,length(theta)),'-r')
hold on
plot(Rc.*(theta+phi/2)*1e3,-zRc(N+1,:)*1e3,'-r')
xlabel('R_c(\theta+\phi/2) [mm]')
ylabel('z [mm]')
axis([0,Rc*phi*1e3,-5,1])
fname=['\Users\Jonathon\Dropbox\Grad Research\2015-03-03-Reynolds_No_Cavitation\Cylindrical\Quiver_Figures\V_Rc_b' num2str(beta*180/pi) '_degrees.jpg'];
saveas(f1,fname,'jpeg')
%axis('equal')

f2=figure('Position', [50, 50, 700, 450]);
quiver(Rmat*1e3,-zTc*1e3,vrtc,vztc,1,'-k')
hold on
plot(r*1e3,zeros(1,length(theta)),'-r')
hold on
plot(r*1e3,-zTc(:,N+1)*1e3,'-r')
xlabel('r [mm]')
ylabel('z [mm]')
axis([R1*1e3,R2*1e3,-5,1])
fname=['\Users\Jonathon\Dropbox\Grad Research\2015-03-03-Reynolds_No_Cavitation\Cylindrical\Quiver_Figures\V_Theta0_b' num2str(beta*180/pi) '_degrees.jpg'];
saveas(f2,fname,'jpeg')
%}
%{
f2=figure('Position', [100, 50, 1000, 1000]);
subplot(2,2,1)
surf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,-Hmat*1e3)
colormap(map6)
view([80,30])
xlabel('--X-- [mm]')
ylabel('--Y-- [mm]')
zlabel('-H(r,\theta) [mm]')
axis([0,20,-10,10,-10,10])

subplot(2,2,2)
surf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,-Hmat*1e3)
hold on
surf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,zeros(N+1))
%colormap([0,0,0])
text(15,6,0,'\rightarrow r\Omega','FontSize',14)
text(15,0,-(h1*1e3)-1,'Fixed','FontSize',14);
view([90,0])
xlabel('--X-- [mm]')
ylabel('--Y-- [mm]')
zlabel('-H(r,\theta) [mm]')
axis([0,20,-10,10,-10,10])

subplot(2,2,4)
quiver(Rc.*(Theta)*1e3,-zRc*1e3,vtRc,vzrc,1.5,'-k')
hold on
plot(Rc.*(theta)*1e3,zeros(1,length(theta)),'-r')
hold on
plot(Rc.*(theta)*1e3,-zRc(N+1,:)*1e3,'-r')
xlabel('R_c\theta [mm]')
ylabel('z [mm]')
axis([-Rc*phi/2*1e3,Rc*phi/2*1e3,-5,1])

subplot(2,2,3)
quiver(Rmat*1e3,-zTc*1e3,vrtc,vztc,1,'-k')
hold on
plot(r*1e3,zeros(1,length(theta)),'-r')
hold on
plot(r*1e3,-zTc(:,N+1)*1e3,'-r')
xlabel('r [mm]')
ylabel('z [mm]')
axis([R1*1e3,R2*1e3,-5,1])

fname=['\Users\Jonathon\Dropbox\Grad Research\2015-03-03-Reynolds_No_Cavitation\Cylindrical\Quiver_Figures\H_b' num2str(beta*180/pi) '_degrees+moving_surf_finer_mesh.jpg'];
saveas(f2,fname,'jpeg')
%}

%{
Nj=1;
X=(Rmat.*cos(Theta)); Y=(Rmat.*sin(Theta));
Ur_t0=zeros(N+1); Ut_rc=zeros(N+1); Z_t0=zeros(N+1); Z_rc=zeros(N+1);
irc=1; while r(irc)<14.25e-3,
    irc=irc+1;
end
figure('Position', [100, 50, 650, 650]);
for i=1:Nj:N+1,
    Z=(1+z(i))/2.*Hmat;
    %[sr,st,sz]=meshgrid(r,-phi/2,(1+z)/2*h0);
    %sx=sr.*cos(st); sy=sr.*sin(st);
    Ur=1./(2*b1).*(2/Rdiff*Dh*PN).*(Z.^2-Z.*Hmat);
    Ut=(1./(2*b1*Rmat).*(2/phi*PN*Dh').*(Z.^2-Z.*Hmat))+(Rmat.*Omega).*((Hmat-Z)./Hmat);
    Uz=zeros(N+1);
    Ux=(Ur.*cos(Theta)-Ut.*sin(Theta));
    Uy=(Ur.*sin(Theta)+Ut.*cos(Theta));
    Ur_t0(:,i)=Ur(:,floor(N/2+1));
    Z_t0(i,:)=Z(floor(N/2+1),:);
    Ut_rc(:,i)=Ut(:,irc);
    Z_rc(i,:)=Z(irc,:);
    %plot3(X(1:Nj:N+1,1:Nj:N+1),Y(1:Nj:N+1,1:Nj:N+1),Z(1:Nj:N+1,1:Nj:N+1),'ro')
    %quiver3(X(1:Nj:N+1,1:Nj:N+1)*1e3,Y(1:Nj:N+1,1:Nj:N+1)*1e3,-Z(1:Nj:N+1,1:Nj:N+1)*1e3,Ux(1:Nj:N+1,1:Nj:N+1),Uy(1:Nj:N+1,1:Nj:N+1),Uz(1:Nj:N+1,1:Nj:N+1),0.75,'k');
    %quiver3(Rmat(1:Nj:N+1,1:Nj:N+1),Theta(1:Nj:N+1,1:Nj:N+1),Z(1:Nj:N+1,1:Nj:N+1),Ur(1:Nj:N+1,1:Nj:N+1),Ut(1:Nj:N+1,1:Nj:N+1),Uz(1:Nj:N+1,1:Nj:N+1),'k');
    %quiver3(X,Y,Z,Ux,Uy,Uz,0.5,'k');
    %stream(X,Y,Z,Ur,Ut,Uz,sx,sy,sz);
    hold on
end
axis([0,20,-10,10,-5,2])
xlabel('X [mm]')
ylabel('Y [mm]')
zlabel('Z [mm]')

%}

%{
X3=zeros(N+1,N+1,N+1); Y3=zeros(N+1,N+1,N+1); Z3=zeros(N+1,N+1,N+1);
Ur3=zeros(N+1,N+1,N+1); Ut3=zeros(N+1,N+1,N+1); Uz3=zeros(N+1,N+1,N+1);
Ux3=zeros(N+1,N+1,N+1); Uy3=zeros(N+1,N+1,N+1);
for i=1:N+1,
    X3(:,:,i)=((Rmat').*cos(Theta'));
    Y3(:,:,i)=((Rmat').*sin(Theta'));
    Z3(:,:,i)=(-(1+z(i))/2.*Hmat);
    Ur3(:,:,i)=1./(2*b1).*(2/Rdiff*Dh*PN).*(Z3(:,:,i).^2-Z3(:,:,i).*Hmat);
    Ut3(:,:,i)=1./(2*b1*Rmat).*(2/phi*Pn*Dh')*(Z3(:,:,i).^2-Z3(:,:,i).*Hmat)+(Rmat.*Omega).*((Hmat-Z3(:,:,i))./Hmat);
    Ux3(:,:,i)=(Ur3(:,:,i).*cos(Theta)-Rmat.*Ut3(:,:,i).*sin(Theta));
    Uy3(:,:,i)=(Ur3(:,:,i).*sin(Theta)+Rmat.*Ut3(:,:,i).*cos(Theta));
end
figure()
for i=1:3:N+1,
    [sr,st,sz]=meshgrid(r,-phi/2,(1+z(1:4:N+1))/2.*Hmat(1:4:N+1,i));
    sx=(sr.*cos(st)); sy=(sr.*sin(st));
    %plot3(sx(:),sy(:),sz(:),'r*')
    streamline(X3,Y3,Z3,Ux3,Uy3,Uz3,sx,sy,sz);
    hold on
end
xlabel('X')
ylabel('Y')
axis(1e-3*[0,20,-10,10,-5,5])
%}
     
    
%{
invgray=zeros(264,3);
for j=1:264,
    invgray(j,3)=1.000000000000000-0.002272727272727*j-0.000005739210285*j^2;
    if invgray(j,3)<0,
        invgray(j,3)=0;
    end
    invgray(j,2)=1.000000000000000-0.002272727272727*j-0.000005739210285*j^2;
    if invgray(j,2)<0,
        invgray(j,2)=0;
    end
    invgray(j,1)=1.000000000000000-0.002272727272727*j-0.000005739210285*j^2;
    if invgray(j,1)<0,
        invgray(j,1)=0;
    end
end


figure()
contourf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,P)
colormap(invgray)
colorbar
caxis([-4000,4000])
axis([0,20,-10,10])
%ax=gca;
%ax.YAxisLocation='right';
%camroll(-90)
xlabel('--X-- [mm]')
ylabel('--Y-- [mm]')
zlabel('--P-- [Pa]')
%}

%{
map6=zeros(264,3);
for j=1:264,
    map6(j,3)=1.003065045133951-0.003065133088643*j+0.000000087954693*j^2;
    if map6(j,3)>1,
        map6(j,3)=1;
    end
    map6(j,2)=-0.007259164659101+0.007279460204459*j-0.000020295545359*j^2;
    map6(j,1)=-0.005358023974690+0.005363916939106*j-0.000005892964416*j^2;
    if map6(j,1)>1,
        map6(j,1)=1;
    end
    %
    if map6(j,3)<0,
        map6(j,3)=0;
    end
    if map6(j,2)<0,
        map6(j,2)=0;
    end
    if map6(j,1)<0,
        map6(j,1)=0;
    end
end

%}



%{
figure('Position', [100, 50, 750, 1000]);
subplot(2,2,1)
surf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,-Hmat*1e3)
colormap(map6)
view([80,30])
text(0,-9,4,'(A)','Fontsize',14);
xlabel('--X-- [mm]')
ylabel('--Y-- [mm]')
zlabel('--h(r,\theta)-- [mm]')
axis([0,20,-10,10,-5,5])

subplot(2,2,2)
surf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,-Hmat*1e3)
hold on
mesh(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,zeros(N+1))
%colormap([0,0,0])
text(15,6,0,'\rightarrow r\Omega','FontSize',14)
text(15,0,-(2.75e-3*1e3)-0.5,'Fixed','FontSize',14);
%text(0,-9,3.5,'(D)','Fontsize',14);
text(0,-9,4.25,'(B)','Fontsize',14);
view([80,5])
xlabel('X [mm]')
ylabel('Y [mm]')
zlabel('-H(r,\theta) [mm]')
axis([0,20,-10,10,-5,5])

subplot(2,2,3)
surf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,P)
colormap(map6)
caxis([-4000,4000])
view([80,30])
axis([0,20,-10,10,-4000,4000])
xlabel('--X-- [mm]')
ylabel('--Y-- [mm]')
zlabel('--P-- [Pa]')
text(0,-9,3250,'(C)','Fontsize',14);

subplot(2,2,4)
contourf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,P)
colormap(map6)
colorbar
caxis([-4000,4000])
axis([0,20,-10,10])
ax=gca;
ax.YAxisLocation='right';
camroll(-90)
xlabel('--X-- [mm]')
ylabel('--Y-- [mm]')
zlabel('--P-- [Pa]')
text(2,-9,'(D)','Fontsize',14);
%}
%{
figure()
surf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,P)
%colormap(map6)
view([80,30])
xlabel('--X-- [mm]')
ylabel('--Y-- [mm]')
zlabel('--P-- [Pa]')
%

figure()
contourf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,P)
colormap(map6)
colorbar
axis([0,20,-10,10])
ax=gca;
ax.YAxisLocation='right';
camroll(-90)
xlabel('--X-- [mm]')
ylabel('--Y-- [mm]')
zlabel('--P-- [Pa]')
%
figure()
surf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,Hmat*1e3)
colormap(map6)
xlabel('--X-- [mm]')
ylabel('--Y-- [mm]')
zlabel('--h(r,\theta)-- [mm]')
axis([0,20,-10,10,0,20])
%}