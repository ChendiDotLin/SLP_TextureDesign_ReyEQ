function []=plotH(H, varargin)

format long
R1=1e-3; R2=20e-3; h0=0.25e-3+19e-6; Rc=14.25e-3;
phi=2*pi/10; Ntex=2*pi/phi;
b2=0;% b2=-5e-6 for 0.5wt%PIB in S6
A=size(H); N=A(1)-1;

[z,w]=zwgll(N);
r=(R2-R1)/2*(z)+(R2+R1)/2;
theta=phi/2*(z);
[Rmat,Theta]=ndgrid(r,theta);
Rdiff=R2-R1;

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

surf(Rmat.*cos(Theta)*1e3,Rmat.*sin(Theta)*1e3,-H*1e3); hold on;
colormap(map6);
xlabel('X [mm]');
ylabel('Y [mm]');
zlabel('-H(r,\theta) [mm]');
axis([0,20,-10,10,-10,10]);

if (nargin == 2)
    nr = length(r);
    nt = length(theta) - 1;
    for i = 1:nr
        for j = 1:nt
            plot3(r(i)*cos(theta(j))*1e3,r(i)*sin(theta(j))*1e3,-H(i,j)*1e3,...
                'o','MarkerEdgeColor',[1 0 0],'LineWidth',2,...
                'MarkerFaceColor','none'); hold on;
        end
    end
end
hold off;