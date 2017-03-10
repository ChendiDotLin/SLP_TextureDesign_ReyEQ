function [Fn,tstar]=Reynolds_Tex2(H,opt)
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

format long

%Define variables used to generate the mesh
R1 = opt.R1;        % Inner radius [1e-3]
R2 = opt.R2;        % Outer radius [20e-3]
phi = opt.phi;      % Angle of sector [2 * pi / 10]
h0 = opt.h0;        % Nominal gap distance [0.25e-3 + 19e-6]
Omega = opt.Omega;  % Angular velocity [10]
b1 = opt.b1;        % Viscosity [1.4]
b2 = opt.b2;        % Second-order fluid parameter, [0]
Ntex = 2*pi/phi;    % Number of periodic sectors, [10]
N = size(H,1)-1;    % Polynomial order, [25]
%N = 3;
%Create the 1-D matrices for the psuedo-spectral method

[Kh,Mh,Ch,Dh,z,w] = semhat(N);
%[~,~,~,Dh2,~,~] = semhat(N-1);
q = zeros(1,N);
q(N) = 1;
Q = [q;eye(N)];

I = eye(N+1);
Res = I(2:N,:);
Pro = Res';

Nuem = I(1:N,:)';
r = (R2-R1)/2*(z) + (R2+R1)/2;
theta = phi/2*(z);
[Rmat,Theta] = ndgrid(r,theta);
Rdiff = R2-R1;

R1d = diag(r);
R = kron(I,R1d);
Rinv1d = diag(1./r);
Rinv = kron(I,Rinv1d);

dhdz = H*Dh';
Hmat = H;

dh = reshape(dhdz,numel(dhdz),1);
h = reshape(H,numel(H),1);

%Create the 2-D matrices used for solving the Reynolds equation

H = diag(h.^3);
Kb = -phi/Rdiff*(kron(I,Dh))'*(kron(Mh,Mh)*R*H)*kron(I,Dh)...
    - Rdiff/phi*(kron(Dh,I))'*(kron(Mh,Mh)*H*Rinv)*kron(Dh,I);
fb = 3*b1*Rdiff*Omega*kron(Mh,Mh)*R*dh;

K = (kron(Q,Nuem))'*Kb*kron(Q,Nuem);
K = sparse(K);
f = kron(Q,Nuem)'*fb;
f = sparse(f);

%Solve for the Newtonian pressure
p=K\f;


%Create the Newtonian Pressure matrix
Pn=zeros(N,N);
for j=1:N,
    for i=1:N,
        Pn(i,j)=p((i)+(j-1)*(N));
    end
end
Pn=Nuem*Pn*Q';


%Calculate the Non-Newtonian Pressure matrix
P=Pn+b2./b1.*Omega.*(2/phi*Pn*Dh')+b2./2*((1./(2*b1).*(2/Rdiff*Dh*Pn).*Hmat)^2+...
    (1./(2*b1.*Rmat).*(2/phi*Pn*Dh').*Hmat).^2+(Omega*Rmat./Hmat)^2+...
    Omega/b1*(2/phi*Pn*Dh'));
%Calculate tauzz
tauzz=b2*((1./(2*b1).*(2/Rdiff*Dh*Pn).*Hmat)^2+...
    (1./(2*b1.*Rmat).*(2/phi*Pn*Dh').*Hmat).^2+(Omega*Rmat./Hmat)^2+...
    Omega/b1*(2/phi*Pn*Dh'));
%Compute Normal Force
Fn=Ntex*(Rdiff*phi/4*w'*(R1d'*(P-tauzz))*w);
%Compute torque and viscosity
tau=-(-1./(2*Rmat).*(2/phi*Pn*Dh').*Hmat-b1.*Rmat.*Omega./(Hmat)-...
    b2.*Omega./(2*b1.*Rmat).*(2/phi*(2/phi*Pn*Dh')*Dh').*Hmat+...
    b2.*Rmat.*Omega^2./(Hmat.^2).*(2/phi*Hmat*Dh'));

M=Ntex*(Rdiff*phi/4*w'*(R1d'.^2*tau)*w);
eta_a=2*h0/(pi*R2^4)*M/Omega;
tstar=(2/(pi*R2^3)*M)/(b1*(Omega*R2)/h0);
%Compute effective friction coefficient

%pltP; drawnow; print('fig_pressure','-depsc'); %caxis([-1.125 +8.149]);
%pltTau; drawnow; print('fig_tau','-depsc'); %caxis([-0.206 +1.072]);
%pltH; drawnow; print('fig_geom','-depsc'); %caxis([-1.5 -0.269]);

mu=(M/R2)/Fn;
