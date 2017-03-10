function [f] = constantf(opt) 

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
N = opt.N;    % Polynomial order, [25]

%Create the 1-D matrices for the psuedo-spectral method

[Kh,Mh,Ch,Dh,z,w] = semhat(N);

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

 %Create the 2-D matrices used for solving the Reynolds equation


f = 3*b1*Rdiff*Omega*kron(Mh,Mh)*R*kron(Dh,I);

%f = kron(Q,Nuem)'*fb;
%f = sparse(f);

