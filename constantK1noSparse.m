function [K] = constantK1noSparse(opt,hguess)
format long

%Define variables used to generate the mesh
R1 = opt.R1;        % Inner radius [1e-3]
R2 = opt.R2;        % Outer radius [20e-3]
phi = opt.phi;      % Angle of sector [2 * pi / 10]
h0 = hguess;        % reference texture
Omega = opt.Omega;  % Angular velocity [10]
b1 = opt.b1;        % Viscosity [1.4]
b2 = opt.b2;        % Second-order fluid parameter, [0]
Ntex = 2*pi/phi;    % Number of periodic sectors, [10]
N = opt.N;    % 

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

h = hguess;
H0 = diag(h.^3);

K = [-phi/Rdiff*(kron(I,Dh))'*(kron(Mh,Mh)*R), - Rdiff/phi*(kron(Dh,I))'*(kron(Mh,Mh)*Rinv)];


%K = (kron(Q,Nuem))'*Kb*kron(Q,Nuem);
%K = sparse(K)

