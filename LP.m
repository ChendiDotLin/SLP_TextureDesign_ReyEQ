% In this program, we still use normal force as the objective function, and
% use the linearized tau as the inequality constraint
% MOSEK Toolbox is needed in this program

function[Fopt,taustaropt1] = LP(alpha,taunow)
% given step size and taustar we want
clc

% setup the path of solver and the optimizer
currentpath = pwd;
if ispc; dirsep = '\'; elseif isunix; dirsep = '/'; end;
path(path,strcat(currentpath, dirsep, 'Solver'));

% make it sequential
% design variable x = [h;dp]
% now h is N*(N+1), so is dpdtheta

% This path is used for mac
%addpath '/Users/apple/mosek/7/toolbox/r2013aom/'

% This path is used for windows
addpath 'C:\Program Files\Mosek\8\toolbox\r2014aom'

format long; 
opt.R1 = 0.01e-3;
opt.R2 = 20e-3;
opt.phi = 2 * pi / 10;
opt.h0 = 0.25e-3 + 19e-6;
opt.hmax = 2e-3;
opt.Omega = 10;
opt.b1 = 1.4;
opt.b2 = 0; % Newtonian
opt.N =25; % dimension of H
N = opt.N;
R1 = opt.R1;        % Inner radius [1e-3]
R2 = opt.R2;        % Outer radius [20e-3]
phi = opt.phi;      % Angle of sector [2 * pi / 10]
h0 = opt.h0;        % Nominal gap distance [0.25e-3 + 19e-6]
hmax = opt.hmax;
Omega = opt.Omega;  % Angular velocity [10]
b1 = opt.b1;        % Viscosity [1.4]
b2 = opt.b2;        % Second-order fluid parameter, [0]
Ntex = 2*pi/phi;    % Number of periodic sectors, [10]
region = 0.1; reg_for_dp =100; % set trust region for both h and dp
error = 0.01;       % initialize error
count = 1;          % Number of iterations

% starting point
h_opt = load('startpoint.mat');
Hopt1 = h_opt.Hopt;
hopt = reshape(Hopt1,26^2,1);
hopt = hopt(27:26*26);

% check if it is out of boundary
hopt = max(h0*ones(N*(N+1),1),hopt);
hopt = min(hmax*ones(N*(N+1),1),hopt);
%because of boundary condition, its dimension is (N+1)*N

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

dr = kron(I,Dh);                                    % take derivative of p with respect to r
dtheta = kron(Dh,I);                                % take derivative of p with respect to theta
dd = [dr;dtheta];                                   % take derivative of lagrange polynomial of the whole p
intedp = (pinv(dd'*dd))*(dd');                      % integrate dp to find out p
intep = ((kron(w,w)).*(kron((ones(N+1,1)),r)))';    % integrate p to find out normal force

h = zeros((N+1)*N);
dpopt = 1.1e6*ones(2*(N+1)^2,1);
ct = [ zeros(1,length(h)),intep*intedp];            % integral of dp twice to find out normal force 
% design variable x = [h;dp]

% find normal force values and taustar values of the starting point
[Fn_initial,tau_initial] = Reynolds_Tex4(reshape(kron(Q,I)*hopt,N+1,N+1),0,10);
taustaropt = tau_initial;
fopt = Fn_initial;

while(error>=0.002)
    count = count+1;
    hopt = max(h0*ones(N*(N+1),1),hopt); hopt = min(hmax*ones(N*(N+1),1),hopt);
    % check boundary condition
    
    hguessstar = hopt; %hguess is the reference texture, updated at the beginning of the loop
    hguess = kron(Q,I)*hguessstar; % apply periodic boundary condition
    fguess = fopt;
    
    
    % equality constraints
    
    %kp = f
    K1 = constantK1noSparse(opt,hguess);
    f1 = constantf(opt);
    Fmat = [f1*kron(Q,I), zeros((N+1)^2,2*(N+1)^2)];% matrix on the right
    DQ1 = diag(dpopt);
    I2 = eye((N+1)^2);
    VAGUE = zeros((N+1)^2,(N+1)^2);
    
    % to find pr and ptheta at this point, we still need the solver
    H0 = reshape(hguess,(N+1),(N+1));                   % original h
    H00 = diag(hguess);                                 % diagonal matrix of h
    [pref2,dp0,~,~,~,~,~,~] = Reynolds_Tex3(H0,0,10);   % here we found p0 and dp0
    
    dpdr0 = dp0(1:(N+1)^2);
    dpdtheta0 = dp0((N+1)^2+1:2*(N+1)^2);
    
    PR = diag(dpdr0);
    PTHETA = diag(dpdtheta0);
    H00SQ = diag(hguess.^2);
    H00CU = diag(hguess.^3);
    HELPER = [3*H00SQ*PR*kron(Q,I),H00CU,VAGUE;3*H00SQ*PTHETA*kron(Q,I),VAGUE,H00CU];
    A1eq = (kron(Q,I))'*[K1*HELPER-Fmat];
    % setting up Kp = f
    
    
    
    % new A2eq
    K2 = zeros((N+1),(N+1)^2);
    for i = 1:N+1
        K2(i,i) = 1;
        K2(i,i+N*(N+1)) = -1;
    end
    
    A2eq = [zeros((N+1),(N+1)*N),K2*intedp];
    % apply periodic boundary to p, which is the integral of dp
    Aeq = [A1eq;A2eq];
    
    %Aeq = [-f,K];
    
    HELPER2 = 3*[hguess.^3;hguess.^3].*[dpdr0;dpdtheta0]; % make a change here
    b1eq = (kron(Q,I))'*K1*HELPER2;
    b2eq = zeros((N+1),1);
    beq = [b1eq;b2eq];
    
    
    % objective function
    c1  =  (h0/(pi*(R2)^4*Omega*b1))*(Rdiff*phi/2)*Ntex;
    rmat = reshape(Rmat,(N+1)^2,1);
    RSQU = diag(rmat.^2);
    A1 = (kron(w,w))'*RSQU;
    B1 = kron(Q,I);
    HELPERTAU = [(diag((1/phi*dpdtheta0./rmat)-b1*Omega*rmat./(hguess.^2)))*B1,VAGUE,(diag(1/phi*hguess./rmat))];
    
    % here is my inequality constraint: taustar should be less than some value
    Ain1 = c1*A1*HELPERTAU;
    %f = zeros(1,N*(N+1)+2*(N+1)^2);
    HELPERConst = -1/phi*dpdtheta0 .* hguess./rmat + 2*b1*Omega*rmat./hguess;
    
    % set the taustar value we want here
    bin1 = taunow-c1*A1*HELPERConst;
    bin3 = -taunow + 0.01 + c1*A1*HELPERConst;

    
    % the objective function should be normal force 
    
    % normal force linear part
    HelperIn = pref2*ones((N+1)^2,1);
    normalforce = Ntex*(Rdiff/2)*(phi/2)*(kron(w,w))'*diag(rmat)*intedp;
    objective = [zeros(1,(N+1)*N),normalforce];
    f = -objective;
    % normalforce*dp0
    
    offsetforce = Ntex*(Rdiff/2)*(phi/2)*(kron(w,w))'*diag(rmat)*HelperIn;
    
    % another inequality constraint
    prefhelper = zeros((N+1),(N+1)^2);
    for i = 1:N+1
        prefhelper(i,i*(N+1)) = 1;
    end
    
    % this inequality constraint make sure the constant term is not smaller
    % than previous one
    Ain2 = [zeros(1,N*(N+1)),-0.5*w'*prefhelper*intedp];
    bin2 = -0.5*w'*prefhelper*intedp*dp0;
    

    
    % boundary condition
    glbforh = h0*ones((N+1)*N,1);
    gubforh = hmax*ones((N+1)*N,1);
    % local bounds for h
    llbforh = (1-region)*hguessstar; lubforh = (1+region)*hguessstar;
    % local bounds for dp
    llbfordp = dp0-reg_for_dp; lubfordp = dp0+reg_for_dp;
    
    lbforh = max(glbforh,llbforh);
    ubforh = min(gubforh,lubforh);
    
    lb = [lbforh;llbfordp];
    ub = [ubforh;lubfordp];
    
    % using MOSEK
    option = mskoptimset('');
    % Turn on diagnostic output
    %option = mskoptimset(option,'Diagnostics','on');
    % Set a MOSEK option, in this case turn basic identification off.
    option = mskoptimset(option,'MSK_IPAR_INTPNT_BASIS','MSK_OFF');
    % Modefy a MOSEK parameter with double value
    option = mskoptimset(option,'MSK_DPAR_INTPNT_TOL_INFEAS',1e-12);
    %}
    [xopt,fval,exitflag,~,lambda] = linprog(f,[Ain1;Ain2],[bin1;bin2],A1eq,b1eq,lb,ub,[],option);
    disp(strcat('exitflag:',num2str(exitflag)));
    disp(strcat('max_lambda (lower, upper):',num2str(max(abs(lambda.lower))),'/',num2str(max(abs(lambda.upper)))));
    disp(strcat('offset:',offsetforce));
     
    dpdropt = zeros((N+1)^2,1);
    dpdthetaopt = zeros((N+1)^2,1);
    hopt = zeros((N+1)*N,1);
    for i = 1:(N+1)*N
        hopt(i) = xopt(i);
        dpdropt(i) = xopt(i+(N+1)*N);
        dpdthetaopt(i) = xopt(i+(N+1)^2+(N+1)*N);
    end
    for i = 1:N+1
        dpdropt(i+(N+1)*N) = xopt(2*(N+1)*N+i);
        dpdthetaopt(i+(N+1)*N) = xopt(2*(N+1)*N+(N+1)^2+i);
    end
    
    [fopt,tauopt] = Reynolds_Tex4(reshape(kron(Q,I)*hopt,N+1,N+1),0,10);
    % get fopt from solver
    
    % if normal force is smaller, then apply step size
    % set up constants of step size
    n1 = 0.01; n2 = 0.9; r1 = 0.5; r2 = 0.5;
    fopt_esti = [zeros(1,(N+1)*N),normalforce]*xopt + Ntex*(Rdiff/2)*(phi/2)*(kron(w,w))'*diag(rmat)*HelperIn;
    tau_esti = Ain1*xopt+c1*A1*HELPERConst;
    pk = (fguess - fopt)/(fguess-fopt_esti);
%{
    if ((pk < n1) && (abs(tauopt-taunow) > 0.05))
        s = hopt - hguessstar;
        hopt = hguessstar + alpha*s;
        [fopt,~] = Reynolds_Tex4(reshape(kron(Q,I)*hopt,N+1,N+1),0,10);
    end
    %}
    
    if (pk < n1)
        s = hopt - hguessstar;
        hopt = hguessstar + alpha*s;
        [fopt,~] = Reynolds_Tex4(reshape(kron(Q,I)*hopt,N+1,N+1),0,10);
    end
    %}
    %{
    if (pk >= n2)
        region = region *1.5; reg_for_dp = reg_for_dp*1.5;
    elseif (pk >=n1)
        region = region*0.75; reg_for_dp = reg_for_dp*0.75;
    else 
        region = r1*region; reg_for_dp = r1*reg_for_dp;
    end
    %}
 

    % Update the error in this iteration
    %error = abs(fopt_esti-fopt)/abs(fopt)
    error = abs(tauopt-tau_esti)
    
    % If maximum iteration reached
    if (count>40)
        break;
    end
    
end
% while loop end

% post process
% hopt = hguessstar;
% apply boundary condition to h
hopt2 = kron(Q,I)*hopt;
Hopt = reshape(hopt2,opt.N+1,opt.N+1);
plotH(Hopt);
%normal force and taustar from the solver
[Fopt,taustaropt1] = Reynolds_Tex4(Hopt,0,10);

end