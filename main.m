clc; close all;
%taustar = [0.25,0.3,0.35,0.4,0.45,0.47,0.5,0.55,0.6,0.63,0.65,0.7];
taustar = linspace(0.2,0.7,15);
alpha = 0.01;
F = [];
tau = [];
for i=1:length(taustar)
    [Fopt,tauopt] = LP(alpha,taustar(i));
    F(i) = Fopt;
    tau(i) = tauopt;
end
plot(tau,F,'-r');
title('Without trust region adjusting');
xlabel('\tau ^*')
ylabel('Normal Force')
figure()
Plotpareto;