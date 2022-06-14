close all;
clear;
clc;

M = 0.112; %Massa da base (Rodas e Motores)
m = 0.438; %Massa do pêndulo
b = 0.1; %Atrito
I = 0.00287; %Momento de Inercia
g = 9.8; %Gravidade
l = 0.03765; %Distância da base ao centro de massa

q = [(M+m)*(I+m*l^2)-(m*l)^2];
s = tf ('s');
H = ((m*l*s)/q)/(s^3+b*(I+m*l^2)*s^2/q-(M+m)*m*g*l*s/q-b*m*g*l/q);
[nH,dH] = tfdata(H,'v');
[A,B,C,D] = tf2ss(nH,dH);
G = ss(A,B,C,D);

x1 = -0.1;
y1 = 1.33;
y2 = -1.33;

P = sdpvar(length(A));
Z = sdpvar(1,length(B));
Q = sdpvar(length(A));
Y = sdpvar(1,length(B));

F0 = [P >= 0, Q >= 0];
F2 = [Q*A'+ A*Q + Y'*B' + B*Y <=0];
F = [F0,F2];

sol = solvesdp(F);

%Z = value(Z);
P = value(P);
Y = value(Y);
Q = value(Q);

K = Y*inv(Q);

Amf = A + B*K;
Gmf = ss(Amf,B,C,D);
U = ss(Amf,B,K,0);

rlocus(Gmf);
figure
[Z,P,K] = ss2zp(Amf,B,C,D);

%% checking solution
if sol.problem == 0
    [primal,~]=check(F); % Checking that the solver returned a proper solution
    if (min(primal)>=0 && all(primal(1)>0))
        disp('Sucessfully solved LMIs without problems');
        
    else
        disp('LMIs not solved');
    end
else
    [primal,~]=check(F);
    if (min(primal)>=0 && all(primal(1)>0))
        disp(['Sucessfully solved LMIs, but solver acused ' yalmiperror(sol.problem)]);
    else
        disp(['LMIs not solved. Solver acused ' yalmiperror(sol.problem)]);
    end
end
opt = stepDataOptions('InputOffset',0,'StepAmplitude',1);
subplot(1,3,1);
step(G,opt);
title('Open Loop');
grid;

subplot(1,3,2);
step(Gmf,opt);
title('Closed Loop');
grid;

subplot(1,3,3);
step(U,opt);
title('Control Signal');
grid;
