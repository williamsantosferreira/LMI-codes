%{
Polytopic + D-stability region
%}
clear;
clc;
H = tf(9,[1 4.8 9]); % process
[nH,dH] = tfdata(H,'v');
[A,B,C,D] = tf2ss(nH,dH);

delta = [0.2 0.4
         0.1  0.2];
A1 = A + delta;
A2 = A - delta;

G1 = ss(A1,B,C,D);
G2 = ss(A2,B,C,D);

Z = sdpvar(1,length(A));
P = sdpvar(length(A));
alfa = -2;
beta =-4; 

F0 = [P >= 0, Z >= 0];
F1 = [-2*alfa*P + A1*P + P*A1' + B*Z + Z'*B' <= 0]; %linha vertical < alfa
F2 = [2*beta*P + A1*P + P*A1' + B*Z + Z'*B' <= 0]; %linha vertical > beta
F3 = [-2*alfa*P + A2*P + P*A2' + B*Z + Z'*B' <= 0];
F4 = [2*beta*P + A2*P + P*A2' + B*Z + Z'*B' <= 0];
F = [F0,F1,F2,F3,F4]; 
sol = solvesdp(F)
Z = value(Z);
P = value(P);

K = Z*inv(P);

Amf1 = A1 + B*K;
Amf2 = A2 + B*K;

Gmf1 = ss(Amf1,B,C,D);
Gmf2 = ss(Amf2,B,C,D);

U1 = ss(Amf1,B,K,0);
U2 = ss(Amf2,B,K,0);

rlocus(Gmf1);
figure;
rlocus(Gmf2);
figure;

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


subplot(2,3,1);
step(G1);
title('G1');

subplot(2,3,2);
step(Gmf1);
title('G1 closed loop');

subplot(2,3,3);
step(U1);
title('U1 control signal');

subplot(2,3,4);
step(G2);
title('G2');

subplot(2,3,5);
step(Gmf2);
title('G2 closed loop');

subplot(2,3,6);
step(U2);
title('U2 control signal');





