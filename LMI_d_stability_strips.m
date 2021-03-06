%{
D-stability code

A região D-stability anilizada abaixo são de duas linhas verticais.
Ou seja, é necessário que:
--> Seja encontrado polos dentro dessa região especificada

Teoria:

%}

clear;
clc;
H = tf(9,[1 4.8 9]); % process
[nH,dH] = tfdata(H,'v');
[A,B,C,D] = tf2ss(nH,dH);
G = ss(A,B,C,D);

Z = sdpvar(1,length(A));
P = sdpvar(length(A));
alfa = -2;
beta =-4; 

F1 = [P >= 0, Z >= 0];
F2 =[-2*alfa*P + A*P + P*A' + B*Z + Z'*B' <= 0]; %linha vertical < alfa
%F3 =[2*beta*P + A*P + P*A' + B*Z + Z'*B' <= 0]; %linha vertical > beta
F = [F1,F2]; 
sol = solvesdp(F)
Z = value(Z);
P = value(P);

K = Z*inv(P);

Amf = A + B*K;
Gmf = ss(Amf,B,C,D);
U = ss(Amf,B,K,0);

[Z,P,K] = ss2zp(Amf,B,C,D);
rlocus(Gmf);
figure
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

subplot(1,3,1);
step(G);
title('Open Loop');
grid;

subplot(1,3,2);
step(Gmf);
title('Closed Loop');
grid;

subplot(1,3,3);
step(U);
title('Control Signal');
grid;





