%{
Polytopic Approach
%}
clear;
clc;
H = tf(1,[1 0.1 2]); % process
[nH,dH] = tfdata(H,'v');
[A,B,C,D] = tf2ss(nH,dH);
G = ss(A,B,C,D);

delta = [0.5 1
         1.1  0.2];
A1 = A + delta;
A2 = A - delta;

G1 = ss(A1,B,C,D);
G2 = ss(A2,B,C,D);

Q = sdpvar(length(A));
Y = sdpvar(1,length(B));

F0 = [Q >= 0];
F1 = [Q*A1' + A1*Q + Y'*B' + B*Y <= 0];
F2 = [Q*A2' + A2*Q + Y'*B' + B*Y <= 0];

F = [F0,F1,F2];
sol = solvesdp(F)
Q = value(Q);
Y = value(Y);

K = Y*inv(Q)

Amf1 = A1 + B*K;
Amf2 = A2 + B*K;

Gmf1 = ss(Amf1,B,C,D);
Gmf2 = ss(Amf2,B,C,D);

U1 = ss(Amf1,B,K,0);
U2 = ss(Amf2,B,K,0);

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
     