%{
State Feedback + Integral Action + Polytopic conditions
%}
H = tf(1,[1 0.1 2]); % process
[nH,dH] = tfdata(H,'v');
[A,B,C,D] = tf2ss(nH,dH);
G = ss(A,B,C,D);

delta = [0.2 0.4
         0.1  0.2];
A1 = A + delta;
A2 = A - delta;

r = 2; %referência

Aa1 = [A1 zeros(2,1)
      -C 0];
Aa2 = [A2 zeros(2,1)
       -C 0];
Bu = [B
      0];

Q = sdpvar(length(Aa1));
N = sdpvar(1,length(Bu));

F0 = [Q >= 0];
F1 = [Q*Aa1'+ Aa1*Q + N'*Bu' + Bu*N <=0];
F2 = [Q*Aa2'+ Aa2*Q + N'*Bu' + Bu*N <=0];

F = [F0,F1,F2];
sol = solvesdp(F)

Q = value(Q);
N = value(N);

Ka = N*inv(Q);

Amf1 = Aa1 + Bu*Ka;
Amf2 = Aa2 + Bu*Ka;

Bmf = [0 0 1]'*r;
Cmf = [C 0];

Gmf1 = ss(Amf1,Bmf,Cmf,D);
Gmf2 = ss(Amf2,Bmf,Cmf,D);

G1 = ss(A1,B,C,D);
G2 = ss(A2,B,C,D);
U1 = ss(Amf1,Bmf,Ka,D);
U2 = ss(Amf2,Bmf,Ka,D);

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


