%LMI State Integral Discrete
Ts = 0.01; %tempo de amostragem
r = 1; %referência

H = tf(1,[1 -1 2]);
[nH,dH] = tfdata(H,'v');
[A,B,C,D] = tf2ss(nH,dH);

G = ss(A,B,C,D);%contínuo
G = c2d(G,Ts);

[A,B,C,D] = ssdata(G);

Aa = [A zeros(2,1)
      -C 1]; %porque a estrutura é diferente?
Bu = [B
      0];

Q = sdpvar(length(Aa));
N = sdpvar(1,length(Bu));

F = [Q >= 0, [-Q Q*Aa'+N'*Bu';Aa*Q+Bu*N -Q] <= 0];

sol = solvesdp(F);
Q = value(Q);
N = value(N);

Ka = N*inv(Q);

Amf = Aa + Bu*Ka;
Bmf = [0 0 1]'*r;
Cmf = [C 0];

Gmf = ss(Amf,Bmf,Cmf,D,Ts);
G = ss(A,B,C,D,Ts);
U = ss(Amf,Bmf,Ka,D,Ts);

subplot(1,3,1);
step(G);
title('Open Loop');
grid;

subplot(1,3,2);
step(Gmf);
title(['Closed Loop with reference equal to ',num2str(r)])
grid;

subplot(1,3,3);
step(U);
title('Control Signal');
grid;

