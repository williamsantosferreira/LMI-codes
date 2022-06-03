H = tf(1,[1 0.1 2]); % process
[nH,dH] = tfdata(H,'v');
[A,B,C,D] = tf2ss(nH,dH);
G = ss(A,B,C,D);

r = 2; %referência

Aa = [A zeros(2,1)
      -C 0];
Bu = [B
      0];

Q = sdpvar(length(Aa));
N = sdpvar(1,length(Bu));

F = [Q >= 0, Q*Aa'+ Aa*Q + N'*Bu' + Bu*N <=0];
solvesdp(F)

Q = value(Q);
N = value(N);

Ka = N*inv(Q);

Amf = Aa + Bu*Ka;
Bmf = [0 0 1]'*r;
Cmf = [C 0];

Gmf = ss(Amf,Bmf,Cmf,D);
G = ss(A,B,C,D);
U = ss(Amf,Bmf,Ka,D);

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


