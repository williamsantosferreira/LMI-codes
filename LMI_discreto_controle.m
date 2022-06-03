%{
Código: LMI Caso Discreto

Aplicar a propriedade de Lyapunov para encontrar um possível controlar que
permita a estabilidade do meu sistema.

A'PA - P < 0
P > 0

Porém, agora meu sistema é do tipo:
dx(t) = Ax(t) + Bu(t)
u(t) = -Kx(t)

[  -Q   QA'+N'B'
  AQ+BN   -Q] < 0

Q = P^-1
N = KQ
K = NQ^-1
%}
Ts = 0.1; %tempo de amostragem
H = tf(1,[1 -1 2]);
[nH,dH] = tfdata(H,'v');
[A,B,C,D] = tf2ss(nH,dH);

G = ss(A,B,C,D);%contínuo
G = c2d(G,Ts);

[A,B,C,D] = ssdata(G);


Q = sdpvar(length(A));
N = sdpvar(1,length(B));

F = [Q >= 0, [-Q Q*A'+N'*B';A*Q+B*N -Q] <= 0];

sol = solvesdp(F);
Q = value(Q);
N = value(N);

K = N*inv(Q);

Amf = A + B*K;
Gmf = ss(Amf,B,C,D,Ts);
U = ss(Amf,B,K,0,Ts);


subplot(1,3,1);
step(G);
title('Open Loop');
set(findall(gcf,'type','line'),'linewidth',3);
grid;

subplot(1,3,2);
step(Gmf);
title('Closed Loop');
grid;

subplot(1,3,3);
step(U);
title('Control Signal');
grid;