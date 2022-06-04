%{
Código: LMI Caso Contínuo

Aplicar a propriedade de Lyapunov para encontrar um possível controlar que
permita a estabilidade do meu sistema.

A'P + PA < 0
P > 0

Porém, agora meu sistema é do tipo:
dx(t) = Ax(t) + Bu(t)
u(t) = -Kx(t)

QA' + AQ +Y'B' + BY < 0
P > 0

Q = P^-1
Y = -KQ
K = -YQ^-1
%}
clear;
clc;
H = tf(1,[1 0.1 2]); % process
[nH,dH] = tfdata(H,'v');
[A,B,C,D] = tf2ss(nH,dH);
G = ss(A,B,C,D);

Q = sdpvar(length(A));
Y = sdpvar(1,length(B));

F = [Q >= 0, Q*A'+ A*Q + Y'*B' + B*Y <=0];
solvesdp(F)
Q = value(Q);
Y = value(Y);

K = Y*inv(Q)

Amf = A + B*K;
Gmf = ss(Amf,B,C,D);
U = ss(Amf,B,K,0);

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
rlocus(Gmf);






