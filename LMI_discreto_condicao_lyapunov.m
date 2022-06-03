%{
Código: LMI Caso Discreto

Neste código está sendo aplicado a condição de Lyapunov para identificar a
estabilidade de um sistema.

A'PA - P < 0
P > 0

Com o sistema estável, espera-se que se encontre um P que torne a condição
verdadeira.

Com o sistema instável, espera-se que não se encontre um P.
%}
Ts = 0.1; %tempo de amostragem

A = [1 2 0;3 4 1;0 0 2]; %sistema instável
%A = [-1 2 0;-3 -4 1;0 0 -2]; %sistema estável
B = [1;0;0];
C = [1;0;0]';
D = 0;

G = ss(A,B,C,D); %contínuo
step(G);

G = c2d(G,Ts); %discreto

[A,B,C,D] = ssdata(G);%obtendo as novas matrizes para o caso discreto

P = sdpvar(length(A));
F = [P>=0, A'*P*A - P <= 0];

sol = solvesdp(F);
P = value(P);

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