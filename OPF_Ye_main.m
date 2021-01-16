clc;clear;

%load case9 data
%[baseMVA, bus, gen, branch, areas, gencost, info] = loadcase(WB2);u = 10; Vmin = [0.95; 0.95]; Vmax = [1.05; 1.009]; % 1-step
%[baseMVA, bus, gen, branch, areas, gencost, info] = loadcase(WB2);u = 10; Vmin = [0.95; 0.95]; Vmax = [1.05; 0.983]; % 2-step
%[baseMVA, bus, gen, branch, areas, gencost, info] = loadcase(WB3);u = 1; Vmin = 0.95; Vmax = 1.05; % 1-step
[baseMVA, bus, gen, branch, areas, gencost, info] = loadcase(WB5);u = 1.40; Vmin = 0.95; Vmax = 1.05; % 7-step
%[baseMVA, bus, gen, branch, areas, gencost, info] = loadcase(LMBM3);u = 1; Vmin = 0.90; Vmax = 1.10; % 1-step
%[baseMVA, bus, gen, branch, areas, gencost, info] = loadcase(case9); u = 1; Vmin = 0.95; Vmax = 1.05; % 2-step;
%[baseMVA, bus, gen, branch, areas, gencost, info] = loadcase(case9mod); u = 8000; Vmin = 0.95; Vmax = 1.05; % failed;
%[baseMVA, bus, gen, branch, areas, gencost, info] = loadcase(case22loop);u = 100; Vmin = 0.95; Vmax = 1.05; % 1-step
%[baseMVA, bus, gen, branch, areas, gencost, info] = loadcase(case39mod1);u = 4500; Vmin = 0.95; Vmax = 1.05; % 661-step;
%[baseMVA, bus, gen, branch, areas, gencost, info] = loadcase(case39mod2);u = 100; Vmin = 0.95; Vmax = 1.05;
%[baseMVA, bus, gen, branch, areas, gencost, info] = loadcase(case118mod);u = 100; Vmin = 0.95; Vmax = 1.05;
%[baseMVA, bus, gen, branch, areas, gencost, info] = loadcase(case300mod);u = 100; Vmin = 0.95; Vmax = 1.05;
%[baseMVA, bus, gen, branch, areas, gencost, info] = loadcase(case30); u = 1; Vmin = 0.95; Vmax = 1.05;
%[baseMVA, bus, gen, branch, areas, gencost, info] = loadcase(case14); u = 1; Vmin = 0.95; Vmax = 1.05;
%[baseMVA, bus, gen, branch, areas, gencost, info] = loadcase(case57); u = 1; Vmin = 0.95; Vmax = 1.05;
% [baseMVA, bus, gen, branch, areas, gencost, info] = loadcase(case30mod); u = 1; Vmin = 0.95; Vmax = 1.05;
% [baseMVA, bus, gen, branch, areas, gencost, info] = loadcase(case14mod); u = 1; Vmin = 0.95; Vmax = 1.05;
% [baseMVA, bus, gen, branch, areas, gencost, info] = loadcase(case57mod); u = 1; Vmin = 0.95; Vmax = 1.05;

% admitance 
[Y, Yf, Yt] = makeYbus(baseMVA, bus, branch);

n_bus = length(bus(:,1)); 
n_gen = length(gen(:,1));
n_branch = length(branch(:,1));
Pg = sdpvar(n_gen,1);
Qg = sdpvar(n_gen,1);
W = sdpvar(n_bus,n_bus,'hermitian','complex');

% define the load bus
if (n_bus - n_gen) > 0
    c = 0;
for k = 1:n_bus
    flag = 0;
    for j = 1:n_gen
        if bus(k,1) == gen(j,1)
            flag = flag + 1;
        end
    end
    if flag ==0
         c = c+1; loadbus(c) = bus(k);
    end
end
n_load = length(loadbus);
end

Pl = 0.01*bus(:,3); Ql = 0.01*bus(:,4);
Pmin = 0.01*gen(:,10); Pmax = 0.01*gen(:,9);
Qmin = 0.01*gen(:,5); Qmax = 0.01*gen(:,4);

F = [Pg >= Pmin , Pg <= Pmax];
F = [F , Qg >= Qmin , Qg <= Qmax];


%for generator bus 
for k = 1:n_gen
    F = [F, Pg(k)-Pl(gen(k,1)) + (Qg(k)-Ql(gen(k,1)))*1i ...
        == sum(W(gen(k,1),:).*conj(Y(gen(k,1),:)))];
end

%for load bus 
if (n_bus - n_gen) > 0
for k = 1:n_load
    F = [F, -Pl(loadbus(k))-Ql(loadbus(k))*1i ...
        == sum(W(loadbus(k),:).*conj(Y(loadbus(k),:)))];
end
end

F = [F , diag(W) >= Vmin.^2 , diag(W) <= Vmax.^2];

% for k = 1:n
%     for m = 1:n
%         F = [F , imag(W(k,m)) <= real(W(k,m))*tan(45*pi/180)];
%     end
% end

F = [F , W >= 0];

obj = sum(gencost(:,5).*(100*Pg).^2 + gencost(:,6).*(100*Pg) + gencost(:,7));

option = sdpsettings('solver','sedumi');
result = solvesdp(F,obj,option),
%optimize(F,obj);
%optimize(F,[],sdpsettings('lmirank.solver','sedumi','sedumi.eps',0))
Wf = value(W); P0 = value(Pg); Q0 = value(Qg); W0 = Wf; lamda0 = eig(W0),
[v0,D0] = eig(Wf);
optimal_initial = value(obj)

[lamda,index] = max(max(D0)); %max eigenvalue
x = v0(:,index); %max eigenvector

step = 1;

while trace(Wf) - lamda >  1e-5
    obj = abs(trace(W)-x'*W*x);
    option = sdpsettings('solver','sedumi','verbose',0);
    result = solvesdp(F,obj,option)
    Wf = value(W);
    Pf = value(Pg);
    Qf = value(Qg);
    [v,D] = eig(Wf);
    [lamda,index] = max(max(D)); %max eigenvalue
    x = v(:,index); %max eigenvector
    step = step + 1,
    value(obj); eig(Wf)
    optimal = sum(gencost(:,5).*(100*Pf).^2 + gencost(:,6).*(100*Pf) + gencost(:,7))
    check3 = trace(Wf) - lamda;
end








