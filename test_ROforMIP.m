clear all
clc
%% variables statement
x=binvar(2,1,'full');
y=binvar(2,1,'full');
r1=sdpvar(2,1,'full'); % uncertain data
r2=[1,1]';
c=[1,2]';
d=[2,1]';
b=[1,1]';

%% Objective
Obj=c'*x+d'*y;

%% Constraints
Cons=[[r1,r2]*y>=b];
Cons_uncertain=[uncertain(r1),0<=r1<=1];

%% Solve the robust problem
ops=sdpsettings('solver','gurobi');
sol=optimize(Cons+Cons_uncertain,Obj,ops);

%% Save and display the results
s_Obj=value(Obj)
s_y=value(y)
s_x=value(x)