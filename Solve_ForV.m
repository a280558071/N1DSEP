clear all
load('Case16_Gurobi_BB_V_nox0_withDE2.mat');

v_min=0.95;
v_max=1.05;
M=1e8;
Sbase=1e6;  % unit:VA
Ubase=10e3;  % unit:V
Ibase=Sbase/Ubase/1.732;  %unit: A
Zbase=Ubase/Ibase/1.732;  %unit: ¦¸
z=ConsInf(:,3).*ConsInf(:,5)/Zbase; % line impedance unit:p.u. 
v=sdpvar(N,L,'full');    %Vars for nodal voltage, unit: p.u.
ops=sdpsettings('solver','gurobi');

%% clean s_f1 by making some values zero
s_f1(find(s_f1<=1e-7))=0;
s_y1(find(s_y1<=1e-7))=0;
s_y1(find(abs(s_y1-1)<=1e-7))=1;

%% Cons for calculate Voltage
Cons_Vol=[];
for C_l=1:L
%     Cons_Vol=[Cons_Vol,v_min<=v(:,C_l)<=v_max];
    Cons_Vol=[Cons_Vol,v(N_Subs,C_l)==v_max];
    Cons_Vol=[Cons_Vol,abs(s_f1(:,C_l).*z-Full_I'*v(:,C_l))<=(1-s_y1(:,C_l))*M];
end
sol1=optimize(Cons_Vol,0,ops);
s_v1=value(v);