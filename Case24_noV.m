%% On the Topology Constraint for N-1 DSEP
% This function solve the 24-node DSEP problem without considering the voltage
% constraints
% 该函数解决24节点配电网N-1扩展规划问题，未考虑节点电压约束
% [1] Lin Z, Hu Z, Song Y. Distribution Network Expansion Planning Considering N-1 Criterion[J]. 
% IEEE Transactions on Power Systems, 2019, 34(3): 2476-2478.
% [2] Lei S, Chen C, Song Y, et al. Radiality Constraints for Resilient Reconfiguration of Distribution Systems: Formulation and Application to Microgrid Formation
% [J]. IEEE Transactions on Smart Grid, 2020, 11(5): 3944-3956.
% [3] Wang Y, Xu Y, Li J, et al. On the Radiality Constraints for Distribution System Restoration and Reconfiguration Problems[J]. 
% IEEE Transactions on Power Systems, 2020, 35(4): 3294-3296.
% Case Data: a 24-node case with 4 substations from [1]

clear all

st=[2	21
6	22
1	21
8	22
1	5
1	9
1	14
2	3
2	12
3	10
3	16
3	23
4	7
4	9
4	15
% 4	16
5	6
5	24
6	13
6	17
7	8
7	11
7	19
7	23
10	16
10	23
11	23
14	18
15	17
15	19
17	22
18	24
20	24
13  24
20  18
12  21
12  10
9   21
9   5
2   16
3   4
];
CXY=[1.31	0.11	-0.96	-0.64	0.90	0.00	-1.60	-1.69	0.55	-0.63	-1.98	0.29	0.80	2.14	-0.10	-0.60	-0.75	2.35	-1.17	2.35	0.86	-1.17	-1.28	1.61
0.07	-0.70	-1.02	-0.44	0.91	1.53	-0.59	0.58	-0.11	-1.99	-1.34	-1.75	2.05	0.70	0.53	-1.41	1.16	1.54	0.07	2.14	-0.97	2.05	-1.61	1.75];
s=st(:,1);
t=st(:,2);
L=length(s);
N=24;

N_Subs=[21,22,23,24];
N_Loads=1:20;
G=digraph(s,t);
idxOut = findedge(G,s,t); % !!!!!!*************** index of edge is not the same with that of in mpc.branch !!!!!!*******
A=adjacency(G);  % 邻接矩阵
Full_A=full(A);
I=incidence(G);  % 节支关联矩阵
fI=abs(I);
Full_I=full(I(:,idxOut));

%% Plot the topology of DN
figure;
p=plot(G,'Layout','force');
p.XData=CXY(1,:);
p.YData=CXY(2,:);
hold on;
labelnode(p,N_Subs,{'Sub_2_1','Sub_2_2','Sub_2_3','Sub_2_4'});
highlight(p,N_Subs,'Marker','s','NodeColor','g');
% G_Subs = subgraph(G,N_Subs);
% p1=plot(G_Subs);

%% Variable statement
ConsInf=xlsread('multistage_node24.xlsx','G3:R42');
Load=xlsread('multistage_node24.xlsx','D2:D21');
f_Max=ConsInf(:,10);
f_Max=ones(L,1)*50;
g_Sub_Max=50;
VOLL=1e8;
Cost=ConsInf(:,12).*ConsInf(:,3);
x=binvar(L,1,'full');     %Vars for line construction, x(i,1)==1 dentoes that line i is constructed.
y=binvar(L,L,'full');   %Vars for line operation flag in different contigencies, y(line operation,Cont_l)==0
be=binvar(L,L,2,'full');  %Vars for line direction flag in different contigencies
f=sdpvar(L,L,'full');    %Vars for power flow in each line
rt=sdpvar(length(N_Loads),L,'full');    %Vars for curtailed load in each load node
g_Sub=sdpvar(length(N_Subs),L,'full');    %Vars for generated power of Subs
Obj=sum(Cost.*x)+VOLL*sum(sum(rt));
% Obj=sum(Cost.*x);

Cons=[];
%% Cons1: Operation logic y<=x in any contigency C_l
Cons_Op=[];
for C_l=1:L
    Cons_Op=[Cons_Op, y(:,C_l)<=x];
    Cons_Op=[Cons_Op, sum(y(:,C_l))==20];
end
Cons=[Cons,Cons_Op];
% size(Cons_Op)
% size(Cons)
%% Cons2: Contigencies happens, when y(i,i)==0, which indicates that outage of line i happens
Cons_Co=[];
for C_l=1:L
    Cons_Co=[Cons_Co, y(C_l,C_l)==0];
    Cons_Op=[Cons_Op, sum(y(:,C_l))==20];
end
Cons=[Cons,Cons_Co];
% size(Cons_Co)
% size(Cons)
%% Cons3: Single Commodity Flow Constr.
%% Cons4: Fencing Constr.
%% Cons5: Spanning Tree Constr. 
% J. A. Taylor and F. S. Hover, “Convex models of distribution system reconfiguration,” IEEE Trans. Power Syst., vol. 27, no. 3, pp. 1407–1413,Aug. 2012.
Cons_ST=[];
for C_l=1:L
    for l=1:L
        Os=Full_I(s(l),:);
        Os1=find(Os==1);   % the set of lines end at node s(l)
        Os2=find(Os==-1);  % the set of lines start from node s(l)
        Ot=Full_I(t(l),:);
        Ot1=find(Ot==1);   % the set of lines end at node t(l)
        Ot2=find(Ot==-1);   % the set of lines start from node t(l)
        Cons_ST=[Cons_ST,be(l,C_l,1)+be(l,C_l,2)==y(l,C_l)];
        if ismember(s(l),N_Subs)
            Cons_ST=[Cons_ST,be(l,C_l,2)==0];  % be(l,C_l,2)==0 means node t(l) cannot be the parent node of s(l)
        end
        if ismember(t(l),N_Subs)
            Cons_ST=[Cons_ST,be(l,C_l,1)==0];  % be(l,C_l,1)==0 means node s(l) cannot be the parent node of t(l)
        end
        if ismember(s(l),N_Loads)
            Cons_ST=[Cons_ST,sum([be(Os1,C_l,1);be(Os2,C_l,2)])==1];  %sum(be(Os2,C_l,2))==1 means only one node start from s(l) can be the parent node of s(l)
        end
        if ismember(t(l),N_Loads)
            Cons_ST=[Cons_ST,sum([be(Ot1,C_l,1);be(Ot2,C_l,2)])==1];  %sum(be(Ot,C_l,1)==1 means only one node in Ot can be the parent node of t(l)
        end
    end
end
Cons=[Cons,Cons_ST];
% size(Cons_ST)
% size(Cons)
%% Cons6: Degree of Each Node Constr.
Cons_De=[];
for i=N_Loads
    Cons_De=[Cons_De,sum(x([find(s==i);find(t==i)]))>=2];
end
% Cons=[Cons,Cons_De];
% size(Cons_De)
% size(Cons)
%% Cons7: Power balance
Cons_Load=[];
for C_l=1:L
    Cons_Load=[Cons_Load,Full_I*f(:,C_l)==[-(Load-rt(:,C_l));g_Sub(:,C_l)]];
    Cons_Load=[Cons_Load,Load>=rt(:,C_l)>=0,rt(:,C_l)==0];
end
Cons=[Cons,Cons_Load];
% size(Cons_Load)
% size(Cons)
%% Cons8: power flow limitation in each line
Cons_Line=[];
for C_l=1:L
    Cons_Line=[Cons_Line,-y(:,C_l).*f_Max<=f(:,C_l)<=y(:,C_l).*f_Max];
end
Cons=[Cons,Cons_Line];
% size(Cons_Line)
% size(Cons)

%% Cons9: Power limits of Subs
Cons_Sub=[0<=g_Sub<=g_Sub_Max];
Cons=[Cons,Cons_Sub];
% size(Cons_Sub)
% size(Cons)

%% Set initial guess of x,y and be_Nodes to values in "Case24_nof.mat"
ops=sdpsettings('solver','cplex','verbose',2,'cplex.mip.display',3,'usex0',0);
% load('Case24_nof.mat','s_x','s_y','s_be');
% assign(x,s_x);
% assign(y,s_y);
% assign(be,s_be);

%% solve the problem
sol1=optimize(Cons,Obj,ops);
%% Save the solution with "s_" as start
s_x1=value(x);
s_y1=value(y);
s_be1=value(be);
s_f1=value(f);
s_rt1=value(rt);
s_g_Sub1=value(g_Sub);
s_Obj1=value(Obj);
save('Case24_noV_nox0');
%% Highlight the lines to be bulit and plot all the operation conditions
for i=1:5 % Contigency i happens
    figure;
    s_temp=s;
    t_temp=t;
    for l=1:L
        if s_be1(l,i,2)==1
            s1=s_temp(l);
            s_temp(l)=t(l);
            t_temp(l)=s1;
        end
    end
    Gi=digraph(s_temp,t_temp);
    pi(i)=plot(Gi,'Layout','force');
    pi(i).XData=CXY(1,:);
    pi(i).YData=CXY(2,:);
    labelnode(pi(i),N_Subs,{'Sub_2_1','Sub_2_2','Sub_2_3','Sub_2_4'});
    highlight(pi(i),N_Subs,'Marker','s','NodeColor','g','MarkerSize',10);
    highlight(pi(i),s_temp(find(s_x1==0)),t_temp(find(s_x1==0)),'LineStyle','--','LineWidth',1);
    highlight(pi(i),s_temp(find(s_x1==1)),t_temp(find(s_x1==1)),'LineWidth',4); % bold line denotes the newly-built line
    highlight(pi(i),s_temp(i),t_temp(i),'EdgeColor','k','LineStyle','-','LineWidth',4);  % black line denotes the outage line
    highlight(pi(i),s_temp(find(s_y1(:,i)==1)),t_temp(find(s_y1(:,i)==1)),'EdgeColor','r','LineWidth',6,'LineStyle','-');
end

% figure;
% pi=plot(G,'Layout','force');
% labelnode(pi,N_Subs,{'Sub_2_1','Sub_2_2','Sub_2_3','Sub_2_4'});
% highlight(pi,N_Subs,'Marker','s','NodeColor','g','MarkerSize',10);
% highlight(pi,s(find(s_x==1)),t(find(s_x==1)),'LineWidth',4); % bold line denotes the newly-built line
% 
% figure;
% pi1=plot(G,'Layout','force');
% labelnode(pi1,N_Subs,{'Sub_2_1','Sub_2_2','Sub_2_3','Sub_2_4'});
% highlight(pi1,N_Subs,'Marker','s','NodeColor','g','MarkerSize',10);
% highlight(pi1,s(find(s_x1==1)),t(find(s_x1==1)),'LineWidth',4); % bold line denotes the newly-built line
