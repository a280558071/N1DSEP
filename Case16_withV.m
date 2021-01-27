%% On the Topology Constraint for N-1 DSEP
% This function solve the 25-node DSEP problem without considering the voltage
% constraints
% 该函数解决25节点配电网N-1扩展规划问题，未考虑节点电压约束
% [1] Lin Z, Hu Z, Song Y. Distribution Network Expansion Planning Considering N-1 Criterion[J]. 
% IEEE Transactions on Power Systems, 2019, 34(3): 2476-2478.
% [2] Lei S, Chen C, Song Y, et al. Radiality Constraints for Resilient Reconfiguration of Distribution Systems: Formulation and Application to Microgrid Formation
% [J]. IEEE Transactions on Smart Grid, 2020, 11(5): 3944-3956.
% [3] Wang Y, Xu Y, Li J, et al. On the Radiality Constraints for Distribution System Restoration and Reconfiguration Problems[J]. 
% IEEE Transactions on Power Systems, 2020, 35(4): 3294-3296.
% Case Data: a 16-node case with 4 substations from [1]

clear all

st=[
    1 13
    1 2
    1 4
    2 5
    2 14
    14 6
    13 3
    3 4
    3 7
    4 5
    4 8
    5 6
    5 9
    6 10
    7 15
    7 8
    8 11
    8 9
    9 12
    9 10
    10 16
    15 11
    11 12
    12 16
    13 4
    4 9
    9 16
    1 5
    3 8
    8 12
    5 10
    2 6
    7 11
];
CXY=[1 4
2 4
0 3
1 3
2 3
3 3
0 2
1 2
2 2
3 2
1 1
2 1
0 4
3 4
0 1
3 1
]';
s=st(:,1);
t=st(:,2);
L=length(s);
N=16;
Sbase=1e6;  % unit:VA
Ubase=10e3;  % unit:V
Ibase=Sbase/Ubase/1.732;  %unit: A
Zbase=Ubase/Ibase/1.732;  %unit: Ω

N_Subs=[13,14,15,16];
N_Loads=1:12;
G=graph(s,t);
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
labelnode(p,N_Subs,{'13','14','15','16'});
highlight(p,N_Loads,'NodeColor','y','Markersize',20,'NodeFontSize',20);
highlight(p,N_Subs,'Marker','s','NodeColor','c','Markersize',30,'NodeFontSize',40);
highlight(p,s,t,'EdgeColor','k','LineStyle','-.','LineWidth',2,'EdgeFontSize',8);
text(p.XData, p.YData, p.NodeLabel,'HorizontalAlignment', 'center','FontSize', 15); % put nodes' label in right position.
p.NodeLabel={};
% G_Subs = subgraph(G,N_Subs);
% p1=plot(G_Subs);

%% Variable statement
ConsInf=xlsread('case data\16bus33lines.xlsx','G3:L35');
Load=xlsread('case data\16bus33lines.xlsx','D3:D14');  % unit: p.u.
f_Max=ConsInf(:,4)*1e6/Ubase/1.732/Ibase;   % unit: p.u.
z=ConsInf(:,3).*ConsInf(:,5)/Zbase; % line impedance unit:p.u. 
v_min=0.95;
v_max=1.05;
% f_Max=ones(L,1)*50;
g_Sub_Max=50;
M=1e8;
Cost=ConsInf(:,6).*ConsInf(:,3);
x=binvar(L,1,'full');     %Vars for line construction, x(i,1)==1 dentoes that line i is constructed.
y=binvar(L,L,'full');   %Vars for line operation flag in different contigencies, y(line operation,Cont_l)==0
be=binvar(L,L,2,'full');  %Vars for line direction flag in different contigencies
f=sdpvar(L,L,'full');    %Vars for power flow in each line, unit: p.u.
v=sdpvar(N,L,'full');    %Vars for nodal voltage, unit: p.u.
F=sdpvar(L,L,'full');    %Fictitious flow in each line, to complete SCF constraints
D=ones(length(N_Loads),L);   %Fictitious demand in each load
rt=sdpvar(length(N_Loads),L,'full');    %Vars for curtailed load in each load node
g_Sub=sdpvar(length(N_Subs),L,'full');    %Vars for generated power of Subs
Obj=sum(Cost.*x)+M*sum(sum(rt));
% Obj=sum(Cost.*x);

Cons=[];
%% Cons1: Operation logic y<=x in any contigency C_l
Cons_Op=[];
for C_l=1:L
    Cons_Op=[Cons_Op, y(:,C_l)<=x];
    Cons_Op=[Cons_Op, sum(y(:,C_l))==N-length(N_Subs)];
end
Cons=[Cons,Cons_Op];
% size(Cons_Op)
% size(Cons)
%% Cons2: Contigencies happens, when y(i,i)==0, which indicates that outage of line i happens
Cons_Co=[];
for C_l=1:L
    Cons_Co=[Cons_Co, y(C_l,C_l)==0];
end
Cons=[Cons,Cons_Co];
% size(Cons_Co)
% size(Cons)
%% Cons3: Single Commodity Flow Constr.
% Cons_SCF=[];
% M=1e3;
% for i=N_Loads
%     for C_l=1:L
%         li=Full_I(i,:);
%         lij=find(li==1);   % the set of lines start from node i
%         lki=find(li==-1);  % the set of lines end at node i
%         Cons_SCF=[Cons_SCF,sum(F(lij,C_l))+D(i,C_l)==sum(F(lki,C_l))];
%     end
% end
% for C_l=1:L
%     for l=1:L
%         Cons_SCF=[Cons_SCF, abs(F(l,C_l))<=y(l,C_l)*M];
%     end
% end
% Cons=[Cons,Cons_SCF];
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
Cons=[Cons,Cons_De];
%% Cons7: Power balance
Cons_Load=[];
for C_l=1:L
    Cons_Load=[Cons_Load,Full_I*f(:,C_l)==[-(Load-rt(:,C_l));g_Sub(:,C_l)]];
    Cons_Load=[Cons_Load,Load>=rt(:,C_l)>=0];
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

%% Cons10: minimum links —— Jan R H, Hwang F J, Chen S T. Topological optimization of a communication network subject to a reliability constraint[J]. IEEE Transactions on Reliability, 1993, 42(1): 63-70.
% Cons_Links=[sum(x)>=a]; % a is a magical number that limits the reliability of distribution network.
% Cons=[Cons,Cons_Links];

%% Cons11: Voltage limits
Cons_Vol=[];
for C_l=1:L
    Cons_Vol=[Cons_Vol,v_min<=v(:,C_l)<=v_max];
    Cons_Vol=[Cons_Vol,v(N_Subs,C_l)==v_max];
    Cons_Vol=[Cons_Vol,-(1-y(:,C_l))*M<=f(:,C_l).*z-Full_I'*v(:,C_l)<=(1-y(:,C_l))*M];
end
Cons=[Cons,Cons_Vol];

%% Set initial guess of x,y and be_Nodes to values in "Case25_noV_nox0_withDE3_realf12.mat"
ops=sdpsettings('solver','cplex','verbose',2,'cplex.mip.limits.cutpasses',-1,'cplex.mip.display',3,'cplex.mip.strategy.heuristicfreq',-1,'usex0',1); %,'cplex.mip.tolerances.mipgap',5e-2);
% ops=sdpsettings('solver','gurobi', 'gurobi.Heuristics',0,'gurobi.Cuts',0,'usex0',0);%'gurobi.MIPGap',5e-2,
load('Case16_Cplex_BB_V_nox0_withDE3.mat','s_x1','s_y1','s_be1');
assign(x,s_x1);
assign(y,s_y1);
assign(be,s_be1);

%% solve the problem
sol1=optimize(Cons,Obj,ops);
%% Save the solution with "s_" as start
s_x1=value(x);
s_y1=value(y);
s_be1=value(be);
s_f1=value(f);
s_v1=value(v);
s_rt1=value(rt);
s_g_Sub1=value(g_Sub);
s_Obj1=value(Obj);
save('Case16_Cplex_BB_V_withx0_withDE2');
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
    Gi=graph(s_temp,t_temp);
    pi(i)=plot(Gi,'Layout','force');
    pi(i).XData=CXY(1,:);
    pi(i).YData=CXY(2,:);
    labelnode(pi(i),N_Subs,{'13','14','15','16'});
    highlight(pi(i),N_Loads,'NodeColor','y','Markersize',20,'NodeFontSize',20);
    highlight(pi(i),N_Subs,'Marker','s','NodeColor','c','Markersize',30,'NodeFontSize',40);
    highlight(pi(i),s_temp(find(round(s_x1)==0)),t_temp(find(round(s_x1)==0)),'LineStyle','-.','LineWidth',1,'EdgeColor','k');
    highlight(pi(i),s_temp(find(round(s_x1)==1)),t_temp(find(round(s_x1)==1)),'EdgeColor','g','LineWidth',2); % bold line denotes the newly-built line
    highlight(pi(i),s_temp(find(round(s_y1(:,i))==1)),t_temp(find(round(s_y1(:,i))==1)),'EdgeColor','b','LineWidth',6,'LineStyle','-');  % blue lines are lines in operation
    l_Open=find(round((1-s_y1(:,i)).*s_x1)==1); 
    highlight(pi(i),s_temp(l_Open),t_temp(l_Open),'EdgeColor','g','LineWidth',2,'LineStyle','-'); % green line denotes the newly-built line NOT in operation
    highlight(pi(i),s_temp(i),t_temp(i),'EdgeColor','r','LineStyle','-','LineWidth',4);  % red line denotes the outage line
    text(pi(i).XData, pi(i).YData, pi(i).NodeLabel,'HorizontalAlignment', 'center','FontSize', 15); % put nodes' label in right position.
    text(0.5*(pi(i).XData(s_temp(i))+pi(i).XData(t_temp(i))), 0.5*(pi(i).YData(s_temp(i))+pi(i).YData(t_temp(i))), 'X','HorizontalAlignment', 'center','FontSize', 15,'Color','r'); % label outage lines with 'x' in the middle
    l_newOpen=setdiff(l_Open,i);
    text(0.5*(pi(i).XData(s_temp(l_newOpen))+pi(i).XData(t_temp(l_newOpen))), 0.5*(pi(i).YData(s_temp(l_newOpen))+pi(i).YData(t_temp(l_newOpen))), 'O','HorizontalAlignment', 'center','FontSize', 15,'Color','g'); % label newly-built and opened lines with 'o' in the middle
    pi(i).NodeLabel={};
end

%% find where rt!=0 and plot the Contigency senario
% [xr,yr]=find(s_rt1~=0);
% j=0;
% for i=yr'
%     j=j+1;
%     figure;
%     s_temp=s;
%     t_temp=t;
%     N_RLoads=find(N_Loads~=xr(j));
%     for l=1:L
%         if s_be1(l,i,2)==1
%             s1=s_temp(l);
%             s_temp(l)=t(l);
%             t_temp(l)=s1;
%         end
%     end
%     Gi=graph(s_temp,t_temp);
%     pir(j)=plot(Gi,'Layout','force');
%     pir(j).XData=CXY(1,:);
%     pir(j).YData=CXY(2,:);
%     labelnode(pir(j),N_Subs,{'13','14','15','16'});
%     highlight(pir(j),N_Loads,'NodeColor','y','Markersize',20,'NodeFontSize',20);
%     highlight(pir(j),N_Subs,'Marker','s','NodeColor','c','Markersize',30,'NodeFontSize',40);
%     highlight(pir(j),s_temp(find(s_x1==0)),t_temp(find(s_x1==0)),'LineStyle','--','LineWidth',1);
%     highlight(pir(j),s_temp(find(s_x1==1)),t_temp(find(s_x1==1)),'LineWidth',4); % bold line denotes the newly-built line
%     highlight(pir(j),s_temp(i),t_temp(i),'EdgeColor','k','LineStyle','-','LineWidth',4);  % black line denotes the outage line
%     highlight(pir(j),s_temp(find(s_y1(:,i)==1)),t_temp(find(s_y1(:,i)==1)),'EdgeColor','r','LineWidth',6,'LineStyle','-');
%     labelnode(pir(j),xr(j),{[num2str(Load(xr(j))) ' Lost ' num2str(s_rt1(xr(j),yr(j)))]});
%     for k=1:length(N_RLoads)
%         labelnode(pir(j),N_RLoads(k),num2str(Load(N_RLoads(k))));
%     end
%     title(['Contigency ' num2str(i) ' happens, shedding load in node ', num2str(xr(j))]);
%     text(pir(i).XData, pir(i).YData, pir(i).NodeLabel,'HorizontalAlignment', 'center','FontSize', 15); % put nodes' label in right position.
%     pir(i).NodeLabel={};
% end
