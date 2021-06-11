%% DNR for NS
% This function solve a distribution network reconfiguration problem in
% Nansha, Guangzhou, China, in which 4 distribution lines are included and
% connected with each other, where the power is supplied by two
% substations.

clear all

%% read data from xlsx file
ConsInf=xlsread('case data\NS4lines.xlsx','G3:L65');
Load=xlsread('case data\NS4lines.xlsx','D3:D49');
L_rate=0.3;         %load rate 30% (负载率)
Load=Load*L_rate;        
s=ConsInf(:,1);
t=ConsInf(:,2);
w=ConsInf(:,3);
L_lines=47;           %number of distribution lines
L_s=16;         %number of switches
N=49;
Sbase=1e6;      % unit:VA
Ubase=10e3;     % unit:V
Ibase=Sbase/Ubase/1.732;  %unit: A
Zbase=Ubase/Ibase/1.732;  %unit: Ω

N_Subs=[48,49];
N_Loads=1:47;
G=digraph(s,t,w);           % !!!!!!*************** use directed graph is very,very important for index of starting nodes and ending nodes, otherwise Spanning Tree Constraints would not work.
idxOut = findedge(G,s,t); % !!!!!!*************** index of edge is not the same with that of in mpc.branch !!!!!!*******
% [s,t]=findedge(G);
A=adjacency(G);  % 邻接矩阵
Full_A=full(A);
I=incidence(G);  % 节支关联矩阵
fI=abs(I);
Full_I=full(I(:,idxOut));

%% Plot the topology of DN
% figure;
% p=plot(G,'Layout','layered','Direction','right');
% % view(-90,90);
% % p.XData=CXY(:,1);
% % p.YData=CXY(:,2);
% hold on;
% labelnode(p,N_Subs,{'横一','横沥'});
% highlight(p,N_Loads,'NodeColor','y','Markersize',20,'NodeFontSize',20);
% highlight(p,N_Subs,'Marker','s','NodeColor','c','Markersize',30,'NodeFontSize',40);
% highlight(p,s,t,'EdgeColor','k','LineStyle','-.','LineWidth',2,'EdgeFontSize',8);
% text(p.XData, p.YData, p.NodeLabel,'HorizontalAlignment', 'center','FontSize', 15); % put nodes' label in right position.
% p.NodeLabel={};
% G_Subs = subgraph(G,N_Subs);
% p1=plot(G_Subs);

%% Variable statement
f_Max=ConsInf(:,4)*1e6/Ubase/1.732/Ibase*10;   % unit: p.u.
z=w.*ConsInf(:,5)/Zbase; % line impedance unit:p.u. 
r=w.*ConsInf(:,6)/Zbase; % line resistance unit:p.u. 
v_min=0.9;
v_max=1.1;
% f_Max=ones(L,1)*50;
g_Sub_Max=100;
M=1e7;
M2=1e3;
y=binvar(L_lines+L_s,1,'full');   %Vars for line operation flag in different contigencies, y(line operation,Cont_l)==0
be=binvar(L_lines+L_s,2,'full');  %Vars for line direction flag in different contigencies
f=sdpvar(L_lines+L_s,1,'full');    %Vars for power flow in each line，unit: p.u.
v=sdpvar(N,1,'full');    %Vars for nodal voltage, unit: p.u.
% F=sdpvar(L,L,'full');    %Fictitious flow in each line, to complete SCF constraints
% D=ones(length(N_Loads),L);   %Fictitious demand in each load
rt=sdpvar(length(N_Loads),1,'full');    %Vars for curtailed load in each load node
g_Sub=sdpvar(length(N_Subs),1,'full');    %Vars for generated power of Subs
% Obj=sum(f.^2.*r)+M*sum(sum(rt));    %the objective is to minimize the power loss
Obj=sum(f.^2.*r);    %the objective is to minimize the power loss

Cons=[];
%% Cons1: Number of lines must be equal to number of load_nodes (actually it is included if ST constraints is applied)
% Cons_Op=[];
% Cons_Op=[Cons_Op, sum(y)==N-length(N_Subs)];
% Cons=[Cons,Cons_Op];
% size(Cons_Op)
% size(Cons)
%% Cons2: Connections within Switch station should be 0 in Normal Condition, i.e. y(i)==0 for i∈L_s, which indicates i is a switch
Cons_Co_normal=[y(L_lines+1:end)==0];
Cons=[Cons,Cons_Co_normal];
% size(Cons_Co)
% size(Cons)
%% Cons3: Substations must be connected with load nodes
% subs_lines1=[1,9];
% subs_lines2=[19,29];
% Cons_Co_subs=[sum(y(subs_lines1))==1,sum(y(subs_lines2))==1];
% Cons=[Cons,Cons_Co_subs];
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
%% Cons4: Spanning Tree Constr. 
% J. A. Taylor and F. S. Hover, “Convex models of distribution system reconfiguration,” IEEE Trans. Power Syst., vol. 27, no. 3, pp. 1407–1413,Aug. 2012.
Cons_ST=[];
for l=1:(L_lines+L_s)
    Os=Full_I(s(l),:);
    Os1=find(Os==1);   % the set of lines end at node s(l)
    Os2=find(Os==-1);  % the set of lines start from node s(l)
    Ot=Full_I(t(l),:);
    Ot1=find(Ot==1);   % the set of lines end at node t(l)
    Ot2=find(Ot==-1);   % the set of lines start from node t(l)
    Cons_ST=[Cons_ST,be(l,1)+be(l,2)==y(l)];
    if ismember(s(l),N_Subs)
        Cons_ST=[Cons_ST,be(l,2)==0];  % be(l,2)==0 means node t(l) cannot be the parent node of s(l)
    end
    if ismember(t(l),N_Subs)
        Cons_ST=[Cons_ST,be(l,1)==0];  % be(l,1)==0 means node s(l) cannot be the parent node of t(l)
    end
    if ismember(s(l),N_Loads)
        Cons_ST=[Cons_ST,sum([be(Os1,1);be(Os2,2)])==1];  %sum(be(Os,2))==1 means only one node in Os can be the parent node of s(l)
    end
    if ismember(t(l),N_Loads)
        Cons_ST=[Cons_ST,sum([be(Ot1,1);be(Ot2,2)])==1];  %sum(be(Ot,1)==1 means only one node in Ot can be the parent node of t(l)
    end
end
Cons=[Cons,Cons_ST];
% size(Cons_ST)
% size(Cons)
%% Cons6: Degree of Substation nodes.
Cons_De=[];
for i=N_Subs
    Cons_De=[Cons_De,sum(y([find(s==i);find(t==i)]))>=1];
end
N_TransferNode=[18,20,24,28];
for i=N_TransferNode
    Cons_De=[Cons_De,sum(y([find(s==i);find(t==i)]))>=1];
end
Cons=[Cons,Cons_De];
% size(Cons_De)
% size(Cons)
%% Cons7: Power balance
Cons_Load=[];
% Cons_Load=[Cons_Load,Full_I*f==[-(Load-rt);g_Sub]];
Cons_Load=[Cons_Load,Full_I*f==[-(Load);g_Sub]];
% Cons_Load=[Cons_Load,Load>=rt>=0];
Cons=[Cons,Cons_Load];
% size(Cons_Load)
% size(Cons)
%% Cons8: power flow limitation in each line
Cons_Line=[];
Cons_Line=[Cons_Line,-y.*f_Max<=f<=y.*f_Max];
Cons=[Cons,Cons_Line];
% size(Cons_Line)
% size(Cons)

%% Cons9: Power limits of Subs
Cons_Sub=[0<=g_Sub<=g_Sub_Max];
% Cons=[Cons,Cons_Sub];
% size(Cons_Sub)
% size(Cons)

%% Cons10: minimum links —— Jan R H, Hwang F J, Chen S T. Topological optimization of a communication network subject to a reliability constraint[J]. IEEE Transactions on Reliability, 1993, 42(1): 63-70.
% Cons_Links=[sum(x)>=a]; % a is a magical number that limits the reliability of distribution network.
% Cons=[Cons,Cons_Links];
%% Cons11: Voltage limits
Cons_Vol=[];
Cons_Vol=[Cons_Vol,v_min<=v<=v_max];
Cons_Vol=[Cons_Vol,v(N_Subs)==v_max];   
Cons_Vol=[Cons_Vol,-(1-y)*M2<=f.*z-Full_I'*v<=(1-y)*M2];  %if y==1, then f.*z-Full_I'*v==0; if y==1, then -1000<=f.*z-Full_I'*v<=1000
Cons=[Cons,Cons_Vol];
%% Set optimization settings
% ops=sdpsettings('solver','cplex','verbose',2,'cplex.mip.display',3,'usex0',1,'cplex.mip.tolerances.mipgap',5e-2);%,'cplex.mip.limits.cutpasses',-1,'cplex.mip.tolerances.integrality',1e-8);
ops=sdpsettings('solver','gurobi','gurobi.MIPGap',1e-4);% ,'gurobi.Heuristics',0,'gurobi.Cuts',0

%% solve the problem
sol1=optimize(Cons,Obj,ops);
%% Save the solution with "s_" as start
s_y1=value(y);
s_be1=value(be);
s_f1=value(f);
% s_rt1=value(rt);
s_v1=value(v);
s_g_Sub1=value(g_Sub);
s_Obj1=value(Obj);

display(['当前配网网损为：  ',num2str(s_Obj1),' MW']);
% display(['最低节点电压： ', num2str(min(s_v1)), ' p.u.']);
% display(['最高线路负载率： ', num2str(max(abs(s_f1)./f_Max)*100), ' %']);
%% Highlight the lines to be bulit and plot all the operation conditions
figure;
s_temp=s;
t_temp=t;
for l=1:(L_lines+L_s)
    if s_be1(l,2)==1
        s1=s_temp(l);
        s_temp(l)=t(l);
        t_temp(l)=s1;
    end
end
% Gi=graph(s_temp,t_temp);
pi=plot(G,'Layout','layered','Direction','right');
% view(-90,90)
% figure
% pi_test=plot(Gi,'Layout','subspace3');
%     pi(i).XData=CXY(1,:);
%     pi(i).YData=CXY(2,:);
labelnode(pi,N_Subs,{'横一','横沥'});
highlight(pi,N_Loads,'NodeColor','y','Markersize',20,'NodeFontSize',20);
highlight(pi,N_Subs,'Marker','s','NodeColor','c','Markersize',30,'NodeFontSize',40);
%     highlight(pi(i),s_temp(find(s_x1==0)),t_temp(find(s_x1==0)),'LineStyle','-.','LineWidth',1,'EdgeColor','k');
%     highlight(pi(i),s_temp(find(s_x1==1)),t_temp(find(s_x1==1)),'EdgeColor','g','LineWidth',2); % bold line denotes the newly-built line
highlight(pi,s(find(s_y1==1)),t(find(s_y1==1)),'EdgeColor','b','LineWidth',6,'LineStyle','-');  % blue lines are lines in operation
%     l_Open=find((1-s_y1(:,i)).*s_x1==1);
%     highlight(pi(i),s_temp(l_Open),t_temp(l_Open),'EdgeColor','g','LineWidth',2,'LineStyle','-'); % green line denotes the newly-built line NOT in operation
% highlight(pi,s_temp,t_temp,'EdgeColor','r','LineStyle','-','LineWidth',4);  % red line denotes the outage line
text(pi.XData, pi.YData, pi.NodeLabel,'HorizontalAlignment', 'center','FontSize', 15); % put nodes' label in right position.
% text(0.5*(pi.XData(s_temp)+pi.XData(t_temp)), 0.5*(pi.YData(s_temp)+pi.YData(t_temp)), 'X','HorizontalAlignment', 'center','FontSize', 15,'Color','r'); % label outage lines with 'x' in the middle
%     l_newOpen=setdiff(l_Open,i);
% text(0.5*(pi.XData(s_temp(l_newOpen))+pi.XData(t_temp(l_newOpen))), 0.5*(pi.YData(s_temp(l_newOpen))+pi.YData(t_temp(l_newOpen))), 'O','HorizontalAlignment', 'center','FontSize', 15,'Color','g'); % label newly-built and opened lines with 'o' in the middle
pi.NodeLabel={};

%% Save all the data
% save('Case_DNR_NS_NormalCondition');
