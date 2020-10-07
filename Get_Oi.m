function Oi=Get_Oi(i,A)
%% 获得节点i相连的节点j的集合Oi
N=size(A,2);
Oi=[];
for j=1:N
    if A(i,j)==1
        Oi=[Oi,j];
    end
end
end