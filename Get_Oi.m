function Oi=Get_Oi(i,A)
%% ��ýڵ�i�����Ľڵ�j�ļ���Oi
N=size(A,2);
Oi=[];
for j=1:N
    if A(i,j)==1
        Oi=[Oi,j];
    end
end
end