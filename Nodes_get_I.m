function I_Matrix=Nodes_get_I(Nodes_Counts,Line_dat)
    %% �˺���Ŀ��Ϊ��ȡ18�ڵ������滮 �� ���ݣ����ɽڵ�-֧·�������� incidence matrix
    %���룺Nodes_Counts�����ڵ���,Line_dat����18�ڵ������滮���ݣ�Ҫ��Nodes From ��Nodes To ��ǰ��
    %�����A_MatrixΪ��֧��������
    
    %%  ��ȡLine_dat�е�����
    F_Bus=Line_dat(:,1);                  
    T_Bus=Line_dat(:,2);                  %���˵�
    %% ������֧��������A_Matrix������Ϊ��A_Matrix(F_Bus(i),T_Bus(i))=1
    I_Matrix=zeros(Nodes_Counts,length(F_Bus));
    for i=1:length(F_Bus)
        I_Matrix(F_Bus(i),i)=1;                          
        I_Matrix(T_Bus(i),i)=1;
    end
end