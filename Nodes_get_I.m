function I_Matrix=Nodes_get_I(Nodes_Counts,Line_dat)
    %% 此函数目的为读取18节点配网规划 的 数据，生成节点-支路关联矩阵 incidence matrix
    %输入：Nodes_Counts——节点数,Line_dat——18节点配网规划数据，要求Nodes From 和Nodes To 在前列
    %输出：A_Matrix为节支关联矩阵
    
    %%  读取Line_dat中的数据
    F_Bus=Line_dat(:,1);                  
    T_Bus=Line_dat(:,2);                  %两端点
    %% 建立节支关联矩阵A_Matrix，规则为：A_Matrix(F_Bus(i),T_Bus(i))=1
    I_Matrix=zeros(Nodes_Counts,length(F_Bus));
    for i=1:length(F_Bus)
        I_Matrix(F_Bus(i),i)=1;                          
        I_Matrix(T_Bus(i),i)=1;
    end
end