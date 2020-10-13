[xr,yr]=find(s_rt1~=0);
j=0;
for i=yr'
    j=j+1;
    figure;
    s_temp=s;
    t_temp=t;
    N_RLoads=find(N_Loads~=xr(j));
    for l=1:L
        if s_be1(l,i,2)==1
            s1=s_temp(l);
            s_temp(l)=t(l);
            t_temp(l)=s1;
        end
    end
    Gi=digraph(s_temp,t_temp);
    pir(j)=plot(Gi,'Layout','force');
    pir(j).XData=CXY(1,:);
    pir(j).YData=CXY(2,:);
    labelnode(pir(j),N_Subs,{'Sub_2_1','Sub_2_2','Sub_2_3','Sub_2_4'});
    highlight(pir(j),N_Subs,'Marker','s','NodeColor','g','MarkerSize',10);
    highlight(pir(j),s_temp(find(s_x1==0)),t_temp(find(s_x1==0)),'LineStyle','--','LineWidth',1);
    highlight(pir(j),s_temp(find(s_x1==1)),t_temp(find(s_x1==1)),'LineWidth',4); % bold line denotes the newly-built line
    highlight(pir(j),s_temp(i),t_temp(i),'EdgeColor','k','LineStyle','-','LineWidth',4);  % black line denotes the outage line
    highlight(pir(j),s_temp(find(s_y1(:,i)==1)),t_temp(find(s_y1(:,i)==1)),'EdgeColor','r','LineWidth',6,'LineStyle','-');
    labelnode(pir(j),xr(j),{[num2str(Load(xr(j))) ' Lost ' num2str(s_rt1(xr(j),yr(j)))]});
    for k=1:length(N_RLoads)
        labelnode(pir(j),N_RLoads(k),num2str(Load(N_RLoads(k))));
    end
    title(['Contigency ' num2str(i) ' happens, shedding load in node ', num2str(xr(j))]);
end