 hfig102=figure (102);
    clf
        set(gcf,'FileName','Trace','PaperPosition',[0.25 0.25 11.19 7.768],...
    'PaperSize',[11.69 8.268 ],'PaperType','a4letter');
 
    Time=time;
ROItraceRaw=ROI_exTrace_Unbleach; 
mTrace=mean(ROItraceRaw, 2);

    set(gcf,'FileName','Sat Fig','PaperPosition',[0.25 0.25 11.19 7.768],...
    'PaperSize',[11.69 8.268 ],'PaperType','a4letter');
    
    
    %%%% standard error bar
    eb1=errorbar(Time,mTrace,std(ROItraceRaw,[],2)/sqrt(size(ROItraceRaw, 2)-1),'Color', [0.5, 0.5, 0.5], ... % [0.5 0.5 0.5], ...
      'LineStyle','none','MarkerSize',1);
    
    %%%% std bar
%     eb1=errorbar(Time,mTrace,std(ROItraceRaw, [], 2),'Color', [0.5, 0.5, 0.5], ... % [0.5 0.5 0.5], ...
%     'LineStyle','none','MarkerSize',1);
    
    
    eb1.LineWidth=0.3;
    eb1.CapSize=2; hold on
    
    p1=plot(Time, mTrace,'k', 'LineWidth',2);
%     axis([-inf 300 -0.5 0.9]);
    
    %title('Ca2+');
%     xlabel('Time(Sec)');
%     ylabel('\DeltaF/F');
%     
 hold on


    ROItraceRaw=ROI_inTrace_Unbleach; 
   mTrace=mean(ROItraceRaw, 2);
    %%% standard error bar
    eb1=errorbar(Time,mTrace,std(ROItraceRaw,[],2)/sqrt(size(ROItraceRaw, 2)-1),'Color', [0.9, 0.5, 0.5], ... % [0.5 0.5 0.5], ...
    'LineStyle','none','MarkerSize',1);

    % %%%% std bar
    % eb1=errorbar(Time,mTrace,std(ROItraceRaw, [], 2),'Color', [0.5, 0.5, 0.5], ... % [0.5 0.5 0.5], ...
    %   'LineStyle','none','MarkerSize',1);
    
    
    eb1.LineWidth=0.3;
    eb1.CapSize=2; hold on
     
    p1=plot(Time, mTrace,'r', 'LineWidth',2);
  
    hold on

 %   p2=plot([min(time) max(time)], [mean(mTrace(90:100))+3*std(mTrace(90:100)) mean(mTrace(90:200))+3*std(mTrace(90:200))], '--r');
%     axis([-inf 400 -0.2 0.6]);
    %title('Ca2+');
%     xlabel('Time(Sec)');
%     ylabel('\Delta F/F');
   %% 
   xlim([0 400])
  % ylim([-0.05 0.15])
   legend({'',' (-) tau','', ' (+) tau'})
 legend boxoff 
 set(gca, 'color', 'none', 'box', 'off');
 

 %%
%File_name_fig102 = ['Fig_' openfile1(1:(lenght_f1-4))  '-MonoPeak_Trace_In_vs_ex' '.eps']; 

File_name_fig102 = ['Fig_' 'CaP301SHalo006_-MonoPeak_Trace_In_vs_ex' '.eps']; 
exportgraphics(hfig102,File_name_fig102,'ContentType','vector')