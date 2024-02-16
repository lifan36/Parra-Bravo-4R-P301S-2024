
%% find KCl mono respone

nThreshod=3;   % setting 2, 3, 05 5 std

hig700=figure(700);

clf
 
pos1=[0.1  0.57  0.8 0.40]; 
subplot('Position',pos1)
imagesc(Ca(:, :, 500), Ca_lim);
colormap('gray');
axis off; axis image;  
colorbar;
hold on
for k = 1 : size(ROI, 1)
	ROIxy =ROI{k};
	plot(ROIxy(:,2), ROIxy(:,1), 'color',CM(k,:), 'LineWidth', 2);
    text(mean(ROIxy(:,2)), mean(ROIxy(:,1)), num2str(k), 'Color','k','FontSize',8)

    clear ROIxy
end
hold off;

for i = 1:size(ROI_index, 2)
  
    mTrace=ROITrace_Unbleach(:, i);
  
    %%%%%%%
    ROIpeak.peakDF(:, i)=max(mTrace);
    ROIpeak.Time2Peak(:, i)=find(mTrace==ROIpeak.peakDF(:, i))/CamFreq;
    Index_halfPeak=find(mTrace>=(ROIpeak.peakDF(:, i)/2));
    ROIpeak.Time_half_rise(:, i)=min(Index_halfPeak)/CamFreq;
    ROIpeak.Time_half_decay(:, i)=max(Index_halfPeak)/CamFreq;


    AmpThreshold=mean(mTrace(90:100))+nThreshod*std(mTrace(90:100));  % using 90-100 as baseline
    indexActivity=find(mTrace>=AmpThreshold);
    %%% skipping the first few frames

    if isempty(indexActivity) %%% no activity than threshold
        ROIpeak.TimeDelay(:, i)=0;
        ROIpeak.duration(:, i)=0;

        %%%%%%%%%%
        ROIpeak.peakDF(:, i)=0;
        ROIpeak.Time2Peak(:, i)=0;
        ROIpeak.Index_halfPeak=0;
        ROIpeak.Time_half_rise(:, i)=0;
        ROIpeak.Time_half_decay(:, i)=0;

    else 
        indexActivity=indexActivity(find(indexActivity>200));
        if isempty(indexActivity)
            ROIpeak.TimeDelay(:, i)=0;
            ROIpeak.duration(:, i)=0;
            %%%%%%%%%%
            ROIpeak.peakDF(:, i)=0;
            ROIpeak.Time2Peak(:, i)=0;
            Index_halfPeak=0;
            ROIpeak.Time_half_rise(:, i)=0;
            ROIpeak.Time_half_decay(:, i)=0;

        else
            ROIpeak.TimeDelay(:, i)=min(indexActivity)/CamFreq;
            offset=max(indexActivity)/CamFreq;
            ROIpeak.duration(:, i)=offset-ROIpeak.TimeDelay(:, i);
        end
    end

    %%%%%%%%%%%%%%%%%%%% plot peak with 'x' marker
    %%%%% 'x' at 0 means no activity than th SD (default=3sd)

    pos1=[0.1  0.01+(i*0.48)/NumROI  0.8 0.48/NumROI];
    subplot('Position',pos1)
    plot(time, ROITrace_Unbleach(:, i), 'color',CM(i,:), 'LineWidth', 2);
    ylabel(['C' num2str(i)])
    if i==1
        box off
        xlabel('Time(s)')
        set(gca, 'box','off',  'YTickLabel',[],'ytick',[]);   

    else
        set(gca, 'box','off','XTickLabel',[],'xtick',[], 'XColor','none', 'YTickLabel',[],'ytick',[]);
    end
      xlim(Ca2_plotXlim) 
      ylim(Ca2_plotYlim)

   hold on
   plot(ROIpeak.Time2Peak(:, i),ROIpeak.peakDF(:, i),'xk')
 
end

 

hold off
%%%%%%%%%%%%%%%%%%%%%%
%% find KCl mono respone ROI_in
 
hig701=figure(701);

clf
 
pos1=[0.1  0.57  0.8 0.40]; 
subplot('Position',pos1)
imagesc(Ca(:, :, 500), Ca_lim);
colormap('gray');
axis off; axis image;  
colorbar;
hold on
for k = 1 : size(ROI_in_index, 2)
	ROIxy =ROI_in{k};
	plot(ROIxy(:,2), ROIxy(:,1), 'color',CM(k,:), 'LineWidth', 2);
    text(mean(ROIxy(:,2)), mean(ROIxy(:,1)), num2str(k), 'Color','k','FontSize',8)

    clear ROIxy
end
hold off;

for k = 1:size(ROI_in_index, 2)
  
    mTrace=ROI_inTrace_Unbleach(:, k);
  
    %%%%%%%
    ROIpeak_in.peakDF(:, k)=max(mTrace);
    ROIpeak_in.Time2Peak(:, k)=find(mTrace==ROIpeak_in.peakDF(:, k))/CamFreq;
    Index_halfPeak=find(mTrace>=(ROIpeak_in.peakDF(:, k)/2));
    ROIpeak_in.Time_half_rise(:, k)=min(Index_halfPeak)/CamFreq;
    ROIpeak_in.Time_half_decay(:, k)=max(Index_halfPeak)/CamFreq;


    AmpThreshold=mean(mTrace(90:100))+nThreshod*std(mTrace(90:100));  % using 90-100 as baseline
    indexActivity=find(mTrace>=AmpThreshold);
    %%% skipping the first few frames

    if isempty(indexActivity) %%% no activity than threshold
        ROIpeak_in.TimeDelay(:, k)=0;
        ROIpeak_in.duration(:, k)=0;

        %%%%%%%%%%
        ROIpeak_in.peakDF(:, k)=0;
        ROIpeak_in.Time2Peak(:, k)=0;
        Index_halfPeak=0;
        ROIpeak_in.Time_half_rise(:, k)=0;
        ROIpeak_in.Time_half_decay(:, k)=0;

    else 
        indexActivity=indexActivity(find(indexActivity>200));
        if isempty(indexActivity)
            ROIpeak_in.TimeDelay(:, k)=0;
            ROIpeak_in.duration(:, k)=0;
            %%%%%%%%%%
            ROIpeak_in.peakDF(:, k)=0;
            ROIpeak_in.Time2Peak(:, k)=0;
            Index_halfPeak=0;
            ROIpeak_in.Time_half_rise(:, k)=0;
            ROIpeak_in.Time_half_decay(:, k)=0;

        else
            ROIpeak_in.TimeDelay(:, k)=min(indexActivity)/CamFreq;
            offset=max(indexActivity)/CamFreq;
            ROIpeak_in.duration(:, k)=offset-ROIpeak_in.TimeDelay(:, k);
        end
    end

    %%%%%%%%%%%%%%%%%%%% plot peak with 'x' marker
    %%%%% 'x' at 0 means no activity than th SD (default=3sd)

    pos1=[0.1  0.01+(k*0.5)/NumROI  0.8 0.5/NumROI];
    subplot('Position',pos1)
    plot(time, ROI_inTrace_Unbleach(:, k), 'color',CM(k,:), 'LineWidth', 2);
    ylabel(['C' num2str(k)])
    if k==1
        box off
        xlabel('Time(s)')
        set(gca, 'box','off',  'YTickLabel',[],'ytick',[]);   

    else
        set(gca, 'box','off','XTickLabel',[],'xtick',[], 'XColor','none', 'YTickLabel',[],'ytick',[]);
    end
      xlim(Ca2_plotXlim) 
      ylim(Ca2_plotYlim)

   hold on
   plot(ROIpeak_in.Time2Peak(:, k),ROIpeak_in.peakDF(:, k),'xk')
 
end

hold off

%% find KCl mono respone ROI_ex
 
hig702=figure(702);

clf
 
pos1=[0.1  0.57  0.8 0.40]; 
subplot('Position',pos1)
imagesc(Ca(:, :, 500), Ca_lim);
colormap('gray');
axis off; axis image;  
colorbar;
hold on
for j = 1 : size(ROI_ex, 1)
	ROIxy =ROI_ex{j};
	plot(ROIxy(:,2), ROIxy(:,1), 'color',CM(k+j,:), 'LineWidth', 2);
    text(mean(ROIxy(:,2)), mean(ROIxy(:,1)), num2str(k+j), 'Color','k','FontSize',8)

    clear ROIxy
end
hold off;

for j = 1:size(ROI_ex, 1)
  
    mTrace=ROI_exTrace_Unbleach(:, j);
  
    %%%%%%%
    ROIpeak_ex.peakDF(:, j)=max(mTrace);
    ROIpeak_ex.Time2Peak(:, j)=find(mTrace==ROIpeak_ex.peakDF(:, j))/CamFreq;
    Index_halfPeak=find(mTrace>=(ROIpeak_ex.peakDF(:, j)/2));
    ROIpeak_ex.Time_half_rise(:, j)=min(Index_halfPeak)/CamFreq;
    ROIpeak_ex.Time_half_decay(:, j)=max(Index_halfPeak)/CamFreq;


    AmpThreshold=mean(mTrace(90:100))+nThreshod*std(mTrace(90:100));  % using 90-100 as baseline
    indexActivity=find(mTrace>=AmpThreshold);
    %%% skipping the first few frames

    if isempty(indexActivity) %%% no activity than threshold
        ROIpeak_ex.TimeDelay(:, j)=0;
        ROIpeak_ex.duration(:, j)=0;

        %%%%%%%%%%
        ROIpeak_ex.peakDF(:, j)=0;
        ROIpeak_ex.Time2Peak(:, j)=0;
        Index_halfPeak=0;
        ROIpeak_ex.Time_half_rise(:, j)=0;
        ROIpeak_ex.Time_half_decay(:, j)=0;

    else 
        indexActivity=indexActivity(find(indexActivity>200));
        if isempty(indexActivity)
            ROIpeak_ex.TimeDelay(:, j)=0;
            ROIpeak_ex.duration(:, j)=0;
            %%%%%%%%%%
            ROIpeak_ex.peakDF(:, j)=0;
            ROIpeak_ex.Time2Peak(:, j)=0;
            Index_halfPeak=0;
            ROIpeak_ex.Time_half_rise(:, j)=0;
            ROIpeak_ex.Time_half_decay(:, j)=0;

        else
            ROIpeak_ex.TimeDelay(:, j)=min(indexActivity)/CamFreq;
            offset=max(indexActivity)/CamFreq;
            ROIpeak_ex.duration(:, j)=offset-ROIpeak_ex.TimeDelay(:, j);
        end
    end

    %%%%%%%%%%%%%%%%%%%% plot peak with 'x' marker
    %%%%% 'x' at 0 means no activity than th SD (default=3sd)

    pos1=[0.1  0.01+(j*0.5)/NumROI  0.8 0.5/NumROI];
    subplot('Position',pos1)
    plot(time, ROI_exTrace_Unbleach(:, j), 'color',CM(k+j,:), 'LineWidth', 2);
    ylabel(['C' num2str(k+j)])
    if j==1
        box off
        xlabel('Time(s)')
        set(gca, 'box','off',  'YTickLabel',[],'ytick',[]);   

    else
        set(gca, 'box','off','XTickLabel',[],'xtick',[], 'XColor','none', 'YTickLabel',[],'ytick',[]);
    end
      xlim(Ca2_plotXlim) 
      ylim(Ca2_plotYlim)

   hold on
     plot(ROIpeak_ex.Time2Peak(:, j), ROIpeak_ex.peakDF(:, j),'xk')
 
end

hold off
%%
hfig703=figure(703)

clf
    set(gcf,'FileName','Trace_in_ex','PaperPosition',[0.25 0.25 11.19 7.768],...
    'PaperSize',[11.69 8.268 ],'PaperType','a4letter');
 
 
pos1=[0.1  0.57  0.4 0.40]; 
subplot('Position',pos1)
imagesc(Ca(:, :, 500), Ca_lim);
colormap('gray');
axis off; axis image;  
colorbar;
hold on
for k = 1 : size(ROI_in_index, 2)
	ROIxy =ROI_in{k};
	plot(ROIxy(:,2), ROIxy(:,1), 'color',CM(k,:), 'LineWidth', 2);
    text(mean(ROIxy(:,2)), mean(ROIxy(:,1)), num2str(k), 'Color','k','FontSize',8)

    clear ROIxy
end
hold off;

pos2=[0.5  0.57  0.4 0.40]; 
subplot('Position',pos2)
imagesc(Ca(:, :, 500), Ca_lim);
colormap('gray');
axis off; axis image;  
colorbar;
hold on
for j = 1 : size(ROI_ex, 1)
	ROIxy =ROI_ex{j};
	plot(ROIxy(:,2), ROIxy(:,1), 'color',CM(k+j,:), 'LineWidth', 2);
    text(mean(ROIxy(:,2)), mean(ROIxy(:,1)), num2str(k+j), 'Color','k','FontSize',8)

    clear ROIxy
end
hold off;



pos3=[0.1  0.1  0.8 0.4];


for k = 1:size(ROI_in_index, 2)
    %%%%%%%%%%%%%%%%%%%% plot peak with 'x' marker
    %%%%% 'x' at 0 means no activity than th SD (default=3sd)

  %  pos1=[0.1  0.01+(k*0.5)/NumROI  0.8 0.5/NumROI];
    subplot('Position',pos3)
    plot(time, ROI_inTrace_Unbleach(:, k), 'color',CM(k,:), 'LineWidth', 2);
   % ylabel(['C' num2str(k)])
   ylabel('\Delta F/F',  'FontSize',12);

%     if k==1
%         box off
        xlabel('Time(s)')
       % set(gca, 'box','off',  'YTickLabel',[],'ytick',[]);   
        set(gca, 'box','off');   

%     else
%         set(gca, 'box','off','XTickLabel',[],'xtick',[], 'XColor','none', 'YTickLabel',[],'ytick',[]);
%     end
      xlim(Ca2_plotXlim) 
      ylim(Ca2_plotYlim)

   hold on
   plot(ROIpeak_in.Time2Peak(:, k),ROIpeak_in.peakDF(:, k),'xk')
end 
hold on




for j = 1:size(ROI_ex, 1)

    plot(time, ROI_exTrace_Unbleach(:, j), 'color',CM(k+j,:), 'LineWidth', 2);
%     ylabel(['C' num2str(k+j)])
%     if j==1
%         box off
%         xlabel('Time(s)')
%         set(gca, 'box','off',  'YTickLabel',[],'ytick',[]);   
% 
%     else
%         set(gca, 'box','off','XTickLabel',[],'xtick',[], 'XColor','none', 'YTickLabel',[],'ytick',[]);
%     end
      xlim(Ca2_plotXlim) 
      ylim(Ca2_plotYlim)

   hold on
     plot(ROIpeak_ex.Time2Peak(:, j), ROIpeak_ex.peakDF(:, j),'xk')
 
end

hold on













%%




File_name7 = ['MonoPeak_' openmkpath((length(openmkpath)-22):(length(openmkpath)-1))...
    '-' openfile1(1:(lenght_f1-4)) openROI((length(openROI)-7):(length(openROI)-4) ) '.mat'] 

File_name_fig700 = ['Fig_' openfile1(1:(lenght_f1-4))  '-MonoPeak' '.fig']; 
File_name_fig701 = ['Fig_' openfile1(1:(lenght_f1-4))  '-MonoPeak_in' '.fig'] ;
File_name_fig702 = ['Fig_' openfile1(1:(lenght_f1-4))  '-MonoPeak_ex' '.fig'] ;


save(File_name7, 'ROIpeak', 'ROIpeak_in', 'ROIpeak_ex', ...
    'RFPTotalFluo_in', 'RFPMaxFluo_in', 'RFPTotalFluo_ex', 'RFPMaxFluo_ex');
savefig(hig700,File_name_fig700);
savefig(hig701,File_name_fig701);
savefig(hig702,File_name_fig702);



File_name_fig703 = ['Fig_' openfile1(1:(lenght_f1-4))  '-MonoPeak_Trace' '.eps']; 
exportgraphics(hfig703,File_name_fig703,'ContentType','vector')