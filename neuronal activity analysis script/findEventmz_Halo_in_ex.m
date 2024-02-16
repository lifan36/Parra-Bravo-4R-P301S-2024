 %function findEventmz

%%
hfig20=figure(600);
clf
% sz=size(I);
% CC=bwconncomp(I);
% CC_index=CC.PixelIdxList;
% for i=1:size(CC_index, 2)
%    [xx yy] =ind2sub(sz, CC_index{i});
%    ROI_index{i}=[xx yy];
% end

for i = 1:size(ROI_index, 2)
    for n=1: size(Ca2, 3)
        ROIxy=ROI_index{i};
        Ca2ROI=Ca(ROIxy(:, 1), ROIxy(:, 2), n);
       T(n)= mean2(Ca2ROI);
       
    end
    
    plot(time, T, 'color',CM(i,:)); 
    hold on     
 
    ROITraceRaw(:, i) = T;
    %%%%%%%%% exponentially weighted moving average for noise filtering
    %%%%%%%%% (photobleech correction)
    
     s=exp2fit(time, T, 1);
     fun = @(s,time) s(1)+s(2)*exp(-time/s(3));
     Trace_ExpFit=fun(s,time);
    
     ROITraceRaw_Unbleach(:, i)=(T'./Trace_ExpFit)*Trace_ExpFit(1);
     ss(:, i)=s;

     %%%%%%%%%% covnert to deltaF trace
     ROITrace_Unbleach(:, i)=ROITraceRaw_Unbleach(:, i)./...
         mean(ROITraceRaw_Unbleach(baseLineIndex, i))-1;

     
end
subplot(2,1,1)
for i = 1:size(ROI_index, 2)
  plot(time, ROITraceRaw(:, i), 'color',CM(i,:)); 
  hold on
end
% plot(time, mean(ROITraceRaw, 2), '-k', 'LineWidth',4)
% hold on
plot([min(time) max(time)], [0 0], '--k');
xlim([min(time) max(time)])
xlabel('Time(sec)', 'FontSize', 12);
ylabel('Fluo Intensity')
title('Raw Trace')
hold off

subplot(2,1,2)
for i = 1:size(ROI_index, 2)
  plot(time, ROITraceRaw_Unbleach(:, i), 'color',CM(i,:)); 
  hold on
end

% plot(time, mean(ROITraceRaw_Unbleach, 2), '-k', 'LineWidth',4)
% hold on
plot([min(time) max(time)], [0 0], '--k');
xlim([min(time) max(time)])
xlabel('Time(sec)', 'FontSize', 12);
ylabel('Fluo Intensity')
title('Raw Trace-Unbleach')
hold off

%%
hfig2=figure(300);
clf
for i = 1:size(ROI_index, 2)

    
    plot(time, ROITrace_Unbleach(:, i), 'color',CM(i,:)); 
    hold on

end
plot(time, mean(ROITrace_Unbleach, 2), '-k', 'LineWidth',4)
hold on
plot([min(time) max(time)], [0 0], '--k');
xlim(Ca2_plotXlim)
ylim(Ca2_plotYlim)
xlabel('Time(sec)', 'FontSize', 12);
ylabel('\Delta F/F',  'FontSize',12);
hold off

%% plot traces1 
hfig=figure(601);
clf
 
 
pos1=[0.1  0.60  0.8 0.38]; 
subplot('Position',pos1)
imagesc(Ca(:, :, 10), Ca_lim);
colormap('gray');
axis off; axis image;  
colorbar;
hold on
for i = 1 : size(ROI, 1)
	ROIxy =ROI{i};
	plot(ROIxy(:,2), ROIxy(:,1), 'color',CM(i,:), 'LineWidth', 2);
    text(mean(ROIxy(:,2)), mean(ROIxy(:,1)), num2str(i), 'Color','k','FontSize',8)

    clear ROIxy
end
 
% pos1=[0.1  (m*0.5)/m   0.8 0.5/m];
% subplot('Position',pos1)
%     plot(time, ROITrace(1, :), 'color',CM(1,:), 'LineWidth', 2);
%     ylabel(['C' num2str(1)])
%     box off
for i = 1:size(ROI_index, 2)
    pos1=[0.1  0.04+(i*0.48)/NumROI  0.8 0.48/NumROI];
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
 
end
%%%%%%%%%%%%%
%% get ROI in trace

% CC1=bwconncomp(I_in);
% CC1_index=CC1.PixelIdxList;
% for k=1:size(CC1_index, 2)
%    [xx yy] =ind2sub(sz, CC1_index{k});
%    ROI_in_index{k}=[xx yy];
% end

hfig301=figure(301);
clf

for k = 1:size(ROI_in_index, 2)
    for n=1: size(Ca2, 3)
         
       ROIxy=ROI_in_index{k};
       Ca2ROI_in=Ca(ROIxy(:, 1), ROIxy(:, 2), n);
       T1(n)= mean2(Ca2ROI_in);
       
       
    end
    
   
    ROI_inTraceRaw(:, k) = T1;
    %%%%%%%%% exponentially weighted moving average for noise filtering
    %%%%%%%%% (photobleech correction)
    
     s=exp2fit(time, T1, 1);
     fun = @(s,time) s(1)+s(2)*exp(-time/s(3));
     Trace_ExpFit=fun(s,time);
    
     ROI_inTraceRaw_Unbleach(:, k)=(T1'./Trace_ExpFit)*Trace_ExpFit(1);
     ss_in(:, k)=s;

     %%%%%%%%%% covnert to deltaF trace
     ROI_inTrace_Unbleach(:, k)=ROI_inTraceRaw_Unbleach(:, k)./...
         mean(ROI_inTraceRaw_Unbleach(baseLineIndex, k))-1;
    plot(time, ROI_inTrace_Unbleach(:, k), 'color',CM(k,:)); 
    hold on   
     
end
plot(time, mean(ROI_inTrace_Unbleach, 2), '-k', 'LineWidth',4)
hold on
plot([min(time) max(time)], [0 0], '--k');
xlim(Ca2_plotXlim)
ylim(Ca2_plotYlim)
xlabel('Time(sec)', 'FontSize', 12);
ylabel('\Delta F/F',  'FontSize',12);
hold off

%% plot traces1 inclusion
hfig401=figure(401);
clf
 
pos1=[0.1  0.50  0.8 0.45]; 
subplot('Position',pos1)
imagesc(Ca(:, :, 10), Ca_lim);
colormap('gray');
axis off; axis image;  
colorbar;
hold on
for k = 1 : size(ROI_in, 1)
	ROIxy_in =ROI_in{k};
	plot(ROIxy_in(:,2), ROIxy_in(:,1), 'color',CM(k,:), 'LineWidth', 2);
    text(mean(ROIxy_in(:,2)), mean(ROIxy_in(:,1)), num2str(k), 'Color','k','FontSize',8)

    clear ROIxy
end
hold off;
% %%
% hfig4=figure(4001);
% clf
%   pos1=[0.1  (NumROI*0.5)/NumROI   0.8 0.5/NumROI];
%   subplot('Position',pos1)
%     plot(time, ROITrace(1, :), 'color',CM(1,:), 'LineWidth', 2);
%     ylabel(['C' num2str(1)])
%     box off
for k = 1:size(ROI_in_index, 2)
   pos1=[0.1  0.04+(k*0.5)/NumROI  0.8 0.5/NumROI];
 %   pos1=[0.1  0.04+(i*0.85)/NumROI  0.8 0.85/NumROI];

    subplot('Position',pos1)
    plot(time, ROI_inTrace_Unbleach(:, k), 'color',CM(k,:), 'LineWidth', 1);
    
    
    ylabel(['C' num2str(k)])
   if k==1
        box off
        xlabel('Time(s)')
    else
        set(gca, 'box','off','XTickLabel',[],'xtick',[])
    end
      % xlim([0 480]) 
      ylim(Ca2_plotYlim)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%
%% get ROI ex trace

% CC2=bwconncomp(I_ex);
% CC2_index=CC2.PixelIdxList;
% for j=1:size(CC2_index, 2)
%    [xx yy] =ind2sub(sz, CC2_index{j});
%    ROI_ex_index{j}=[xx yy];
% end

hfig302=figure(302);
clf

for j = 1:size(ROI_ex_index, 2)
    for n=1: size(Ca2, 3)
         
       ROIxy=ROI_ex_index{j};
       Ca2ROI_ex=Ca(ROIxy(:, 1), ROIxy(:, 2), n);
       T2(n)= mean2(Ca2ROI_ex);
       
    end
    
   
    ROI_exTraceRaw(:, j) = T2;
    %%%%%%%%% exponentially weighted moving average for noise filtering
    %%%%%%%%% (photobleech correction)
    
     s=exp2fit(time, T2, 1);
     fun = @(s,time) s(1)+s(2)*exp(-time/s(3));
     Trace_ExpFit=fun(s,time);
    
     ROI_exTraceRaw_Unbleach(:, j)=(T2'./Trace_ExpFit)*Trace_ExpFit(1);
     ss_ex(:, j)=s;

     %%%%%%%%%% covnert to deltaF trace
     ROI_exTrace_Unbleach(:, j)=ROI_exTraceRaw_Unbleach(:, j)./...
         mean(ROI_exTraceRaw_Unbleach(baseLineIndex, j))-1;
    plot(time, ROI_exTrace_Unbleach(:, j), 'color',CM(k+j,:)); 
    hold on   
     
end
plot(time, mean(ROI_exTrace_Unbleach, 2), '-k', 'LineWidth',4)
hold on
plot([min(time) max(time)], [0 0], '--k');
xlim(Ca2_plotXlim)
ylim(Ca2_plotYlim)
xlabel('Time(sec)', 'FontSize', 12);
ylabel('\Delta F/F',  'FontSize',12);
hold off

%% plot traces1 inclusion
hfig402=figure(402);
clf
 
pos1=[0.1  0.50  0.8 0.45]; 
subplot('Position',pos1)
imagesc(Ca(:, :, 10), Ca_lim);
colormap('gray');
axis off; axis image;  
colorbar;
hold on
for j = 1 : size(ROI_ex, 1)
	ROIxy_ex =ROI_ex{j};
	plot(ROIxy_ex(:,2), ROIxy_ex(:,1), 'color',CM(k+j,:), 'LineWidth', 2);
    text(mean(ROIxy_ex(:,2)), mean(ROIxy_ex(:,1)), num2str(k+j), 'Color','k','FontSize',8)

    clear ROIxy
end
hold off;
% %%
% hfig4=figure(4001);
% clf
%   pos1=[0.1  (NumROI*0.5)/NumROI   0.8 0.5/NumROI];
%   subplot('Position',pos1)
%     plot(time, ROITrace(1, :), 'color',CM(1,:), 'LineWidth', 2);
%     ylabel(['C' num2str(1)])
%     box off
for j = 1:size(ROI_ex, 1)
   pos1=[0.1  0.04+(j*0.5)/NumROI  0.8 0.5/NumROI];
 %   pos1=[0.1  0.04+(i*0.85)/NumROI  0.8 0.85/NumROI];

    subplot('Position',pos1)
    plot(time, ROI_exTrace_Unbleach(:, j), 'color',CM(k+j,:), 'LineWidth', 1);
    
    
    ylabel(['C' num2str(k+j)])
   if j==1
        box off
        xlabel('Time(s)')
    else
        set(gca, 'box','off','XTickLabel',[],'xtick',[])
    end
      % xlim([0 480]) 
      ylim(Ca2_plotYlim)
end
















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find events total

 
    [T,m] = size(ROITraceRaw) 


    metric='mean'
    tau=[1 2]


    Fbar = ROITraceRaw';
    for t = 1:T
       if strcmp(metric, 'mean')
          Fbar(:,t) = mean(Fbar(:,max(1,t-tau(1)/2):min(T,t+tau(1)/2)),2);
       elseif strcmp(metric, 'median')
          Fbar(:,t) = median(Fbar(:,max(1,t-tau(1)/2):min(T,t+tau(1)/2)),2);
       end
    end
    F0 = Fbar;
    for t =  1:T
       F0(:,t) = min(Fbar(:,max(1,t-tau(2)):t),[],2);
    end
    smoothTraces = (Fbar-F0)./F0;
    clear Fbar F0;
    
    dff = smoothTraces;

 

%% plot traces event

hfig602=figure(602);
clf
 
% pos1=[0.1  0.57  0.8 0.45]; 
% subplot('Position',pos1)
% imagesc(Ca(:, :, 500), Ca_lim);
% colormap('gray');
% axis off; axis image;  
% colorbar;
% hold on
% for k = 1 : size(ROI, 1)
% 	ROIxy =ROI{k};
% 	plot(ROIxy(:,2), ROIxy(:,1), 'color',CM(k,:), 'LineWidth', 2);
%     clear ROIxy
% end
% hold off;
% pos1=[0.1  (m*0.5)/m   0.8 0.5/m];
% subplot('Position',pos1)
%     plot(time, ROITrace(1, :), 'color',CM(1,:), 'LineWidth', 2);
%     ylabel(['C' num2str(1)])
%     box off
 
for i = 1:size(ROI_index,2)
    pos1=[0.1  0.04+(i*0.8)/NumROI  0.85 0.8/NumROI];
    subplot('Position',pos1)
    plot(time, dff(i, :), 'color',CM(i,:), 'LineWidth', 0.2);
    ylabel(['C' num2str(i)])
    if i==1
        box off
        xlabel('Time(s)')
    else
        set(gca, 'box','off','XTickLabel',[],'xtick',[])
    end
      %  xlim([0 480]) 
end
%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find events ROI in

 
    [T,m] = size(ROI_inTraceRaw) 


    metric='mean'
    tau=[1 2]


    Fbar = ROI_inTraceRaw';
    for t = 1:T
       if strcmp(metric, 'mean')
          Fbar(:,t) = mean(Fbar(:,max(1,t-tau(1)/2):min(T,t+tau(1)/2)),2);
       elseif strcmp(metric, 'median')
          Fbar(:,t) = median(Fbar(:,max(1,t-tau(1)/2):min(T,t+tau(1)/2)),2);
       end
    end
    F0 = Fbar;
    for t =  1:T
       F0(:,t) = min(Fbar(:,max(1,t-tau(2)):t),[],2);
    end
    smoothTraces = (Fbar-F0)./F0;
    clear Fbar F0;
    
    dff_in = smoothTraces;
%%%%%%%%%%%%%%%%%%%%%
 

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find events ROI ex

 
    [T,m] = size(ROI_exTraceRaw) 


    metric='mean'
    tau=[1 2]


    Fbar = ROI_exTraceRaw';
    for t = 1:T
       if strcmp(metric, 'mean')
          Fbar(:,t) = mean(Fbar(:,max(1,t-tau(1)/2):min(T,t+tau(1)/2)),2);
       elseif strcmp(metric, 'median')
          Fbar(:,t) = median(Fbar(:,max(1,t-tau(1)/2):min(T,t+tau(1)/2)),2);
       end
    end
    F0 = Fbar;
    for t =  1:T
       F0(:,t) = min(Fbar(:,max(1,t-tau(2)):t),[],2);
    end
    smoothTraces = (Fbar-F0)./F0;
    clear Fbar F0;
    
    dff_ex= smoothTraces;
%%%%%%%%%%%%%%%%%%%%%













%%% conver ddf to spikes
 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %% setp 5 run correlation maps
%         tic
%         load('spikes.mat')
% %  upsampling from 10Hz to 30Hz
%         for i=1:23
%             spikes{:,i}=interp(spikes{:, i}, 3);
%         end
%         % 
%  %         addpath(genpath(toolbox)) %Add the toolbox to the matlab working directory when ever you begin a new session
% 
%          fps=30
%         height=1 % background intensity
% %         dff=dff';
% 
%          [m,T]=size(dff) 
%        
%         for ii = 1:m
%             x = dff(ii,:);
%          %   x = interp1(1:T,x,1:fps/30:T);
%             dff1(ii,:) = x;
%             x = dff1(ii,:);
%             parfor i=1:length(spikes)
%                 snippet = spikes{i};
%                 L = length(snippet);
%                 C = zeros(size(x));       
%                                
%                         
%                 for j=1:length(x)-(L-1)
%                     x_snippet = x(j:j+L-1);
%                     if(range(x_snippet)>height)
%                         R = corrcoef(x_snippet,snippet);
%                         C(j) = R(1,2);
%                         
%                         if j == length(x)-(L-1)
%                             for j1 = length(x)-(L-2):length(x)-round(L/2)
%                                 x_snippet = x(j1:end);
%                                 R = corrcoef(x_snippet,snippet(1:length(x_snippet)));
%                                 C(j1) = R(1,2);
%                             end
%                         end
%                     end
%                  end
%              Call(i,:) = C;
%             end
%  
% 
% 
% 
% 
% 
% 
%             Ca3{:,ii} = Call;
%         end
%         %%
%         [m,T]=size(dff1);
%         
%         h = figure (603); %;set(h, 'Visible', 'off');
%         set(h, 'PaperPosition',[0.25 0.25 11 7],'PaperSize',[11 7],...
%            'PaperType','a4letter');
% 
%          clf
%         subplot(1,3,1)
%         plot_counter = 0;
%         for ii = 1:m
%           plot((1:T)./fps,dff1(ii,:)+plot_counter)
%           hold on
%           plot_counter = plot_counter+max(max(dff1))-min(min(dff1));
%         end
%         
%         ylim([min(min(dff1)) plot_counter])
%         ticc = max(max(dff1))-min(min(dff1));
%         yticks(max(max(dff1)):ticc:ticc*m)
%         yticklabels(1:1:m)
%         xlim([0 (size(dff1,2))./fps])
%         set(gca,'Fontsize',10) 
%         title(['height: ', num2str(height)], 'Fontsize', 20)
%         
%         subplot(1,3,2)
%         imagesc(flipud(cat(1,Ca3{:})))
%         tmp = flipud(cat(1,Ca3{:}));
%         yticks(1:size(Ca3{1},1):size(tmp,1))
%         yticklabels(sort((1:1:m),'descend'))
%         set(gca,'Fontsize',10) 
%         title([num2str(size(Ca3{1},1)),' motifs'], 'Fontsize', 20)
%         h =colorbar;
%         t=get(h,'Limits');
%         colorbar   %off
%         
%         subplot(1,3,3)
%         for i = 1:m, Caa(i,:) = max(Ca3{i}); end
%         Caa(Caa<0.6) = 0;
%         imagesc(flipud(Caa))
%         yticks(1:1:m)
%         yticklabels(sort((1:1:m),'descend'))
%         set(gca,'Fontsize',10) 
%         title(['corr of ', num2str(numel(spikes)),' motifs >0.6'], 'Fontsize', 20)
%         %h = colorbar;
%         caxis([t(1) t(2)])
%          colorbar
%         
%          set(gcf,'InvertHardCopy', 'off')
%          set(gcf,'Position',1.0e+03 *[0.0010    0.0410    2.5600    1.3273])
%          set(gcf,'PaperPosition',[0 0 20 m])
%         %%
%   %       save([fileName(1:end-7),num2str(height),'corr.mat'],'Ca3','dff1');
%   %    saveas(gcf,[fileName(1:end-7),num2str(height),'corr.jpg'])	
% 		        
%     %    close all
%         toc
        %end
        
        %%

% %% step 6A synchroniccity
%





%% Results 
% based on percentile based threshold (mean+1.5-2 std): not robust enouth
% to correctly idnetify both small and large amplitude transients)
% 
NumStd=3

hfig604=figure(604);
clf
for i=1:size(dff, 1)

    threshold_Peak=mean(dff(i,:))+NumStd*std(dff(i,:), 1);

    [spks, index] = findpeaks(dff(i,:),'MinPeakHeight',threshold_Peak);
    pos1=[0.1  0.04+(i*0.9)/NumROI  0.85 0.9/NumROI];
    subplot('Position',pos1)
    locs=time(index);
    amp=ROITrace_Unbleach(index, i); %%%%% amp from Ca_Unbleach trace
    
    plot(time, dff(i, :), 'color',CM(i,:), 'LineWidth', 0.2);
    hold on

    plot(locs,spks,'xk')
%       ylim([-0.001 0.08]);



    
    ylabel(['C' num2str(i)])
    if i==1
        box off
        xlabel('Time(s)')
    else
        set(gca, 'box','off','XTickLabel',[],'xtick',[])
    end
    hold on
    spks=spks';

    Cell_spks(i).spk=spks;
    Cell_spks(i).locs=locs;
    Cell_spks(i).index=index';
    Cell_spks(i).amp=amp;
    Cell_spks(i).Ave_amp=mean(amp,1);
    Cell_spks(i).Ave_spk=mean(spks, 1);
    Cell_spks(i).Max_amp=max(amp);
    Cell_spks(i).Max_spk=max(spks);
    Cell_spks(i).Num_spk=size(spks, 1);
    Cell_spks(i).Freq_spk=size(spks, 1)./max(time);

end
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Journal of Neuroscience Methods
% Volume 349, 1 February 2021, 109041
% synchronous spikes are identified using a threshold based on predefined percentile (mean + 1.5-2 X standard deviation, s.d.). 
% Synchronicity rate is quantified as the number of the detected synchronous spikes in one minute. 
% Quantification of single-neuron amplitude and frequency. The amplitude of peaks in each trace is defined as the mean value of △F/F0 of individual peak. 
% The frequency of peaks in each trace is defined as the number of detected peaks in one minute.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Synchronous fire rate (min to /s)
Total_trace=mean(dff, 1);
threshold_Peak=mean(Total_trace)+NumStd*std(Total_trace, 1)
[Nw_spks, Nw_index] = findpeaks(Total_trace,'MinPeakHeight',threshold_Peak);
%%
hfig605=figure(605)
clf
plot(time, dff', 'color',[0.5 0.5 0.5], 'LineWidth', 0.2);
hold on
plot(time, Total_trace, '-b', 'LineWidth', 0.5);
Nw_locs=time(Nw_index);
hold on

plot([min(time) max(time)], [threshold_Peak threshold_Peak], 'color',[0.5 0.5 0.5], 'LineWidth', 0.2);
hold on
plot(Nw_locs,Nw_spks,'xk');
hold off
xlabel('Time(sec)', 'FontSize', 12);
ylabel('\Delta F/F',  'FontSize',12);
title('Synchronous Peaks')

NetworkAna.spks=Nw_spks;
NetworkAna.index=Nw_index;
NetworkAna.trace=Total_trace;
NetworkAna.Num_spks=size(Nw_spks, 2);
NetworkAna.Freq_spks=size(Nw_spks, 2)./max(time);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   result ROI in
hfig606=figure(606);
clf
NumROI1=size(dff_in, 1)
for k=1:size(dff_in, 1)

    threshold_Peak=mean(dff_in(k,:))+NumStd*std(dff_in(k,:), 1);

    [spks, index] = findpeaks(dff_in(k,:),'MinPeakHeight',threshold_Peak);
    pos1=[0.1  0.005+(k*0.9)/NumROI1  0.85 0.9/NumROI1];
    subplot('Position',pos1)
    locs=time(index);
    amp=ROI_inTrace_Unbleach(index, k); %%%%% amp from Ca_Unbleach trace
    
    plot(time, dff_in(k, :), 'color',CM(k,:), 'LineWidth', 0.2);
    hold on

    plot(locs,spks,'xk')
%       ylim([-0.001 0.08]);



    
    ylabel(['C' num2str(k)])
    if k==1
        box off
        xlabel('Time(s)')
    else
        set(gca, 'box','off','XTickLabel',[],'xtick',[])
    end
    hold on
    spks=spks';

    Cell_spks_in(k).spk=spks;
    Cell_spks_in(k).locs=locs;
    Cell_spks_in(k).index=index';
    Cell_spks_in(k).amp=amp;
    Cell_spks_in(k).Ave_amp=mean(amp,1);
    Cell_spks_in(k).Ave_spk=mean(spks, 1);
    Cell_spks_in(k).Max_amp=max(amp);
    Cell_spks_in(k).Max_spk=max(spks);
    Cell_spks_in(k).Num_spk=size(spks, 1);
    Cell_spks_in(k).Freq_spk=size(spks, 1)./max(time);

end
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Journal of Neuroscience Methods
% Volume 349, 1 February 2021, 109041
% synchronous spikes are identified using a threshold based on predefined percentile (mean + 1.5-2 X standard deviation, s.d.). 
% Synchronicity rate is quantified as the number of the detected synchronous spikes in one minute. 
% Quantification of single-neuron amplitude and frequency. The amplitude of peaks in each trace is defined as the mean value of △F/F0 of individual peak. 
% The frequency of peaks in each trace is defined as the number of detected peaks in one minute.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Synchronous fire rate (min to /s)
Total_trace_in=mean(dff_in, 1);
threshold_Peak=mean(Total_trace_in)+NumStd*std(Total_trace_in, 1)
[Nw_spks, Nw_index] = findpeaks(Total_trace_in,'MinPeakHeight',threshold_Peak);

hfig607=figure(607)
clf
plot(time, dff_in', 'color',[0.5 0.5 0.5], 'LineWidth', 0.2);
hold on
plot(time, Total_trace_in, '-b', 'LineWidth', 0.5);
Nw_locs=time(Nw_index);
hold on

plot([min(time) max(time)], [threshold_Peak threshold_Peak], 'color',[0.5 0.5 0.5], 'LineWidth', 0.2);
hold on
plot(Nw_locs,Nw_spks,'xk');
hold off
xlabel('Time(sec)', 'FontSize', 12);
ylabel('\Delta F/F',  'FontSize',12);
title('Synchronous Peaks')

NetworkAna_in.spks=Nw_spks;
NetworkAna_in.index=Nw_index;
NetworkAna_in.trace=Total_trace;
NetworkAna_in.Num_spks=size(Nw_spks, 2);
NetworkAna_in.Freq_spks=size(Nw_spks, 2)./max(time);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   result ROI ex
hfig608=figure(608);
clf
NumROI2=size(dff_ex, 1)
for j=1:size(dff_ex, 1)

    threshold_Peak=mean(dff_ex(j,:))+NumStd*std(dff_ex(j,:), 1);

    [spks, index] = findpeaks(dff_ex(j,:),'MinPeakHeight',threshold_Peak);
    pos1=[0.1  0.005+(j*0.8)/NumROI2  0.85 0.8/NumROI2];
    subplot('Position',pos1)
    locs=time(index);
    amp=ROI_exTrace_Unbleach(index, j); %%%%% amp from Ca_Unbleach trace
    
    plot(time, dff_ex(j, :), 'color',CM(k+j,:), 'LineWidth', 0.2);
    hold on

    plot(locs,spks,'xk')
%       ylim([-0.001 0.08]);



    
    ylabel(['C' num2str(k+j)])
    if j==1
        box off
        xlabel('Time(s)')
    else
        set(gca, 'box','off','XTickLabel',[],'xtick',[])
    end
    hold on
    spks=spks';

    Cell_spks_ex(j).spk=spks;
    Cell_spks_ex(j).locs=locs;
    Cell_spks_ex(j).index=index';
    Cell_spks_ex(j).amp=amp;
    Cell_spks_ex(j).Ave_amp=mean(amp,1);
    Cell_spks_ex(j).Ave_spk=mean(spks, 1);
    Cell_spks_ex(j).Max_amp=max(amp);
    Cell_spks_ex(j).Max_spk=max(spks);
    Cell_spks_ex(j).Num_spk=size(spks, 1);
    Cell_spks_ex(j).Freq_spk=size(spks, 1)./max(time);

end
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Journal of Neuroscience Methods
% Volume 349, 1 February 2021, 109041
% synchronous spikes are identified using a threshold based on predefined percentile (mean + 1.5-2 X standard deviation, s.d.). 
% Synchronicity rate is quantified as the number of the detected synchronous spikes in one minute. 
% Quantification of single-neuron amplitude and frequency. The amplitude of peaks in each trace is defined as the mean value of △F/F0 of individual peak. 
% The frequency of peaks in each trace is defined as the number of detected peaks in one minute.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Synchronous fire rate (min to /s)
Total_trace_ex=mean(dff_ex, 1);
threshold_Peak=mean(Total_trace_ex)+NumStd*std(Total_trace_ex, 1)
[Nw_spks, Nw_index] = findpeaks(Total_trace_ex,'MinPeakHeight',threshold_Peak);
hfig609=figure(609);
plot(time, dff_ex', 'color',[0.5 0.5 0.5], 'LineWidth', 0.2);
hold on
plot(time, Total_trace_ex, '-b', 'LineWidth', 0.5);
Nw_locs=time(Nw_index);
hold on

plot([min(time) max(time)], [threshold_Peak threshold_Peak], 'color',[0.5 0.5 0.5], 'LineWidth', 0.2);
hold on
plot(Nw_locs,Nw_spks,'xk');
hold off
xlabel('Time(sec)', 'FontSize', 12);
ylabel('\Delta F/F',  'FontSize',12);
title('Synchronous Peaks')

NetworkAna_ex.spks=Nw_spks;
NetworkAna_ex.index=Nw_index;
NetworkAna_ex.trace=Total_trace;
NetworkAna_ex.Num_spks=size(Nw_spks, 2);
NetworkAna_ex.Freq_spks=size(Nw_spks, 2)./max(time);



%% run RFP analysis

disp('Open RFP file---- Captured RFP 60x -xx-20ms.nd2 file')
%%%%%%%%%%%%%%%%%%%%%%%
 cMapR = interp1([0;1],[0 0 0; 1 0 0],linspace(0,1,256));
%cMapR = 'gray';

 
RFP_lim=[200 20000];   % [200 4000]  ; 202209 [200 8000]; 202210 [200 30000]
RFP2_lim=[0.05 0.3];
RFP3_lim=[0.05 0.15]; % [-0.00 0.10]
%%%%%%%%%%%%%%%%spont

%%
[openfileRFP,openmkpath]=uigetfile('*.nd2','Please select Raw Nikon Ca files');
 
disp(openfileRFP);

fidRFP = fopen(openfileRFP, 'r');
RFPRaw= bfopen(openfileRFP); 
%%
RFPRaw=RFPRaw{1, 1}; 
RFPRaw=RFPRaw{1};
%%
RFP=rescale(RFPRaw);
AveRFP=mean2(RFP)
StdRFP=std2(RFP)
RFP_th=2
BsRFP=AveRFP+RFP_th*StdRFP
RFP_df=RFP-BsRFP;
%%%%% replace negative value (<=mean+2sd) with zero
RFP_df=max(RFP_df, 0);



%%
hfig800=figure(800);

 
clf
  
pos1=[0.05  0.50  0.45 0.45]; 
subplot('Position',pos1)

imagesc(RFPRaw, RFP_lim);
colormap(cMapR);
axis off; axis image;  
colorbar;
hold on
for k = 1 : size(ROI_in, 1)
	ROIxy_in =ROI_in{k};
	plot(ROIxy_in(:,2), ROIxy_in(:,1), 'color',CM(k,:), 'LineWidth', 2);
    text(mean(ROIxy_in(:,2)), mean(ROIxy_in(:,1)), num2str(k), 'Color','k','FontSize',8)
    hold on
    clear ROIxy
end
title('Inclusion')
hold off;

pos2=[0.05  0.01  0.45 0.45]; 
subplot('Position',pos2)
%%%%% normalize image [0 1]

imagesc(RFP, RFP2_lim);
colormap(cMapR);
axis off; axis image;  
colorbar;
hold on

for k = 1 : size(ROI_in, 1)
	ROIxy_in =ROI_in{k};
	plot(ROIxy_in(:,2), ROIxy_in(:,1), 'color',CM(k,:), 'LineWidth', 2);
    text(mean(ROIxy_in(:,2)), mean(ROIxy_in(:,1)), num2str(k), 'Color','k','FontSize',8)
    hold on
    clear ROIxy
end
hold off;
title('Inclusion (Normalized [0 1])')

%%%%%% exclusion 
pos3=[0.52  0.50  0.45 0.45]; 
subplot('Position',pos3)

imagesc(RFPRaw, RFP_lim);
colormap(cMapR);
axis off; axis image;  
colorbar;
hold on
for j = 1 : size(ROI_ex, 1)
	ROIxy_ex =ROI_ex{j};
	plot(ROIxy_ex(:,2), ROIxy_ex(:,1), 'color',CM(k+j,:), 'LineWidth', 2);
    text(mean(ROIxy_ex(:,2)), mean(ROIxy_ex(:,1)), num2str(k), 'Color','w','FontSize',8)
    hold on
    clear ROIxy
end
title('Exclusion')
hold off;

pos4=[0.52  0.01  0.45 0.45]; 
subplot('Position',pos4)
%%%%% normalize image [0 1]
 
imagesc(RFP, RFP2_lim);
colormap(cMapR);
axis off; axis image;  
colorbar;
hold on

for j = 1 : size(ROI_ex, 1)
	ROIxy_ex =ROI_ex{j};
	plot(ROIxy_ex(:,2), ROIxy_ex(:,1), 'color',CM(k+j,:), 'LineWidth', 2);
    text(mean(ROIxy_ex(:,2)), mean(ROIxy_ex(:,1)), num2str(k+j), 'Color','w','FontSize',8)
    hold on
    clear ROIxy
end
hold off;
title('Exclusion (Normalized [0 1])')

hfig801=figure(801);
clf

pos2=[0.01  0.05  0.48 0.9]; 
subplot('Position',pos2)
%%%%% normalize image [0 1]

imagesc(RFP_df, RFP3_lim);
colormap(cMapR);
axis off; axis image;  
colorbar;
hold on

for k = 1 : size(ROI_in, 1)
	ROIxy_in =ROI_in{k};
	plot(ROIxy_in(:,2), ROIxy_in(:,1), 'color',CM(k,:), 'LineWidth', 2);
    text(mean(ROIxy_in(:,2)), mean(ROIxy_in(:,1)), num2str(k), 'Color','w','FontSize',8)
    hold on
    clear ROIxy
end
hold off;
title('Inclusion (F-(mean+2SD))')

pos4=[0.51  0.05  0.48 0.9]; 
subplot('Position',pos4)
%%%%% normalize image [0 1]
 
imagesc(RFP_df, RFP3_lim);
colormap(cMapR);
axis off; axis image;  
colorbar;
hold on

for j = 1 : size(ROI_ex, 1)
	ROIxy_ex =ROI_ex{j};
	plot(ROIxy_ex(:,2), ROIxy_ex(:,1), 'color',CM(k+j,:), 'LineWidth', 2);
    text(mean(ROIxy_ex(:,2)), mean(ROIxy_ex(:,1)), num2str(k+j), 'Color','w','FontSize',8)
    hold on
    clear ROIxy
end
hold on;
title('Exclusion (F-(mean+2SD))')

%% find index inside ROI_in
sz=size(I);
CC=bwconncomp(I_in);
CC_index=CC.PixelIdxList;
for i=1:size(CC_index, 2)
   [xx yy] =ind2sub(sz, CC_index{i});
   ROI_in_index{i}=[xx yy];
   ROIxy=ROI_in_index{i};
   RFP_temp=RFP_df(ROIxy(:, 1), ROIxy(:, 2));
   RFPRawinfo{i}.in=RFP_temp;
   RFPTotalFluo_in(i)=sum(sum(RFP_temp));
   RFPMaxFluo_in(i)=max(max(RFP_temp));
end
RFPTotalFluo_in=RFPTotalFluo_in';
RFPMaxFluo_in=RFPMaxFluo_in';


%% find index inside ROI_ex
sz=size(I);
CC2=bwconncomp(I_ex);
CC_index2=CC2.PixelIdxList;
for i=1:size(CC_index2, 2)
   [xx yy] =ind2sub(sz, CC_index2{i});
   ROI_ex_index{i}=[xx yy];
   ROIxy=ROI_ex_index{i};
   RFP_temp=RFP_df(ROIxy(:, 1), ROIxy(:, 2));
   RFPRawinfo{i}.ex=RFP_temp;
   RFPTotalFluo_ex(i)=sum(sum(RFP_temp));
   RFPMaxFluo_ex(i)=max(max(RFP_temp));
end
RFPTotalFluo_ex=RFPTotalFluo_ex';
RFPMaxFluo_ex=RFPMaxFluo_ex';

% Ca2ROI=Ca2(ROIxy(:, 1), ROIxy(:, 2), n);

%%




%%%%%%%%%%%%%%%%%%%%%%

%% save files
save(File_name , 'time', 'ROITrace_Unbleach', 'CamFreq', 'Ca_lim', 'Ca2_lim', 'baseLineIndex', 'mvSpeed', 'NumROI',...
    'ROI_inTrace_Unbleach','ROI_exTrace_Unbleach');

File_name6 = ['Event_' openmkpath((length(openmkpath)-22):(length(openmkpath)-1))...
    '-' openfile1(1:(lenght_f1-4)) openROI((length(openROI)-7):(length(openROI)-4) ) '.mat'] 

File_name_fig604=['Fig_' openfile1(1:(lenght_f1-4))  '-spk' '.fig'];
File_name_fig605=['Fig_' openfile1(1:(lenght_f1-4))  '-ImSpk' '.fig'];

File_name_fig606=['Fig_' openfile1(1:(lenght_f1-4))  '-spk_in' '.fig'];
File_name_fig607=['Fig_' openfile1(1:(lenght_f1-4))  '-ImSpk_in' '.fig'];

File_name_fig608=['Fig_' openfile1(1:(lenght_f1-4))  '-spk_ex' '.fig'];
File_name_fig609=['Fig_' openfile1(1:(lenght_f1-4))  '-ImSpk_ex' '.fig'];

File_name_fig800=['Fig_' openfile1(1:(lenght_f1-4))  '-RFP' '.fig'];
File_name_fig801=['Fig_' openfile1(1:(lenght_f1-4))  '-RFP2' '.fig'];




save(File_name6 , 'time', 'dff', 'Cell_spks', 'NetworkAna',...
    'dff_in', 'Cell_spks_in', 'NetworkAna_in',...
    'dff_ex', 'Cell_spks_ex', 'NetworkAna_ex', 'RFPRawinfo', ...
    'RFPTotalFluo_in', 'RFPMaxFluo_in', 'RFPTotalFluo_ex', 'RFPMaxFluo_ex');
 
  

savefig(hfig,File_name2);  
savefig(hfig2,File_name3);
savefig(hfig301,File_name_fig301);
savefig(hfig302,File_name_fig302);
savefig(hfig401,File_name_fig401);
savefig(hfig402,File_name_fig402);

savefig(hfig604,File_name_fig604);
savefig(hfig605,File_name_fig605);
savefig(hfig606,File_name_fig606);
savefig(hfig607,File_name_fig607);
savefig(hfig608,File_name_fig608);
savefig(hfig609,File_name_fig609);


%%
savefig(hfig800,File_name_fig800);
savefig(hfig801,File_name_fig801);


% signalMatrix=a;
% [signalPeaks, signalPeaksArray, signalSigmas] = computeSignalPeaks(signalMatrix, varargin)
	% Binarize [0,1] input analog signals based on peaks in the signal.
	% Biafra Ahanonu
	% started: 2013.10.28
	% inputs
	  % signalMatrix: [nSignals frame] matrix
	% outputs
		% signalPeaks: [nSignals frame] matrix. Binary matrix with 1 = peaks.
		% signalPeaksArray: {1 nSignals} cell array. Each cell contains [1 nPeaks] vector that stores the frame locations of each peak.
	% options
		% See below.
		% % make a plot?
		% options.makePlots = 0;
		% % show waitbar?
		% options.waitbarOn = 1;
		% % make summary plots of spike information
		% options.makeSummaryPlots = 0;
		% % number of standard deviations above the threshold to count as spike
		% options.numStdsForThresh = 3;
		% % minimum number of time units between events
		% options.minTimeBtEvents = 8;
		% % shift peak detection
		% options.nFramesShift = 0;
		% % should diff and fast oopsi be done?
		% options.addedAnalysis = 0;
		% % use simulated oopsi data
		% options.oopsiSimulated = 0;
disp('Find Event Done')