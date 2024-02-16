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
%%
NumROI=size(ROI, 1);
%%
for i = 1:size(ROI, 1)
    for n=1: size(Ca2, 3)
        ROIxy=ROI_index{i};
        Ca2ROI=Ca(ROIxy(:, 1), ROIxy(:, 2), n);
       T(:, n)= mean2(Ca2ROI);
       
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
    
    
  %  ROITraceRaw_Unbleach= ROITraceRaw;

    ss(:, i)=s;

     %%%%%%%%%% covnert to deltaF trace
     ROITrace_Unbleach(:, i)=ROITraceRaw_Unbleach(:, i)./...
         mean(ROITraceRaw_Unbleach(baseLineIndex, i))-1;

     
end
subplot(2,1,1)
for i = 1:size(ROI, 1)
  plot(time, ROITraceRaw(:, i), 'color',CM(i,:)); 
  hold on
end
plot(time, mean(ROITraceRaw, 2), '-k', 'LineWidth',4)
hold on
plot([min(time) max(time)], [0 0], '--k');
xlim([min(time) max(time)])
xlabel('Time(sec)', 'FontSize', 12);
ylabel('Fluo Intensity')
title('Raw Trace')
hold off

subplot(2,1,2)
for i = 1:size(ROI, 1)
  plot(time, ROITraceRaw_Unbleach(:, i), 'color',CM(i,:)); 
  hold on
end

plot(time, mean(ROITraceRaw_Unbleach, 2), '-k', 'LineWidth',4)
hold on
plot([min(time) max(time)], [0 0], '--k');
xlim([min(time) max(time)])
xlabel('Time(sec)', 'FontSize', 12);
ylabel('Fluo Intensity')
title('Raw Trace-Unbleach')
hold off

%%
hfig2=figure(300);
clf
for i = 1:size(ROI, 1)

    
    plot(time, ROITrace_Unbleach(:, i), 'color',CM(i,:)); 
    hold on

end
plot(time, mean(ROITrace_Unbleach, 2), '-k', 'LineWidth',4)
hold on
plot([min(time) max(time)], [0 0], '--k');
xlim(Ca2_plotXlim)
ylim(Ca2_plotYlim)
ytickformat('%.2f');
xlabel('Time(sec)', 'FontSize', 12);
ylabel('\Delta F/F',  'FontSize',12);
hold off

%% plot traces1 
hfig=figure(601);
clf
 
 
pos1=[0.05  0.59  0.8 0.40]; 
subplot('Position',pos1)
imagesc(Ca(:, :, 10), Ca_lim);
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
 
% pos1=[0.1  (m*0.5)/m   0.8 0.5/m];
% subplot('Position',pos1)
%     plot(time, ROITrace(1, :), 'color',CM(1,:), 'LineWidth', 2);
%     ylabel(['C' num2str(1)])
%     box off
for i = 1:NumROI
    pos1=[0.05  0.05+(i*0.48)/NumROI  0.9 0.48/NumROI];
    subplot('Position',pos1)
    plot(time, ROITrace_Unbleach(:, i), 'color',CM(i,:), 'LineWidth', 0.5);
    hold on
    %%%%% draw a line at 60s for CNO exp
     %xline(60, '--k', 'LineWidth', 0.5);
    hold off
    ylabel([num2str(i)], 'FontSize',8)
    if i==1
        box off
        xlabel('Time(s)')
        set(gca, 'box','off',  'YTickLabel',[],'ytick',[]);   

    else
        set(gca, 'box','off','XTickLabel',[],'xtick',[], 'XColor','none', 'YTickLabel',[],'ytick',[]);
    end
      xlim(Ca2_plotXlim) ;
      % xlim([0 55]) ;
        ylim(Ca2_plotYlim)
    ytickformat('%.2f');
 
end
 
%% find events

 
    [T,m] = size(ROITraceRaw) 

%%
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
 
for i = 1:NumROI
    pos1=[0.1  0.01+(i*0.9)/NumROI  0.85 0.9/NumROI];
    subplot('Position',pos1)
    plot(time, dff(i, :), 'color',CM(i,:), 'LineWidth', 0.2);
    ylabel(['C' num2str(i)], 'FontSize',8)
    if i==1
        box off
        xlabel('Time(s)')
    else
        set(gca, 'box','off','XTickLabel',[],'xtick',[])
    end
    ytickformat('%.2f');
    ylim([-0.001 0.01])
      %  xlim([0 480]) 
end
 

%%% conver ddf to spikes



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% setp 5 run correlation maps
        tic
        load('spikes.mat')
%  upsampling from 10Hz to 30Hz
        for i=1:23
            spikes{:,i}=interp(spikes{:, i}, 3);
        end
        % 
 %         addpath(genpath(toolbox)) %Add the toolbox to the matlab working directory when ever you begin a new session

         fps=30
        height=1 % background intensity
%         dff=dff';

         [m,T]=size(dff) 
       
        for ii = 1:m
            x = dff(ii,:);
         %   x = interp1(1:T,x,1:fps/30:T);
            dff1(ii,:) = x;
            x = dff1(ii,:);
            parfor i=1:length(spikes)
                snippet = spikes{i};
                L = length(snippet);
                C = zeros(size(x));       
                               
                        
                for j=1:length(x)-(L-1)
                    x_snippet = x(j:j+L-1);
                    if(range(x_snippet)>height)
                        R = corrcoef(x_snippet,snippet);
                        C(j) = R(1,2);
                        
                        if j == length(x)-(L-1)
                            for j1 = length(x)-(L-2):length(x)-round(L/2)
                                x_snippet = x(j1:end);
                                R = corrcoef(x_snippet,snippet(1:length(x_snippet)));
                                C(j1) = R(1,2);
                            end
                        end
                    end
                 end
             Call(i,:) = C;
            end
 






            Ca3{:,ii} = Call;
        end
        %%
        [m,T]=size(dff1);
        
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
    pos1=[0.08  0.01+(i*0.9)/NumROI  0.9 0.85/NumROI];
    subplot('Position',pos1)
    locs=time(index);
    amp=ROITrace_Unbleach(index, i); %%%%% amp from Ca_Unbleach trace
    
    plot(time, dff(i, :), 'color',CM(i,:), 'LineWidth', 0.2);
    
    set(gca, 'box','off','XTickLabel',[],'xtick',[])
    
    hold on

    plot(locs,spks,'xk', 'LineWidth', 0.2)
 %       ylim([-0.001 0.08]);
   ylim([-0.001 0.01])



    
    ylabel(['C' num2str(i)],'FontSize',8);
    ytickformat('%.2f');
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
    Cell_spks(i).Freq_spk=size(spks, 1)./(max(time)+T0/CamFreq);

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
hfig605=figure(605);
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
ytickformat('%.2f');


NetworkAna.spks=Nw_spks;
NetworkAna.index=Nw_index;
NetworkAna.trace=Total_trace;
NetworkAna.Num_spks=size(Nw_spks, 2);
NetworkAna.Freq_spks=size(Nw_spks, 2)./(max(time)+T0/CamFreq);


CaIm1=Ca(:, :, 10);

%% save files

save(File_name , 'CaIm1', 'ROITraceRaw', 'time', 'ROITrace_Unbleach', 'CamFreq', 'Ca_lim', 'Ca2_lim', 'baseLineIndex', 'mvSpeed', 'NumROI');
 
File_name6 = ['Event_' openmkpath((length(openmkpath)-22):(length(openmkpath)-1))...
    '-' openfile1(1:(lenght_f1-4)) openROI((length(openROI)-7):(length(openROI)-4) ) '_' fileName.timestamp '.mat'] 
save(File_name6 , 'time', 'ROITrace', 'ROI', 'ROI_index', 'ROITrace_Unbleach', 'CamFreq', 'Ca_lim', 'Ca2_lim', 'baseLineIndex', 'mvSpeed', 'NumROI', ...
    'dff', 'Cell_spks', 'NetworkAna');
 
File_name_fig604=['Fig_' openfile1(1:(lenght_f1-4))  '-spk' '.fig'];
File_name_fig605=['Fig_' openfile1(1:(lenght_f1-4))  '-NtSpk' '.fig'];

savefig(hfig,File_name2);  
savefig(hfig2,File_name3);

savefig(hfig604,File_name_fig604);
savefig(hfig605,File_name_fig605);


%%
File_name2_1 = ['Fig_' openfile1(1:(lenght_f1-4))  '-1' '.eps'] 

exportgraphics(hfig,File_name2_1 ,'BackgroundColor','none','ContentType','vector');  

%
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
