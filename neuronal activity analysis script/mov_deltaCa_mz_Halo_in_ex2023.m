% function mov_deltaCa_mz
% open Raw Nikon file captured by Nikon NIS-Elements
% use  a Nikon ND2 file reader. https://www.mathworks.com/matlabcentral/fileexchange/71345-nd2read
% by Mingrui Zhao of Gan Lab  6/2021
clear all

cMapG = interp1([0;1],[0 0 0; 0 1 0],linspace(0,1,256));

%% 
%%%%%
%%%%%
%%%%%% set yes or no for photobleaching
Type_photobeach=0; 
%%% 0= no need for photobleaching correction
%%% 1= yes for photobeaching correction
Type_Recording=2; 
%%% 1= Spont Recording
%%% 2= KCl Recording
%%% 3= AP stimulation 


%%
Ca_lim=[1500 10000]   % [200 4000]  ; 202209 [200 8000]; 202210 [200 30000]
 
Ca2_lim=[-0.00 0.30]  % [-0.00 0.10]
%%%%%%%%%%%%%%%%spont
switch Type_Recording 
    case 1
    CamFreq=30    %CamFreq=30
    Ca2_plotYlim=[-0.04 0.06]; % for kcl [-1 2.5]
    Ca2_plotXlim=[0 120]
%%%% for kCL
    case 2
    CamFreq=10 
    Ca2_plotXlim=[0 400]
    Ca2_plotYlim=[-0.1 1] 
%%%% for filed stimulation
    case 3
    CamFreq=30 
    Ca2_plotXlim=[0 10]
    Ca2_plotYlim=[-0.04 0.1] 
end


 
%baseLineIndex=[300 309] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% baseLineIndex=[5:10];  

baseLineIndex=[5 14] 
T0=1

mvSpeed=10; % set up mov skip frame to speed up procession.   
%%
%


[openfile1,openmkpath]=uigetfile('*.nd2','Please select Raw Nikon Ca files');
cd(openmkpath)
disp(openfile1);
disp(['1) Open ROI.tif file ..... 2) Open ROI-in.tif file(inclusion)......', ...
    '3) Open ROI-ex.tif file (exclusion)...........']);
%%
[openROI,openmkpath]=uigetfile('*.tif','Please select Nikon ROI files');
disp(openROI)
%%
[openROI_in,openmkpath]=uigetfile('*.tif','Please select Nikon ROI inclusion files');
disp(openROI_in)
[openROI_ex,openmkpath]=uigetfile('*.tif','Please select Nikon ROI exclusion files');
disp(openROI_ex)
% [~,index_ROI_in]=ismember(openROI,openROI_in)
% [~,index_ROI_ex]=ismember(openROI,openROI_ex)
%%




% finfo=nd2finfo(openfile1); 
tic
% disp(['analyzing file structure used ', sprintf('%0.2f', toc), ' seconds'])

fid = fopen(openfile1, 'r');
 
 
 
   CaRaw= bfopen(openfile1); 
% open nd2 file as a cell array
% returns an n-by-4 cell array, where n is the number of series in the dataset. If s is the series index between 1 and n:
% 
% The data{s, 1} element is an m-by-2 cell array, where m is the number of planes in the s-th series. If t is the plane index between 1 and m:
% The data{s, 1}{t, 1} element contains the pixel data for the t-th plane in the s-th series.
% The data{s, 1}{t, 2} element contains the label for the t-th plane in the s-th series.
% The data{s, 2} element contains original metadata key/value pairs that apply to the s-th series.
% The data{s, 3} element contains color lookup tables for each plane in the s-th series.
% The data{s, 4} element contains a standardized OME metadata structure, which is the same regardless of the input file format, and contains common metadata values such as physical pixel sizes - see OME metadata below for examples.   
%    
  
   CaRaw=CaRaw{1,1};  % get the data set
   %whos CaRaw
   for i=1: size(CaRaw, 1) %12000 % 6000  % just see first 4 min size(CaRaw, 1) %finfo.img_seq_count
       temp_plane=CaRaw{i,1}; % get ith plane
       Ca(:, :, i)=temp_plane;
   end
%%
 








T0=1

mvSpeed=10; % set up mov skip frame to speed up procession.   
%%
figure(100)
clf
%whos CaRaw Ca
imagesc(Ca(:, :, 1), Ca_lim)
colormap('gray');
axis image; axis off;
colorbar;

%%
%%%%%%%%%%%%%%%%%%%%%
fileName.timestamp = datestr(now, 'yyyy-mm-dd');

lenght_f1=length(openfile1);
% Mov_name=['mov_' openfile1(1:(lenght_f1-4)) '_' fileName.timestamp '.avi']
% Mov_name = ['mp4' openfile1(1:(lenght_f1-4)) '_' fileName.timestamp  '.mp4'];
Mov_name = ['mp4_' openfile1(1:(lenght_f1-4))  '.mp4'] 
%%%%
File_name = ['Trace_' openmkpath((length(openmkpath)-27):(length(openmkpath)-1)) '-' openfile1(1:(lenght_f1-4))  '.mat'] 
File_name2 = ['Fig_' openfile1(1:(lenght_f1-4))  '-1' '.fig'] 

File_name3 = ['Fig_' openfile1(1:(lenght_f1-4))  '-2' '.fig'] 
%%%%%%%%
File_name_fig301=['Fig_' openfile1(1:(lenght_f1-4))  '-in1' '.fig'];
File_name_fig302=['Fig_' openfile1(1:(lenght_f1-4))  '-ex1' '.fig'];
%%%%%
File_name_fig401=['Fig_' openfile1(1:(lenght_f1-4))  '-in2' '.fig'];
File_name_fig402=['Fig_' openfile1(1:(lenght_f1-4))  '-ex2' '.fig'];

%%%%
I=imread(openROI);
% ROI=bwboundaries(I);
[ROI, ROIL, NumROI, ROIa] =bwboundaries(I);
%% ROI_region=bwboundaries
CM = jet(size(ROI, 1));
NumROI=size(ROI, 1)

for i=1:size(ROI,1)
   [xx yy] =find(ROIL==i); 
   ROI_index{i}=[yy xx];  % ???? draw image x vs y
end
 NumROI=size(ROI, 1);



%% read ROI_in ROI_ex
I_in=imread(openROI_in);
I_ex=imread(openROI_ex);
% ROI_in=bwboundaries(I_in);
% ROI_ex=bwboundaries(I_ex);
[ROI_in, ROIL_in, NumROI_in, ROIa_in] =bwboundaries(I_in);
for k=1:size(ROI_in,1)
   [xx yy] =find(ROIL_in==k); 
   ROI_in_index{k}=[yy xx];  % ???? draw image x vs y
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[ROI_ex, ROIL_ex, NumROI_ex, ROIa_ex] =bwboundaries(I_ex);
 
for j=1:size(ROI_ex,1)
   [xx yy] =find(ROIL_ex==j); 
   ROI_ex_index{j}=[yy xx];  % ???? draw image x vs y
end
 

%%

Ca = imfilter(Ca,ones(3,3)/9);

%%
time=((1:size(Ca, 3))-T0)./CamFreq;
time=time';

%%  photon bleach correction
switch Type_photobeach 
    case 0 %%%% no need for correction
         Ca_Unbleach=double(Ca);
    case 1 %%% do correction
        [s, Ca_Unbleach]=BleachCorr(time, Ca, Ca_lim);
end
  
% Ca=Ca_Unbleach;

% Time=(1:finfo.img_seq_count).*0.040-0.040; 

%%%%%%%%%me=(1:finfo.img_seq_count).*0.040-0.040; 

%%%%%%%%%
%Ba_Ca = squeeze(mean(Ca(:,:,baseLineIndex),3));
Ba_Ca = squeeze(mean(Ca_Unbleach(:,:,baseLineIndex),3));
%  Ca2 = Ca1./(repmat(Ba_Ca,[1 1 l]));
Ca2 = Ca_Unbleach./Ba_Ca-1;
 
 

%%  whole ROI
hfig2=figure(300);
clf
%  sz=size(I);
% CC=bwconncomp(I);
% CC_index=CC.PixelIdxList;
% for i=1:size(CC_index, 2)
%    [xx yy] =ind2sub(sz, CC_index{i});
%    ROI_index{i}=[xx yy];
% end

for i = 1:size(ROI_index, 2)
    for n=1: size(Ca2, 3)
        ROIxy=ROI_index{i};
        Ca2ROI=Ca2(ROIxy(:, 1), ROIxy(:, 2), n);
       T(:, n)= mean2(Ca2ROI);
       
    end
    
    plot(time, T, 'color',CM(i,:)); 
    hold on
    ROITrace(:, i) = T;
     
end
 
plot(time, mean(ROITrace, 2), '-k', 'LineWidth',4)
hold on
plot([min(time) max(time)], [0 0], '--k');
xlim([min(time) max(time)])
ylim(Ca2_plotYlim)
%ylim()
xlabel('Time(sec)', 'FontSize', 12);
ylabel('\Delta F/F',  'FontSize',12);
hold off

%% plot traces1 
hfig=figure(400);
clf
 
pos1=[0.1  0.6  0.8 0.38]; 
subplot('Position',pos1)
imagesc(Ca(:, :, 10), Ca_lim);
colormap('gray');
axis off; axis image;  
colorbar;
hold on
for k = 1 : size(ROI_index, 2)
	ROIxy =ROI{k};
	plot(ROIxy(:,2), ROIxy(:,1), 'color',CM(k,:), 'LineWidth', 2);
    text(mean(ROIxy(:,2)), mean(ROIxy(:,1)), num2str(k), 'Color','k','FontSize',8)
 
    hold on
    clear ROIxy
end
hold off;
 
for i = 1:size(ROI_index, 2)
   pos1=[0.1  0.04+(i*0.48)/NumROI  0.8 0.48/NumROI];
 %   pos1=[0.1  0.04+(i*0.85)/NumROI  0.8 0.85/NumROI];

    subplot('Position',pos1)
    plot(time, ROITrace(:, i), 'color',CM(i,:), 'LineWidth', 1);
    ylabel(['C' num2str(i)], 'FontSize', 6)
    if i==1
        box off
        xlabel('Time(s)')
    else
        set(gca, 'box','off','XTickLabel',[],'xtick',[])
    end
      xlim(Ca2_plotXlim)  
      % ylim(Ca2_plotYlim)
end
%%%%%%%%%%
%%




%%%%%
%%  whole ROI_in
hfig301=figure(301);
clf
%%%%%%%%% index inside ROI_in
% CC=bwconncomp(I_in);
% CC_index=CC.PixelIdxList;
% for k=1:size(CC_index, 2)
%    [xx yy] =ind2sub(sz, CC_index{k});
%    ROI_in_index{k}=[xx yy];
% end


for k = 1:size(ROI_in_index, 2)
    for n=1: size(Ca2, 3)
        ROIxy_in=ROI_in_index{k};
        Ca2ROI_in=Ca2(ROIxy_in(:, 1), ROIxy_in(:, 2), n);
       T_in(n)= mean2(Ca2ROI_in);
       
    end
    
    plot(time, T_in, 'color',CM(k,:), 'LineWidth', 1); 
    hold on
    ROITrace_in(:, k) = T_in;
     
end
plot(time, mean(ROITrace_in, 2), '-k', 'LineWidth',4)
hold on
plot([min(time) max(time)], [0 0], '--k');
xlim(Ca2_plotXlim)
%ylim(Ca2_plotYlim)
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
for k = 1 : size(ROI_in_index, 2)
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
    plot(time, ROITrace_in(:, k), 'color',CM(k,:), 'LineWidth', 1);
    
    
    ylabel(['C' num2str(k)])
   if k==1
        box off
        xlabel('Time(s)')
    else
        set(gca, 'box','off','XTickLabel',[],'xtick',[])
    end
      % xlim([0 480]) 
     % ylim(Ca2_plotYlim)
end





%%   ROI_ex
hfig302=figure(302);
clf
%%%%%%%% index inside ROI_ex
% CC=bwconncomp(I_ex);
% CC_index=CC.PixelIdxList;
% for j=1:size(CC_index, 2)
%    [xx yy] =ind2sub(sz, CC_index{j});
%    ROI_ex_index{j}=[xx yy];
% end 

for j = 1:size(ROI_ex_index, 2)
    for n=1: size(Ca2, 3)
        ROIxy_ex=ROI_ex_index{j};
        Ca2ROI_ex=Ca2(ROIxy_ex(:, 1), ROIxy_ex(:, 2), n);
       T_ex(n)= mean2(Ca2ROI_ex);
       
    end
    
    plot(time, T_ex, 'color',CM(k+j,:)); 
    hold on
    ROITrace_ex(:, j) = T_ex;
     
end
plot(time, mean(ROITrace_ex, 2), '-k', 'LineWidth',4)
hold on
plot([min(time) max(time)], [0 0], '--k');
xlim([min(time) max(time)])
ylim(Ca2_plotYlim)
xlabel('Time(sec)', 'FontSize', 12);
ylabel('\Delta F/F',  'FontSize',12);
hold off

%% plot traces1 exclusion
hfig402=figure(402);
clf
 
pos1=[0.1  0.57  0.8 0.45]; 
subplot('Position',pos1)
imagesc(Ca(:, :, 10), Ca_lim);
colormap('gray');
axis off; axis image;  
colorbar;
hold on
for j = 1 : size(ROI_ex_index, 2)
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
for j = 1:size(ROI_ex_index, 2)
   pos1=[0.1  0.04+(j*0.5)/NumROI  0.8 0.5/NumROI];
 %   pos1=[0.1  0.04+(i*0.85)/NumROI  0.8 0.85/NumROI];

    subplot('Position',pos1)
    plot(time, ROITrace_ex(:, j), 'color',CM(k+j,:), 'LineWidth', 1);
    ylabel(['C' num2str(k+j)])
    if j==1
        box off
        xlabel('Time(s)')
    else
        set(gca, 'box','off','XTickLabel',[],'xtick',[])
    end
     %  xlim([0 120]) 
      ylim(Ca2_plotYlim)
end




%%







%save(File_name , 'time', 'ROITrace', 'CamFreq', 'Ca_lim', 'Ca2_lim', 'baseLineIndex', 'mvSpeed', 'NumROI',...
%    'ROITrace_ex','ROITrace_in');
% savefig(hfig,File_name2);  
% savefig(hfig2,File_name3);
% %%
% savefig(hfig301,File_name_fig301);
% savefig(hfig302,File_name_fig302);
% savefig(hfig401,File_name_fig401);
% savefig(hfig402,File_name_fig402);
%%

%%%%% verion fig 2021B
% savefig(File_name2, hfig);
% 
% savefig(File_name3, hfig2);



%%

%% movie maker
% 
% lenght_f1=length(openfile1);
% tic
% % writerObj = VideoWriter(Mov_name,'Uncompressed AVI');
% writerObj = VideoWriter(Mov_name,'MPEG-4');
% writerObj.FrameRate = CamFreq;
% % writerObj.Quality=100;
% open(writerObj);
%  
% h = figure(10);
% clf
% % colormap('jet');
% % set(h,'position',[100 100 826 364]);
% set(h,'position',[100 100 220 300]);
% % for i = Fr_mark_line(1):Fr_mark_line(2)  %1697  %
% for i = 1:mvSpeed:size(Ca2, 3)   % 10:6000 %1100:10:2400    % =4min at 10Hz   size(Ca2, 3)  %1697  %
%     
%     clf(h);
% % h1 = axes('Parent',h,'position',[10/826  182/364  262/826  152/364]);
% 
% 
% h1 = axes('Parent',h,'position',[2/220 160/300  220/220 120/300]);
% 
% imagesc(Ca(:,:,i), Ca_lim);
% axis image; axis off;
% colorbar;
% colormap(gca, cMapG)
% title('\fontsize{14}Ca Raw')
% 
% % plot(x, LFP1, 'k','LineWidth',0.5);
% % xlim(Ax_lim);
% % % xlim([-2 50]);
% % hold on
% % plot(Time_T1,T1,'k','LineWidth',0.5);
% % 
% % % plot(Time_T1,T2,'G','LineWidth',1);
% % % plot(Time_T1,T3,'R','LineWidth',1);
% % % plot(Time_T1,T4,'B','LineWidth',1);
% % 
% % xlim(Ax_lim);axis off
% % ht1 = text(5,1.8,'LFP','FontSize',10);
% % set(ht1,'Rotation',90)
% % ht2 = text(5,0.5,'iGlu','FontSize',10);
% % set(ht2,'Rotation',90)
% % % text(550,0.5,'Hbo','FontSize',12);
% % % text(550,-0.4,'Hbt','FontSize',12);
% % % text(550,-1.5,'Hbr','FontSize',12);
% % plot([i*0.016794 i*0.016794],[-0.5 2],':r')
% % title(['\fontsize{10}' num2str((fix((i-Fr_sz_onset)*0.016794*1000))./1000) 's']);  % 120HZ to 60 Hz per chan around 16.7ms per frame
% 
% 
% % h2 = axes('Parent',h,'position',[10/826 10/364 262/826  152/364]);
% %h2 = axes('Parent',h,'position',[10/252  10/300  242/252 152/300]);
% 
% h2 = axes('Parent',h,'position',[2/220  10/300   220/220 120/300]);
% 
% % imagesc(Ca2(:,:,i).*MaskImage.*(CoefImage.^3).*1.2,[-0.02 0.02]);  % use bone mask and blood vasual mask
% % imagesc(Ca2(:,:,i).*MaskImage,Glu_lim);
% imagesc(Ca2(:,:,i), Ca2_lim);
% axis off; axis image
% 
% 
% colormap(gca, cMapG)
% colorbar;
%  
% title(['\fontsize{14}\Delta F/F' ': '  num2str(time(i), '%.3f') ' sec' ])
%  %%%%% 35Hz=0.028571435; 50Hz=0.02
%  
% % h3 = axes('Parent',h,'position',[10/826  10/364  262/826  152/364]);
% % imagesc(Hbo2(:,:,i).*MaskImage,Hbx_lim);
% % axis off;
% % colorbar;
% % title('\fontsize{14}dHbo')
% % % imagesc(Ca3(:,:,i),[-0.8 0.8]);
% % % axis off;
% % % colorbar;
% % % title('Ca Normalize')
% % 
% % h4 = axes('Parent',h,'position',[282/826 10/364  262/826  152/364]);
% % imagesc(Hbt2(:,:,i).*MaskImage,Hbx_lim);
% % axis off;
% % colorbar;
% % title('\fontsize{14}dHbt')
% % 
% % h5 = axes('Parent',h,'position',[554/826 10/364 262/826  152/364]);
% % imagesc(Hbr2(:,:,i).*MaskImage,Hbx_lim);
% % axis off;
% % colorbar;
% % title('\fontsize{14}dHbr')
% % h5 = axes('Parent',h,'position',[310/600  20/600  270/600  170/600]);
% 
% 
% % pause(0.005)
% % end
%          frame =getframe(h);
%     writeVideo(writerObj,frame);
% end
% close(writerObj);
% %clearvars h h1 h2 h3 h4 h5 i i l w x b1 b2 b3 b4 
%  toc
%  
%  
%  
%  
% 
% %fclose(fid);

%%
 disp('Done')