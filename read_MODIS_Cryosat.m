clear all; close all; clc;


[filename, pathname] = uigetfile('*.hdf');

modis_raw = hdfread(filename, ...
    '/MODIS_SWATH_Type_L1B/Data Fields/EV_250_RefSB', ...
    'Index', {[1  1  1],[1  1  1],[2  8120  5416]});
lat = hdfread(filename, '/MODIS_SWATH_Type_L1B/Geolocation Fields/Latitude',...
    'Index', {[1  1],[1  1],[2030  1354]});lat = double(lat);
lon = hdfread(filename, '/MODIS_SWATH_Type_L1B/Geolocation Fields/Longitude',...
    'Index', {[1  1],[1  1],[2030  1354]});lon = double(lon);


band1=double(squeeze(modis_raw(1,:,:)));


lat = double(lat);
lon = double(lon);

lat1 = imresize(lat,[8120 5416]);
lon1 = imresize(lon, [8120 5416]);


%% 라인그리고 값까지 추출하는거
[HDR, CS]=Cryo_L2_read('CS_OFFL_SIR_SAR_2A_20140411T134439_20140411T135025_B001.DBL');
% cslon = CS.MEA.LON_20Hz(:); % 이것을 먼저 하시오!
% 
last = 7609; % 140411 스발바르오른쪽섬

 cslon = CS.MEA.LON_20Hz(:); cslon(cslon == 0) = nan; cslon=cslon(isfinite(cslon));
 cslat = CS.MEA.LAT_20Hz(:); cslat(cslat == 0) = nan; cslat=cslat(isfinite(cslat));
 len = length(cslon);
 
 bks = CS.MEA.backsc_sig_20Hz; bks = bks(1:last); 
 
[HDR1, CS1]=Cryo_L1b_read('CS_OFFL_SIR_SAR_1B_20140411T134439_20140411T135025_B001.DBL'); 

 std = CS1.SAR.beam_param(1,:,:); std = std(1:last); %std(std == 0) = nan;
 
 kur = CS1.SAR.beam_param(5,:,:); kur = kur(1:last); 
 
 skew = CS1.SAR.beam_param(4,:,:); skew = skew(1:last); %skew(skew == 0) = nan;
%  skew = skew(isfinite(skew));
 
 wf = CS1.SAR.data; [a, b, c] = size(wf);
 wf  = reshape(wf,a,b*c);
%  wf(wf == 0) = nan; wf = wf(isfinite(wf));

for i = 1:len
    line = wf(:,i);
    pp(i) = (max(line)*1)./sum(line);
%     max_bin = find(line == max(line)); max_bin = max_bin(1);
%     pp_r(i) = (max(line*3))./mean([line(max_bin-3),line(max_bin-1)]);
%     pp_l(i) = (max(line*3))./mean([line(max_bin+1),line(max_bin+3)]);
end
 for x = 1:7609
     if kur(x) > 650
         kur(x) = kur(x) - 655.36;
     else
        kur(x) = kur(x);
     end
 end
 
  for x = 1:7609
     if skew(x) > 650
         skew(x) = skew(x) - 655.36;
     else
        skew(x) = skew(x);
     end
  end

kk = [std' skew' kur' pp' bks'];


 %% 라인 분류에 따라 색깔 넣기
hold on;

type = zeros(7755,1); type(:,1) = NaN;

%% Machine learning criteria (New lead detection) : "Total_April_L2I_add_ocean_rf"
load('RF_prediction_0924')
% load('aa')
lead_pos = find(RFprediction0924==1);
ice_pos = find(RFprediction0924==0);    
ocean_pos = find(RFprediction0924==2);    

% m_proj('lambert','long',[60 69],'lat',[76 81.3]); %with ocean
m_proj('lambert','long',[65 69],'lat',[80 81.15]); %without ocean
m_pcolor(lon1,lat1,band1); shading flat
m_grid('box','fancy','tickdir','in','fontsize',35,'fontname','arial');
m_gshhs_i('patch',[.6 .6 .6]);
colormap('bone');
set(gca,'clim',[0 3000]);
hold on;



m_scatter(cslon(ice_pos),cslat(ice_pos),15,type(ice_pos),'markerfacecolor',[0/255 179/255 255/255]); % Blue=ice
m_scatter(cslon(lead_pos),cslat(lead_pos),15,type(lead_pos),'markerfacecolor',[255/255 0/255 0/255]); % Orange=Lead
m_scatter(cslon(ocean_pos),cslat(ocean_pos),15,type(ocean_pos),'markerfacecolor',[11/255 75/255 178/255]); % Orange=Lead

set(gca,'clim',[0 2000]);
 set(gcf,'paperposition',[0.01 0.01 45 29]);
print('-dtiff','-r300','RF_result_250');


