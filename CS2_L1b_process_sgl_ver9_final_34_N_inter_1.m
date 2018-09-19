clear all;close all;clc

list = dir('CS_OFFL_SIR_SAR_1B*.DBL');
list1 = dir('CS_OFFL_SIR_SAR_2A*.DBL');

[filename, pathname] = uigetfile('*.hdf');

modis_raw = hdfread(filename, ...
    '/MODIS_SWATH_Type_L1B/Data Fields/EV_250_RefSB', ...
    'Index', {[1  1  1],[1  1  1],[2  8120  5416]});
lat = hdfread(filename, '/MODIS_SWATH_Type_L1B/Geolocation Fields/Latitude',...
    'Index', {[1  1],[1  1],[2030  1354]});lat = double(lat);
lon = hdfread(filename, '/MODIS_SWATH_Type_L1B/Geolocation Fields/Longitude',...
    'Index', {[1  1],[1  1],[2030  1354]});lon = double(lon);
band1=double(squeeze(modis_raw(1,:,:)));
band2=band1(1:4:end,1:4:end);

lat = double(lat);
lon = double(lon);

lat1 = imresize(lat,[8120 5416]);
lon1 = imresize(lon, [8120 5416]);

figure; pcolor(lon1,lat1,band1); shading flat
colormap('bone');
set(gca,'clim',[0 2500]);
hold on
scatter(cs2_lon,cs2_lat,'.');


for list_l=1%:1:length(list)
    
[HDR,CS] = Cryo_L1b_read(list(list_l).name);
[HDR2,CS2] =Cryo_L2_read(list1(list_l).name); 
 cs2_sar_wf = CS.SAR.data; 
 cs2_sar_lon = CS.GEO.LON;
 cs2_sar_lat = CS.GEO.LAT;
 cs2_sar_alt = CS.GEO.H;
 cs2_sar_win_delay = CS.MEA.win_delay;
 
 [bin_l,freq_l,trk_l] = size(cs2_sar_wf);
 
 surf_elev = zeros(freq_l,trk_l);surf_elev(:,:) = nan;
 tl_pos = zeros(freq_l,trk_l);tl_pos(:,:) = nan;
 pp = zeros(freq_l,trk_l);pp(:,:) = nan;
 pp_r = zeros(freq_l,trk_l);pp_r(:,:) = nan;
 pp_l = zeros(freq_l,trk_l);pp_l(:,:) = nan;

 
 
 for i=1:1:freq_l
 clear cs2_wf cs2_lon cs2_lat cs2_alt cs2_win_delay 
 cs2_wf = squeeze(cs2_sar_wf(:,i,:));
 cs2_lon = squeeze(cs2_sar_lon(i,:));
 cs2_lat = squeeze(cs2_sar_lat(i,:));
 cs2_alt = squeeze(cs2_sar_alt(i,:));
 cs2_win_delay = squeeze(cs2_sar_win_delay(i,:));
  
 cs2_lon(cs2_lon == 0) = nan;
 cs2_lat(cs2_lat == 0) = nan;
 cs2_alt(cs2_alt == 0) = nan;
 cs2_win_delay(cs2_win_delay == 0) = nan;
% figure(1)
% imagesc(cs2_wf)
 
%============================================================
%           Lowpass filtering of Wave form
%============================================================
 
 for j=1:1:trk_l-1

 line = cs2_wf(:,j);
 
  
 if line(:) == 0
 continue
 else

filtlow = 20;
[f1,f2] = butter(1,filtlow/(length(line)/2),'low'); %40 and type
fline = filtfilt(f1,f2,line);
fline = fline./max(fline);
 
%=========================================================================
%                          First Maximum bin
%=========================================================================
 
 fv10 = find(fline > 0.1);
 devs_fline = gradient(gradient(fline(fv10(1):end)));
 
 idx_dev = zeros(length(devs_fline),1);
 for dev_fv=1:1:length(devs_fline)-2;
     fwdv = devs_fline(dev_fv) - devs_fline(dev_fv+1);
     backv = devs_fline(dev_fv+2) - devs_fline(dev_fv+1);
     
     if fwdv < 0 && backv > 0
         idx_dev(dev_fv) = 1;
     else 
         idx_dev(dev_fv) = 0;
     end
 end
 fv0 = find(idx_dev == 0);
 idx_dev(1:fv0(1))=0;
 maxbin = find(idx_dev == 1)+(fv10(1)-1);
 
 if isempty(maxbin) || maxbin(1) < 40
    continue
 else
 k= 1;
 pp(i,j) = (max(line)*k)/sum(line); 
%  ppr_min = maxbin(1) - 3; ppr_max = maxbin(1) - 1;
%  ppl_min = maxbin(1) + 1; ppl_max = maxbin(1) + 3;
%  pp_r(i,j) = fline(maxbin(1))*3./mean([fline(ppr_min),fline(ppr_max)]);
%  pp_l(i,j) = fline(maxbin(1))*3./mean([fline(ppl_min),fline(ppl_max)]);
 
 %=========================================================================
 %                      Determine DC level
 %=========================================================================
 
 %dc_samp = fline(10:maxbin(1));
 %dc_samp_diff = dc_samp - min(dc_samp);
 %dc = mean(dc_samp(dc_samp_diff < max(dc_samp_diff)*0.01));
 dc_samp = fline(10:15);
 dc = mean(dc_samp);

 %==============================================================
 %                Determine  Threshold level
 %==============================================================
 
 t_coeff = 0.4;
 tl = dc + t_coeff*(fline(maxbin(1)) - dc); tl = tl(1);
 
 p1 = find(fline > tl);
 
 if isnan(dc) == 1 || p1(1) == 1
     tl_pos(i,j) = nan;
 else
 tl_pos(i,j) = (p1(1)-1) + (tl - fline(p1(1)-1))/(fline(p1(1))-fline(p1(1)-1));
 end
%  
%  figure(2);set(gcf,'w');set(gca,'fontsize',30,'fontweight','bold','linewidth',1.5);hold on
%  plot(1:length(fline),fline,'b-','linewidth',2);
%  plot(1:length(fline),dc,'k-','linewidth',1.5);
%  plot(1:length(fline),fline(maxbin(1)),'k','linewidth',1.5);
%  plot(1:length(fline),tl,'g-','linewidth',5);
%  plot([tl_pos(i) tl_pos(i)],[0 1],'r-','linewidth',1.5);
%  plot([64 64],[0 1],'k-','linewidth',1.5);
%  set(gca,'xlim',[0 128],'ylim',[0 1]);grid on
%  box on
%  ylabel('Nomalized power','fontsize',35,'fontname','georgia');
%  xlabel('Range bine','fontsize',35,'fontname','georgia');
%  
%  
% set(gcf,'paperposition',[0.1 0.1 35 25]);
% print('-dtiff','-r600',['wf' num2str(i)]);
 end
 end

%==========================================================================
%                  surface elevation from ellipsoid
%==========================================================================

c = 299792458; %[m/s]
rb_resol = 0.2342; %[m]

rg1 = (0.5*c.*cs2_win_delay);

rg2 = rb_resol.*(tl_pos(i,:) - length(fline)/2);

rg3 = (CS.COR.dry_trop+CS.COR.wet_trop+CS.COR.inv_bar+CS.COR.dac+...
     CS.COR.gim_ion+CS.COR.model_ion+CS.COR.ocean_equilibrium_tide+...
     CS.COR.ocean_longperiod_tide+CS.COR.ocean_loading_tide+...
     CS.COR.solidearth_tide+CS.COR.geocentric_polar_tide);

 
surf_elev(i,:) = cs2_alt - (rg1+rg2+rg3);
var_mean_surf_elev = nanvar(surf_elev);
var_mean_surf_elev1 = repmat(var_mean_surf_elev,20,1);
var_mean_surf_elev1 = var_mean_surf_elev1(:);
a_surf_elev = surf_elev(:);
a_surf_elev(var_mean_surf_elev1 >= 0.5) = nan;
a_surf_elev(a_surf_elev >= 70) = nan;
a_surf_elev = a_surf_elev';

% test = smooth(a_surf_elev,10);

% figure(3)
% plot(surf_elev(i,:));
end
 end
 
%==========================================================================
%                 lead detection [Ricker, 2014]
%==========================================================================

a_pp = pp(:); a_pp(isnan(a_surf_elev) ==1) = nan;

std = squeeze(CS.SAR.beam_param(1,:,:)); std = std(:);
kur = squeeze(CS.SAR.beam_param(5,:,:)); kur = kur(:);
skew = squeeze(CS.SAR.beam_param(4,:,:)); skew = skew(:);
std(isnan(a_surf_elev) ==1) = nan;
kur(isnan(a_surf_elev) ==1) = nan;
skew(isnan(a_surf_elev) ==1) = nan;
msig = CS2.MEA.backsc_sig_20Hz(:);
msig(isnan(a_surf_elev) ==1) = nan;

for x= 1:freq_l*trk_l
        if kur(x) > 650
            kur(x) = kur(x)-655.36;
          else 
            kur(x) = kur(x);
         end
      
end


    for x= 1:freq_l*trk_l
        if skew(x) > 650
            skew(x) = skew(x)-655.36;
          else 
            skew(x) = skew(x);
         end
    end


type = zeros(1,freq_l*trk_l);

for l = 1:freq_l*trk_l

%Ocean
if skew(l) <= 0.73 && msig(l) > 14.89 && pp(l) <= 0.043
    type(l) = 3;
 else if skew(l) <= 0.73 && msig(l) > 14.89 && pp(l) > 0.043 && msig(l) > 16.48;
   type(l) = 3;
   else if  skew(l) > 0.73 && msig(l) > 27.8 && kur(l) <= 17.53 && std(l) > 25.6;
    type(l) = 3;
       end
     end
end

%Lead 
if skew(l) > 0.73 && msig(l) > 27.8 && kur(l) > 17.53
    type(l) = 1;
 else if skew(l) > 0.73 && msig(l) > 27.8 && kur(l) <= 17.53 && std(l) <= 25.6 && msig(l) > 31.47
    type(l) = 1;
     end
end    
%Ice
if skew(l) <= 0.73 && msig(l) <= 14.89
    type(l) =2;
else if skew(l) <= 0.73 && msig(l) > 14.89 && pp(l) > 0.0431 && msig(l) <= 16.48
        type(l) = 2;
  else if skew(l) > 0.73 && msig(l) <= 27.8
         type(l) = 2;
      else if  skew(l) > 0.73 && msig(l) > 27.8 && kur(l) <= 17.53 && std(l) <= 25.6 && msig(l) <= 31.47;
              type(l) = 2;
          end  
      end
end
end
    
end

load dtu10_arctic

 cs2_sar_lon(cs2_sar_lon == 0) = nan; cs2_lon = cs2_sar_lon(:);
cs2_sar_lat(cs2_sar_lat == 0) = nan;  cs2_lat = cs2_sar_lat(:);


R = georasterref('RasterInterpretation','postings');
R.RasterSize = [1401 10801];
mxlon = max(max(dtu_lon));
mxlat = max(max(dtu_lat));
mlon = min(min(dtu_lon));
mlat = min(min(dtu_lat));
R.Latlim = [mlat mxlat];
R.Lonlim = [mlon mxlon];
   
mssh = ltln2val(dtu_mssh,R,cs2_lat,cs2_lon,'cubic');
mssh = mssh';
a_surf_elev(isnan(mssh) == 1) = nan;


lead_x = find(type==1);

% if isempty(lead_x) == 1
% ssha = [];
% end

if length(lead_x) == 1 %리드가 하나도 없을 때
ssha = a_surf_elev(type == 1) - mssh(type ==1);
end
if length(lead_x) == 3
    ssha = a_surf_elev(lead_x) - mssh(lead_x);
   pre_ssha = mean(ssha);
   radar_fb = a_surf_elev - (mssh + pre_ssha);
 cs2_sar_lon(cs2_sar_lon == 0) = nan; cs2_lon = cs2_sar_lon(:);
cs2_sar_lat(cs2_sar_lat == 0) = nan;  cs2_lat = cs2_sar_lat(:);
else
if length(lead_x) == 1
    pre_ssha = a_surf_elev(lead_x) - mssh(lead_x);
    radar_fb = a_surf_elev - (mssh+pre_ssha);
 cs2_sar_lon(cs2_sar_lon == 0) = nan; cs2_lon = cs2_sar_lon(:);
cs2_sar_lat(cs2_sar_lat == 0) = nan;  cs2_lat = cs2_sar_lat(:); 
else
    if length(lead_x) == 2
   ssha = a_surf_elev(lead_x) - mssh(lead_x);
   pre_ssha = mean(ssha);
   radar_fb = a_surf_elev - (mssh + pre_ssha);
 cs2_sar_lon(cs2_sar_lon == 0) = nan; cs2_lon = cs2_sar_lon(:);
cs2_sar_lat(cs2_sar_lat == 0) = nan;  cs2_lat = cs2_sar_lat(:);
else

if isempty(lead_x) == 1
    ssha = [];
    radar_fb = a_surf_elev - (mssh+nan);
 cs2_sar_lon(cs2_sar_lon == 0) = nan; cs2_lon = cs2_sar_lon(:);
cs2_sar_lat(cs2_sar_lat == 0) = nan;  cs2_lat = cs2_sar_lat(:);
else
if length(lead_x) >= 2
lead_ssha1 = zeros(size(a_surf_elev));
lead_ssha1(:) = nan;
lead_ssha1(type == 1) = ((a_surf_elev(type == 1)) - (mssh(type == 1)));
inter_x = find(type == 1);
inter_v = lead_ssha1(inter_x);
inter_xq = 1:length(lead_ssha1);
inter = interp1(inter_x,inter_v,inter_xq,'linear');

inter_find = isfinite(inter);
inter_find1 = find(inter_find == 1);
inter_idx1 = inter_find1(1);
inter_idx_end = inter_find1(end);
inter1 = inter(inter_idx1:inter_idx_end);
inter_length = length(inter); inter_length1 = 1:inter_length;

 inter(1:inter_idx1) = inter(inter_idx1);
 inter(inter_idx_end:end) = inter(inter_idx_end);

end
ssha = inter;
filt_idx =  length(type)*0.9;
filtlow1 = filt_idx;

[f11,f22] = butter(1,filtlow1/(length(type)),'low'); %40 and type
fssha = filtfilt(f11,f22,ssha);

% cs2_sar_lon(cs2_sar_lon == 0) = nan;
% cs2_sar_lat(cs2_sar_lat == 0) = nan;
% lon = nanmean(cs2_sar_lon);
% lat = nanmean(cs2_sar_lat);
% lon = lon(inter_idx1:inter_idx_end);
% lat = lat(inter_idx1:inter_idx_end);


% idx = ~isnan(lead_ssha1);
% lead_ssha2 = lead_ssha1(idx);
% lead_pos = cumsum(idx);
% fv1 = find(lead_pos == 0);
% lead_pos(fv1) = 1;
% fv2 = find(lead_pos == length(lead_ssha2)+1);
% if isempty(fv2) == 0
% lead_pos(fv2) = length(lead_ssha2);
% ssha = lead_ssha2(lead_pos);
% end
% mssh(isnan(mean_surf_elev)==1)=nan;
if isempty(ssha) == 1
 radar_fb = nan;
else
radar_fb = a_surf_elev - (mssh+fssha);
fv_final = find(isnan(radar_fb)==1);
% radar_fb = smooth(radar_fb,3);
radar_fb(fv_final) = nan;
radar_fb(type == 3) = nan;
radar_fb(type == 0) = nan;
radar_fb(radar_fb < 0) = nan; %더 높아짐
% mean_surf_elev = mean_surf_elev(inter_idx1:inter_idx_end);
end 
end
end
end
end
% radar_fb(type < 2) = nan;





radar_fb(type == 3) = nan;
radar_fb(type == 0) = nan;
radar_fb(radar_fb < 0) = nan; %더 높아짐
radar_fb(radar_fb > 1) = nan;

radar_fb = smooth(radar_fb,35);
radar_fb(radar_fb < 0) = nan; %더 높아짐
radar_fb(radar_fb > 1) = nan;
fv_final = find(isnan(radar_fb)==1);
radar_fb(fv_final) = nan;


% % Anomalies
% figure(1)%;set(gcf,'pos',[-1800 80 1800 600],'color','w');
% subplot(131)
% imagesc(cs2_wf)
% ylabel('Range bins','fontsize',25)
% xlabel('Numner of tracks','fontsize',25)
% 
% % figure(2)
% subplot(132)
% plot(mean_surf_elev,'linewidth',2)
% hold on
% plot(mssh,'r','linewidth',2)
% hold on
% plot(ssha,'k','linewidth',2)
% ylabel('Meter [m]','fontsize',25)
% xlabel('Numner of tracks','fontsize',25)
% title(char(list(list_l).name(1:55)),'fontsize',35)
% legend('Mean.Surf.Elev','Mssh','Ssha','fontsize',40)
% 
% % figure(3)
% subplot(133)
% plot(radar_fb,'r','linewidth',2)
% ylabel('Meter [m]','fontsize',25)
% xlabel('Numner of tracks','fontsize',25)
% legend('Fb-Smoothing','fontsize',40)
% cd plot2
% set(gcf,'paperposition',[0.1 0.1 40 10]);
% print('-dtiff','-r300',char(list(list_l).name(1:55)));
% cd ..
% close all
cd nc_34
filename = list(list_l).name;filename = filename(1:end-4);
nccreate([filename '.nc'],'long','Dimensions',{'x' length(cs2_lon)},'Format','classic');
nccreate([filename '.nc'],'lat','Dimensions',{'x' length(cs2_lon)});
nccreate([filename '.nc'],'surf','Dimensions',{'x' length(cs2_lon)});
nccreate([filename '.nc'],'fb','Dimensions',{'x' length(cs2_lon)});
nccreate([filename '.nc'],'type','Dimensions',{'x' length(cs2_lon)});

ncwriteatt([filename '.nc'],'long','_FillValue', 9999);
ncwriteatt([filename '.nc'],'long','long_name','longitude');
ncwriteatt([filename '.nc'],'long','standard_name','lon');
ncwriteatt([filename '.nc'],'long','units','degrees_east');
ncwriteatt([filename '.nc'],'long','valid_max','180');
ncwriteatt([filename '.nc'],'long','valid_min','-180');

ncwriteatt([filename '.nc'],'lat','_FillValue', 9999);
ncwriteatt([filename '.nc'],'lat','long_name','latitude');
ncwriteatt([filename '.nc'],'lat','standard_name','lat');
ncwriteatt([filename '.nc'],'lat','units','degrees_north');
ncwriteatt([filename '.nc'],'lat','valid_min','-90');
ncwriteatt([filename '.nc'],'lat','valid_max','90');

ncwriteatt([filename '.nc'],'surf','_FillValue', 9999);
ncwriteatt([filename '.nc'],'surf','long_name','surface elevation');
ncwriteatt([filename '.nc'],'surf','standard_name','h');
ncwriteatt([filename '.nc'],'surf','units','m');
ncwriteatt([filename '.nc'],'surf','valid_max',max(a_surf_elev));
ncwriteatt([filename '.nc'],'surf','valid_min',min(a_surf_elev));

ncwriteatt([filename '.nc'],'fb','_FillValue', 9999);
ncwriteatt([filename '.nc'],'fb','long_name','radar freeboard');
ncwriteatt([filename '.nc'],'fb','standard_name','fb');
ncwriteatt([filename '.nc'],'fb','units','m');
ncwriteatt([filename '.nc'],'fb','valid_max',max(radar_fb));
ncwriteatt([filename '.nc'],'fb','valid_min',min(radar_fb));

ncwrite([filename '.nc'],'long',cs2_lon);
ncwrite([filename '.nc'],'lat',cs2_lat);
ncwrite([filename '.nc'],'surf',a_surf_elev);
ncwrite([filename '.nc'],'fb',radar_fb);
ncwrite([filename '.nc'],'type',type);
cd ..
end

