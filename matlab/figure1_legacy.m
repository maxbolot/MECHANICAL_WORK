warning off all

% data location
folder_map_work_control = 'D:/Projects/PIRE.2020/PIRE/mechanical_work/work';
folder_map_water_lift_control = 'D:/Projects/PIRE.2020/PIRE/mechanical_work/water_lift';
folder_map_dissip_control = 'D:/Projects/PIRE.2020/PIRE/frictional_dissipation/dissipation_total';
folder_map_prmax_control = 'D:/Projects/PIRE.2020/PIRE/maps/prmax';

folder_map_work_warming = 'D:/Projects/PIRE.2020/PIRE.4KCO2/mechanical_work/work';
folder_map_water_lift_warming = 'D:/Projects/PIRE.2020/PIRE.4KCO2/mechanical_work/water_lift';
folder_map_dissip_warming = 'D:/Projects/PIRE.2020/PIRE.4KCO2/frictional_dissipation/dissipation_total';
folder_map_prmax_warming = 'D:/Projects/PIRE.2020/PIRE.4KCO2/maps/prmax';

%% COMPUTE AVERAGE OVER ZONES

% zone 1: global

lift_control_zone1 = ncread(fullfile(folder_map_water_lift_control,"lift.2020.m90p90.nc"), 'water_lift');
work_control_zone1 = ncread(fullfile(folder_map_work_control,"work.2020.m90p90.nc"), 'work');
dissip_control_zone1 = ncread(fullfile(folder_map_dissip_control,"dissip.2020.m90p90.nc"), 'dissip');
prmax_control_zone1 = squeeze(ncread(fullfile(folder_map_prmax_control,'prmax.2020.m90p90.nc'), 'pr'));
ke_control_zone1 = work_control_zone1 - lift_control_zone1;
scale_height_control_zone1 = dissip_control_zone1 ./ prmax_control_zone1 / 9.81 / 1000;     % in km

lift_warming_zone1 = ncread(fullfile(folder_map_water_lift_warming,"lift.2020.m90p90.nc"), 'water_lift');
work_warming_zone1 = ncread(fullfile(folder_map_work_warming,"work.2020.m90p90.nc"), 'work');
dissip_warming_zone1 = ncread(fullfile(folder_map_dissip_warming,"dissip.2020.m90p90.nc"), 'dissip');
prmax_warming_zone1 = squeeze(ncread(fullfile(folder_map_prmax_warming,'prmax.2020.m90p90.nc'), 'pr'));
ke_warming_zone1 = work_warming_zone1 - lift_warming_zone1;
scale_height_warming_zone1 = dissip_warming_zone1 ./ prmax_warming_zone1 / 9.81 / 1000;     % in km

change_prmax_zone1 = (prmax_warming_zone1 - prmax_control_zone1) ./ prmax_control_zone1 / 4 * 100;
change_scale_height_zone1 = (scale_height_warming_zone1 - scale_height_control_zone1) ./ scale_height_control_zone1 / 4 * 100;


% zone 2: tropics -30 to +30

lift_control_zone2 = ncread(fullfile(folder_map_water_lift_control,"lift.2020.m30p30.nc"), 'water_lift');
work_control_zone2 = ncread(fullfile(folder_map_work_control,"work.2020.m30p30.nc"), 'work');
dissip_control_zone2 = ncread(fullfile(folder_map_dissip_control,"dissip.2020.m30p30.nc"), 'dissip');
prmax_control_zone2 = squeeze(ncread(fullfile(folder_map_prmax_control,'prmax.2020.m30p30.nc'), 'pr'));
ke_control_zone2 = work_control_zone2 - lift_control_zone2;
scale_height_control_zone2 = dissip_control_zone2 ./ prmax_control_zone2 / 9.81 / 1000;     % in km

lift_warming_zone2 = ncread(fullfile(folder_map_water_lift_warming,"lift.2020.m30p30.nc"), 'water_lift');
work_warming_zone2 = ncread(fullfile(folder_map_work_warming,"work.2020.m30p30.nc"), 'work');
dissip_warming_zone2 = ncread(fullfile(folder_map_dissip_warming,"dissip.2020.m30p30.nc"), 'dissip');
prmax_warming_zone2 = squeeze(ncread(fullfile(folder_map_prmax_warming,'prmax.2020.m30p30.nc'), 'pr'));
ke_warming_zone2 = work_warming_zone2 - lift_warming_zone2;
scale_height_warming_zone2 = dissip_warming_zone2 ./ prmax_warming_zone2 / 9.81 / 1000;     % in km

change_prmax_zone2 = (prmax_warming_zone2 - prmax_control_zone2) ./ prmax_control_zone2 / 4 * 100;
change_scale_height_zone2 = (scale_height_warming_zone2 - scale_height_control_zone2) ./ scale_height_control_zone2 / 4 * 100;

% zone 3: maritime continent

lift_control_zone3 = ncread(fullfile(folder_map_water_lift_control,"lift.2020.mc.nc"), 'water_lift');
work_control_zone3 = ncread(fullfile(folder_map_work_control,"work.2020.mc.nc"), 'work');
dissip_control_zone3 = ncread(fullfile(folder_map_dissip_control,"dissip.2020.mc.nc"), 'dissip');
prmax_control_zone3 = squeeze(ncread(fullfile(folder_map_prmax_control,'prmax.2020.mc.nc'), 'pr'));
ke_control_zone3 = work_control_zone3 - lift_control_zone3;
scale_height_control_zone3 = dissip_control_zone3 ./ prmax_control_zone3 / 9.81 / 1000;     % in km

lift_warming_zone3 = ncread(fullfile(folder_map_water_lift_warming,"lift.2020.mc.nc"), 'water_lift');
work_warming_zone3 = ncread(fullfile(folder_map_work_warming,"work.2020.mc.nc"), 'work');
dissip_warming_zone3 = ncread(fullfile(folder_map_dissip_warming,"dissip.2020.mc.nc"), 'dissip');
prmax_warming_zone3 = squeeze(ncread(fullfile(folder_map_prmax_warming,'prmax.2020.mc.nc'), 'pr'));
ke_warming_zone3 = work_warming_zone3 - lift_warming_zone3;
scale_height_warming_zone3 = dissip_warming_zone3 ./ prmax_warming_zone3 / 9.81 / 1000;     % in km

change_prmax_zone3 = (prmax_warming_zone3 - prmax_control_zone3) ./ prmax_control_zone3 / 4 * 100;
change_scale_height_zone3 = (scale_height_warming_zone3 - scale_height_control_zone3) ./ scale_height_control_zone3 / 4 * 100;


% zone 4: NH +30 to +60

lift_control_zone4 = ncread(fullfile(folder_map_water_lift_control,"lift.2020.p30p60.nc"), 'water_lift');
work_control_zone4 = ncread(fullfile(folder_map_work_control,"work.2020.p30p60.nc"), 'work');
dissip_control_zone4 = ncread(fullfile(folder_map_dissip_control,"dissip.2020.p30p60.nc"), 'dissip');
prmax_control_zone4 = squeeze(ncread(fullfile(folder_map_prmax_control,'prmax.2020.p30p60.nc'), 'pr'));
ke_control_zone4 = work_control_zone4 - lift_control_zone4;
scale_height_control_zone4 = dissip_control_zone4 ./ prmax_control_zone4 / 9.81 / 1000;     % in km

lift_warming_zone4 = ncread(fullfile(folder_map_water_lift_warming,"lift.2020.p30p60.nc"), 'water_lift');
work_warming_zone4 = ncread(fullfile(folder_map_work_warming,"work.2020.p30p60.nc"), 'work');
dissip_warming_zone4 = ncread(fullfile(folder_map_dissip_warming,"dissip.2020.p30p60.nc"), 'dissip');
prmax_warming_zone4 = squeeze(ncread(fullfile(folder_map_prmax_warming,'prmax.2020.p30p60.nc'), 'pr'));
ke_warming_zone4 = work_warming_zone4 - lift_warming_zone4;
scale_height_warming_zone4 = dissip_warming_zone4 ./ prmax_warming_zone4 / 9.81 / 1000;     % in km

change_prmax_zone4 = (prmax_warming_zone4 - prmax_control_zone4) ./ prmax_control_zone4 / 4 * 100;
change_scale_height_zone4 = (scale_height_warming_zone4 - scale_height_control_zone4) ./ scale_height_control_zone4 / 4 * 100;

%% PLOT

Y = [lift_control_zone1 ke_control_zone1; ...
    lift_warming_zone1 ke_warming_zone1; ...
    lift_control_zone2 ke_control_zone2; ...
    lift_warming_zone2 ke_warming_zone2; ...
    lift_control_zone3 ke_control_zone3; ...
    lift_warming_zone3 ke_warming_zone3; ...
    lift_control_zone4 ke_control_zone4; ...
    lift_warming_zone4 ke_warming_zone4];

newY = reshape([reshape(Y,2,[]); zeros(1,numel(Y)/2)],[],2);  % Add zeroes for spacing

Y2 = [change_prmax_zone1  change_scale_height_zone1; ...
    change_prmax_zone2  change_scale_height_zone2; ...
    change_prmax_zone3  change_scale_height_zone3; ...
    change_prmax_zone4  change_scale_height_zone4];


fig = figure('units','inch','position',[0,0,8,3.8]);
set(gcf,'color','w')

t = tiledlayout(1,2);

ax1 = nexttile([1 1]);
set(ax1,'FontSize',13)
b1=bar(newY,'stacked','FaceColor','flat');
set(gca,'XLim',[0 12],'XTick',[1.5 4.5 7.5 10.5],'XTickLabel',{'global', 'tropics', 'maritime continent', 'north. midlatitudes'});
ylabel('W m^{-2}')
legend('control','warming','location','northwest')
legend('boxoff')
title('a')

% set colors
clr1 = [rgb('blue');...
    rgb('red');...
    rgb('white');...
    rgb('blue');...
    rgb('red');...
    rgb('white');...
    rgb('blue');...
    rgb('red');...
    rgb('white');...
    rgb('blue');...
    rgb('red');...
    rgb('white');];
clr2 = [rgb('pale blue');...
    rgb('pale red');...
    rgb('white');...
    rgb('pale blue');...
    rgb('pale red');...
    rgb('white');...
    rgb('pale blue');...
    rgb('pale red');...
    rgb('white');...
    rgb('pale blue');...
    rgb('pale red');...
    rgb('white');];
b1(1).CData = clr1;
b1(2).CData = clr2;


ax2 = nexttile([1 1]);
set(ax2,'FontSize',13)
b2=bar(Y2,'stacked','FaceColor','flat');
set(gca,'XLim',[0.5 4.5],'XTick',[1 2 3 4],'XTickLabel',{'global', 'tropics', 'maritime continent', 'north. midlatitudes'});
ylabel('% K^{-1}')
ylim([0 9])
% legend('control','warming','location','northwest')
% legend('boxoff')
title('b')
text(0.05,0.9,'precip. dissipation change','Units','normalized')

% set colors
clr1 = [rgb('grey');...
    rgb('grey');...
    rgb('grey');...
    rgb('grey');];
clr2 = [rgb('very pale grey');...
    rgb('very pale grey');...
    rgb('very pale grey');...
    rgb('very pale grey');];
b2(1).CData = clr1;
b2(2).CData = clr2;
b2(1).BarWidth = 0.26;


t.Padding = 'compact';
t.TileSpacing = 'compact';
ax1.TitleHorizontalAlignment = 'left';
ax2.TitleHorizontalAlignment = 'left';


% export_fig ../figures_paper.v6/figure1.eps -painters



