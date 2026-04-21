warning off all

% data location 
folder_map_dissip_control = 'D:/Projects/PIRE.2020/PIRE/frictional_dissipation/dissipation_total';
folder_map_prmax_control = 'D:/Projects/PIRE.2020/PIRE/maps/prmax';

folder_map_dissip_warming = 'D:/Projects/PIRE.2020/PIRE.4KCO2/frictional_dissipation/dissipation_total';
folder_map_prmax_warming = 'D:/Projects/PIRE.2020/PIRE.4KCO2/maps/prmax';

%% DATA
% (keep only annual, skip JJA and DJF)

% control

% annual
p1.dissip = squeeze(ncread(fullfile(folder_map_dissip_control,'dissip.2020.nc'), 'dissip'));
p1.prmax = squeeze(ncread(fullfile(folder_map_prmax_control,'prmax.annual.nc'), 'pr'));
p1.lon = ncread(fullfile(folder_map_dissip_control,'dissip.2020.nc'), 'lon');
p1.lat = ncread(fullfile(folder_map_dissip_control,'dissip.2020.nc'), 'lat');
p1.scale_height = p1.dissip ./ p1.prmax / 9.81 / 1000;     % in km

% warming

% annual
p4.dissip = squeeze(ncread(fullfile(folder_map_dissip_warming,'dissip.2020.nc'), 'dissip'));
p4.prmax = squeeze(ncread(fullfile(folder_map_prmax_warming,'prmax.annual.nc'), 'pr'));
p4.lon = ncread(fullfile(folder_map_dissip_warming,'dissip.2020.nc'), 'lon');
p4.lat = ncread(fullfile(folder_map_dissip_warming,'dissip.2020.nc'), 'lat');
p4.scale_height = p4.dissip ./ p4.prmax / 9.81 / 1000;     % in km


%% PLOT

fig = figure('units','inch','position',[0,0,12,6]);
set(gcf,'color','w')

t = tiledlayout(3,1);


% dissipation 

ax1 = nexttile([1 1]);
[X,Y] = ndgrid(p1.lon,p1.lat);
axesm('mercator', 'Frame', 'on', 'Grid', 'on', 'MapLatLimit', [-60 60],...  
    'meridianlabel', 'off', 'parallellabel', 'on', 'labelformat', 'signed', ...
    'mlabellocation', 60, 'plinelocation', 20, 'mlinelocation', 45, 'flinewidth', 1,...
    'MLabelParallel','south','FontSize',10)
geoshow(Y,X,p1.dissip,'DisplayType','texturemap')
shading flat
geoshow('landareas.shp','FaceColor','None')
axis off
set(gca,'Color','none')
set(gca,'layer','top')
colormap(gca,'turbo')
clim([0 15])
title('a. precip. dissipation')
ax1.TitleHorizontalAlignment = 'left';
% text(0.13,0.82,'a. D_p','Units','normalized','Color','w','FontSize',12,'FontWeight','bold')
cb1 = colorbar('eastoutside','Ticks',[0 5 10 15]);
cb1.FontSize = 16;
cb1.FontWeight = 'bold';
cb1.Label.String = 'W m^{-2}';
cb1.Label.FontSize = 16;
cb1.Label.FontWeight = 'bold';
% cb1.Color = 'w';


%  precipitation

ax2 = nexttile([1 1]);
[X,Y] = ndgrid(p1.lon,p1.lat);
axesm('mercator', 'Frame', 'on', 'Grid', 'on', 'MapLatLimit', [-60 60],...  
    'meridianlabel', 'off', 'parallellabel', 'on', 'labelformat', 'signed', ...
    'mlabellocation', 60, 'plinelocation', 20, 'mlinelocation', 45, 'flinewidth', 1,...
    'MLabelParallel','south','FontSize',10)
geoshow(Y,X,p1.prmax,'DisplayType','texturemap')
shading flat
geoshow('landareas.shp','FaceColor','None')
axis off
set(gca,'Color','none')
set(gca,'layer','top')
colormap(gca,'turbo')
clim([0 20]*1e-5)
title('b. precip.')
ax2.TitleHorizontalAlignment = 'left';
% text(0.13,0.82,'b. pr_{surf}','Units','normalized','Color','w','FontSize',12,'FontWeight','bold')
cb2 = colorbar('eastoutside','Ticks',[0 5 10 15 20]*1e-5);
cb2.FontSize = 16;
cb2.FontWeight = 'bold';
cb2.Label.String = 'kg m^{-2} s^{-1}';
cb2.Label.FontSize = 16;
cb2.Label.FontWeight = 'bold';
% cb2.Color = 'w';


% scale height

ax5 = nexttile([1 1]);
[X,Y] = ndgrid(p1.lon,p1.lat);
axesm('mercator', 'Frame', 'on', 'Grid', 'on', 'MapLatLimit', [-60 60],...  
    'meridianlabel', 'on', 'parallellabel', 'on', 'labelformat', 'signed', ...
    'mlabellocation', 45, 'plinelocation', 20, 'mlinelocation', 45, 'flinewidth', 1,...
    'MLabelParallel','south','FontSize',10)
geoshow(Y,X,p1.scale_height,'DisplayType','texturemap')
shading flat
geoshow('landareas.shp','FaceColor','None')
axis off
set(gca,'Color','none')
set(gca,'layer','top')
colormap(gca,'turbo')
clim([0 10])
title('c. precip. scale height')
ax5.TitleHorizontalAlignment = 'left';
% text(0.13,0.82,'d. H_p','Units','normalized','Color','k','FontSize',12,'FontWeight','bold')
cb5 = colorbar('eastoutside','Ticks',[0 2 4 6 8 10]);
cb5.FontSize = 16;
cb5.FontWeight = 'bold';
cb5.Label.String = 'km';
cb5.Label.FontSize = 16;
cb5.Label.FontWeight = 'bold';



% 

t.Padding = 'compact';
t.TileSpacing = 'compact';

% export_fig ../figures_paper.v6/figure2.v6.eps -painters





