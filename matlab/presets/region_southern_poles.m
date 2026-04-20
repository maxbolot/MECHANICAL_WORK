function region = region_southern_poles()
    region.name = 'southern_poles';
    region.title_label = 'Southern Poles';
    region.result_label = 'Southern Poles (90°S to 60°S, all longitudes)';
    region.average_label = 'southern poles';
    region.lat_bounds = [-90, -60];
    region.lon_bounds = [];
end