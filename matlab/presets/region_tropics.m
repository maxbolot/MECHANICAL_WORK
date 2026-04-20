function region = region_tropics()
    region.name = 'tropics';
    region.title_label = 'Tropical';
    region.result_label = 'Tropics (30°S to 30°N)';
    region.average_label = 'tropical band';
    region.lat_bounds = [-30, 30];
    region.lon_bounds = [];
end