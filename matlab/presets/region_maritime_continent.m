function region = region_maritime_continent()
    region.name = 'maritime_continent';
    region.title_label = 'Maritime Continent';
    region.result_label = 'Maritime Continent (15°S to 15°N, 90°E to 150°E)';
    region.average_label = 'maritime continent';
    region.lat_bounds = [-15, 15];
    region.lon_bounds = [90, 150];
end