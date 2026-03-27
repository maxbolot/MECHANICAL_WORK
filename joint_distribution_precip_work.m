function[] = joint_distribution_precip_work(file_precip, file_work, varname, file_grid_area, code, results_folder, lon1, lon2, lat1, lat2)
% VARBIN is cell area weight binned by dissipation
% VARBIN_UNORM is cell area (unnormalized) binned by dissipation

    precip_edges  = logspace(-6,-1,400);
    precip_bin    = (precip_edges(1:end-1) + precip_edges(2:end)) / 2;
    nbin_precip   = length(precip_bin);
    
    work_edges   = linspace(-2500,2500,5001);
    work_bin     = (work_edges(1:end-1) + work_edges(2:end)) / 2;
    nbin_work    = length(work_bin);
    
    lon   = ncread(file_grid_area, 'grid_xt_coarse');
    ilon1 = find(lon>=lon1,1,'first');
    ilon2 = find(lon<=lon2,1,'last');
    nlon  = ilon2-ilon1+1;
    lat   = ncread(file_grid_area, 'grid_yt_coarse');
    ilat1 = find(lat>=lat1,1,'first');
    ilat2 = find(lat<=lat2,1,'last');
    nlat  = ilat2-ilat1+1; 
    time  = ncread(file_precip, 'time');
    ntime = length(time);
    
    VARBIN = zeros(nbin_precip, nbin_work);

    grid_area = ncread(file_grid_area, 'cell_area', [ilon1 ilat1], [nlon nlat]);
    
    for i=1:ntime
        
        disp(['accumulation over step: ' num2str(i)])
        
        precip = ncread(file_precip, 'PRATEsfc_coarse', [ilon1 ilat1 i], [nlon nlat 1]);
        precip(precip<=min(precip_edges)) = min(precip_edges);
        precip(precip>=max(precip_edges)) = max(precip_edges);
        precip_index = discretize(precip, precip_edges);
        
        work = ncread(file_work, varname, [ilon1 ilat1 i], [nlon nlat 1]);
        work(work<=min(work_edges)) = min(work_edges);
        work(work>=max(work_edges)) = max(work_edges);
        work_index = discretize(work, work_edges);

        VARBIN = VARBIN + accumarray([precip_index(:) work_index(:)], grid_area(:), [nbin_precip nbin_work], @sum);
        
    end
    
    data = struct(...
        'file_precip',file_precip,...
        'file_work',file_work,...
        'file_grid_area',file_grid_area,...
        'precip_edges',precip_edges,...
        'precip_bin',precip_bin,...
        'work_edges',work_edges,...
        'work_bin',work_bin,...
        'variable',varname,...
        'VARBIN',VARBIN./(ntime*sum(grid_area,[1 2])),...
        'VARBIN_UNORM',VARBIN./ntime);
    
    if ~isfolder(results_folder)
        mkdir(results_folder)
    end
    
    save([results_folder '/results_' code '.mat'],'-struct','data','-v7.3');
    
end
    
