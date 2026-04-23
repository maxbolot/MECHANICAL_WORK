function out = read_tropical_histogram_diagnostics(ncfile)
    try
        out.time = ncread(ncfile, 'time');
        out.pr_edges = double(ncread(ncfile, 'pr_edges'));
        out.work_edges = double(ncread(ncfile, 'work_edges'));

        % ncread reverses NetCDF dimension order (Fortran/column-major):
        %   hist_area:  NetCDF (time, nbin_pr) -> MATLAB (nbin_pr, time)
        %   hist2d_*:   NetCDF (time, nbin_work, nbin_pr) -> MATLAB (nbin_pr, nbin_work, time)
        hist_area_raw = double(ncread(ncfile, 'hist_area'));
        hist2d_lift_raw = double(ncread(ncfile, 'hist2d_lift'));
        hist2d_work_raw = double(ncread(ncfile, 'hist2d_work'));

        fprintf('time size: %d\n', numel(out.time));
        fprintf('hist_area size (raw): %s\n', mat2str(size(hist_area_raw)));
        fprintf('hist2d_lift size (raw): %s\n', mat2str(size(hist2d_lift_raw)));
        fprintf('hist2d_work size (raw): %s\n', mat2str(size(hist2d_work_raw)));

        out.hist_area = hist_area_raw';
        out.hist2d_lift = permute(hist2d_lift_raw, [3, 2, 1]);
        out.hist2d_work = permute(hist2d_work_raw, [3, 2, 1]);

        fprintf('hist_area size (permuted): %s\n', mat2str(size(out.hist_area)));
        fprintf('hist2d_lift size (permuted): %s\n', mat2str(size(out.hist2d_lift)));

        out.ntime = numel(out.time);
        out.pr_centers = sqrt(out.pr_edges(1:end-1) .* out.pr_edges(2:end));
        out.work_centers = 0.5 * (out.work_edges(1:end-1) + out.work_edges(2:end));
    catch ME
        error('Failed to read required histogram variables: %s', ME.message);
    end
end
