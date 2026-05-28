function data_root = get_globalfv3_data_root()
    % See matlab/README.md for environment setup and override examples.
    % Resolve the GLOBALFV3 data root from env var or OS-specific defaults.
    env_root = getenv('GLOBALFV3_DATA_ROOT');
    if ~isempty(env_root)
        data_root = env_root;
        return;
    end

    if ispc
        % Windows local/offline analysis.
        data_root = 'C:\climate_processed_data\GLOBALFV3';
    elseif isunix
        % Princeton HPC analysis environment.
        data_root = '/scratch/gpfs/mbolot/results/GLOBALFV3';
    else
        error('Unknown operating system. Set GLOBALFV3_DATA_ROOT explicitly.');
    end
end