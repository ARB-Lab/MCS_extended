%% SLURM parpool (not relevant for computations on local machine)
function pdir = start_parallel_pool_on_SLURM_node()
tempdir = getenv('SLURM_TEMP');
if ~isempty(tempdir) && isempty(gcp('nocreate'))
    pdir = [tempdir 'prefs'];
    jdir = [tempdir 'jobs'];
    mkdir(pdir); % locate preferences-directory to tmp path
    mkdir(jdir); % locate job/worker-exchange-directory to tmp path
    [~,mem] = unix('sacct -j $SLURM_JOB_ID --format=reqmem -P -n --noconvert'); % set allocated memory
    setenv('SLURM_REQMEM',strtok(mem,'Mn'));                                    % for the job, as env-variable
    try copyfile([prefdir() '/matlabprefs.mat'],pdir), catch, end
    setenv('SLURM_JOB_STORAGE_LOCATION',jdir);
    disp(['slurm job storage location ''' jdir ''' was created and is used for parpool']);
    clstr = parcluster('local');
    clstr.NumWorkers = feature('numcores');
    clstr.JobStorageLocation = jdir;
    parpool(clstr);
    wait(parfevalOnAll(@startcna,0,1)); % startcna on all workers
end
end