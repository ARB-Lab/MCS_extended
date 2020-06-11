% This script compares the new geneMCSEnumerator2 with the existing gMCS 
% approach by Apaolaza et al. 2017 (https://doi.org/10.1038/s41467-017-00555-y)
% by computing synthetic lethals up to the size of 4 gene deletions in the genome-scale
% E. coli model iML1515. At the end of the script, the computation time and the 
% number of solutions is returned for both cases.
%
% -Jun 2020
%

% start CNA
if ~exist('cnan','var')
    startcna(1)
end
% Add helper functions to matlab path
function_path1 = [fileparts(mfilename('fullpath') ) '/../functions'];
function_path2 = [fileparts(mfilename('fullpath') ) '/../e_coli'];
addpath(function_path1);
addpath(function_path2);

% start Cobra Toolbox
global CBTDIR
if isempty(CBTDIR)
    initCobraToolbox(0);
end
clear('CBTDIR');

% Start Parallel pool
if ~isempty(getenv('SLURM_JOB_ID')) && isempty(gcp('nocreate'))
    % start parpool and locate preferences-directory to tmp path
    prefdir = start_parallel_pool_on_SLURM_node();
% If running on local machine, start parallel pool and keep compression
% flags as defined above.
elseif license('test','Distrib_Computing_Toolbox') && isempty(getCurrentTask()) && ...
       (~isempty(ver('parallel'))  || ~isempty(ver('distcomp'))) && isempty(gcp('nocreate'))
    parpool(); % remove this line if MATLAB Parallel Toolbox is not available
    wait(parfevalOnAll(@startcna,0,1)); % startcna on all workers
end

load('iML1515.mat')

% cnap = CNAsetGenericReactionData_with_array(cnap,'geneProductAssociation',cellstr(strcat('gene_',cnap.reacID)));
cnap.reacMin(cnap.reacMin>=0) = 0;
cnap.reacMin(cnap.reacMin<0)  = -inf;
cnap.reacMax(cnap.reacMax<=0) =  0;
cnap.reacMax(cnap.reacMax>0) =  inf;

[~,~,~,gpr_Rules] = CNAgenerateGPRrules(cnap);
maxSize = 4;
% Compare results from gmcs (Apaolaza 2018) with the gene-extension MCS

%% new MCS
idx_bm = find(cnap.objFunc);
T = {sparse(1,idx_bm,-1,1,cnap.numr)};
t = {-1e-3}; % this is the default threshold used by Apaolaza
maxSolutions = inf;
options.mcs_search_mode     = 2; % bottom-up stepwise enumeration of MCS.
options.preproc_check_feas  = false;
options.preproc_D_violations= 0;
options.postproc_verify_mcs = false;
tic;
[rmcs, new_mcs, new_cnap, cmp_mcs, cmp_cnap, mcs_idx_cmp_full] = CNAgeneMCSEnumerator2(cnap,T,t,{},{},[],[],maxSolutions,maxSize,[],[],gpr_Rules,options,1);
time_new = toc;

%% Apaolaza
model = CNAcna2cobra(cnap);
model.grRules = CNAgetGenericReactionData_as_array(cnap,'geneProductAssociation');

options.timelimit = inf;
% options.forceLength = 0;

model = generateRules(model,0);
model = buildRxnGeneMat(model);
global CBT_MILP_SOLVER;
CBT_MILP_SOLVER = 'ibm_cplex';
delete(['G_' getenv('SLURM_JOB_ID') '_.mat']);
tic;
[gmcs, gmcs_time] = calculateGeneMCS('', model, inf, maxSize, options);
time_apa = toc;
delete(['G_' getenv('SLURM_JOB_ID') '_.mat']);

apa_gmcs = text2num_mcs(gmcs,new_cnap);

%% validate and compare

[valid_mcs,max_mue_mcs] = verify_lethals(new_cnap,new_mcs);
[valid_apa,max_mue_apa] = verify_lethals(new_cnap,apa_gmcs);

global time_apa_milp;
global time_new_milp;
disp([{'new'} {'apa'}; num2cell([time_new_milp time_apa_milp ; time_new time_apa; length(valid_mcs) length(valid_apa)])]);

[~,~,~,compare_mat] = compare_mcs_sets(new_mcs,apa_gmcs);
if all(sum(compare_mat,1) == 3) && all(sum(compare_mat,2) == 3)
    disp('The solutions found by both approaches are identical.')
else
    disp('The solutions found by both approaches are not identical.')
end
rmpath(function_path1);
rmpath(function_path2);
