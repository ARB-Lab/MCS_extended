% This script reproduces the results from S5 Table, Sheet P. putida, scenario 2.
%
% MCS computation for the strain design of a 2,3-butanediol production host
%
% We enumerate all Minimal Gene Cut Sets up to the size of 7 for the
% strongly growth coupled production of 2,3-butanediol from Glucose with 
% P. putida. Compared to the first scenario ("benchmark"), a second
% desired region is introduced to make sure all strain designs support higher 
% ATP demands (18 mM/gBDW/h). The enumerated MCS for scenario 2 
% (for 2 desired, 1 target) must therefore represent a subset of the MCS 
% found in the scenario 1 (1 desired, 1 Target). GPR rule compression and 
% network compression are enabled to speed up the MCS computation.
%
% % required files/models:
%   or iJN746.mat 
%
% % key variables:
%   max_num_interv: defines the maximum number of possible gene cuts
%
% % process:
%   0) Start parallel pool to speed up FVAs
%   1) Setup model, add heterologous  reactions
%   2) Define Target and Desired regions for MCS computation
%   3) Run MCS computation
%   4) Validate MCS
%   5) Characterize and Rank results
%
% Correspondence: cellnetanalyzer@mpi-magdeburg.mpg.de
% -Mar 2020
%

%% 0) Starting CNA and Parallel pool (for faster FVA), defining compression setting
if ~exist('cnan','var')
    startcna(1)
end
% Add helper functions to matlab path
function_path = [fileparts(mfilename('fullpath') ) '/../functions'];
addpath(function_path);

max_num_interv             = 6;
options.mcs_search_mode    = 2; % bottom-up stepwise enumeration of MCS.
max_solutions              = inf;
options.milp_time_limit    = inf;

% If runnning on a system with a SLURM workload manager:
% Use directory on internal memory to share data between the workers. 
if ~isempty(getenv('SLURM_JOB_ID')) && isempty(gcp('nocreate'))
    % start parpool and locate preferences-directory to tmp path
    prefdir = start_parallel_pool_on_SLURM_node();
% If running on local machine, start parallel pool if available
elseif license('test','Distrib_Computing_Toolbox') && isempty(getCurrentTask()) && ...
       (~isempty(ver('parallel'))  || ~isempty(ver('distcomp'))) && isempty(gcp('nocreate'))
    parpool();
    wait(parfevalOnAll(@startcna,0,1)); % startcna on all workers
end

%% 1) Model setup
% load model from file and prepare it for MCS computation
model = 'iJN746';
load('iJN746.mat')
cnap = CNAcobra2cna(iJN746);
reac_in_core_metabolism = true(cnap.numr,1);

cnap.reacMin(cnap.reacMin <= -1000) = -1000;
cnap.reacMax(cnap.reacMax >=  1000) =  1000;
cnap.specNotes = strcat('[', iJN746.metFormulas, ']')';

DM_SK_reacs = find(~cellfun(@isempty,regexp(cellstr(cnap.reacID),'(DM_|SK_)')));
DM_SK_min = cnap.reacMin(DM_SK_reacs);
DM__SK_max = cnap.reacMax(DM_SK_reacs);
cnap = block_non_standard_products(cnap);
cnap.reacMin(DM_SK_reacs(DM_SK_min<0)) = -0.5;
cnap.reacMin(DM_SK_reacs(DM_SK_min==0)) = 0;
cnap.reacMax(DM_SK_reacs) = 0.5;
cnap.reacMin(ismember(cnap.reacID,{'EX_glc__D_e'})) = -10;
cnap.reacMax(ismember(cnap.reacID,{'EX_gly_e'})) = 1000;
cnap.reacMax(ismember(cnap.reacID,{'EX_glyclt_e'})) = 1000;
cnap.reacMin(logical(cnap.objFunc)) = 0; % make lower bound for biomass reaction zero
cnap = CNAsetGenericReactionData_with_array(cnap,'geneProductAssociation',iJN746.grRules);

cnap = CNAaddReactionMFN(cnap,'EX_o2_anaer_e','1 o2_e = ' ,-0.2,0,0,nan,0,...
'//START_GENERIC_DATA{;:;deltaGR_0;#;num;#;NaN;:;uncertGR_0;#;num;#;NaN;:;geneProductAssociation;#;str;#;;:;}//END_GENERIC_DATA',0,0,0,0);
cnap = CNAaddReactionMFN(cnap,'ATPM','atp_c + 1 h2o_c = 1 h_c + 1 adp_c + 1 pi_c' ,1,1000,0,nan,0,...
'//START_GENERIC_DATA{;:;deltaGR_0;#;num;#;NaN;:;uncertGR_0;#;num;#;NaN;:;geneProductAssociation;#;str;#;;:;}//END_GENERIC_DATA',0,0,0,0);

% Add pathway from DOI 10.1186/s12934-018-1038-0 Erian, Pfluegl 2018
cnap = CNAaddSpeciesMFN(cnap,'actn_c',0,'3-Hydroxybutan-2-one');
cnap = CNAaddSpeciesMFN(cnap,'23bdo_c',0,'3-Hydroxybutan-2-one');
cnap = CNAaddReactionMFN(cnap,'ACLDC','1 alac__S_c + 1 h_c = 1 co2_c + 1 actn_c' ,0,1000,0,nan,0,...
'//START_GENERIC_DATA{;:;deltaGR_0;#;num;#;NaN;:;uncertGR_0;#;num;#;NaN;:;geneProductAssociation;#;str;#;;:;}//END_GENERIC_DATA',0,0,0,0);
cnap = CNAaddReactionMFN(cnap,'AR','1 h_c + 1 nadh_c + 1 actn_c = 1 nad_c + 1 23bdo_c' ,0,1000,0,nan,0,...
'//START_GENERIC_DATA{;:;deltaGR_0;#;num;#;NaN;:;uncertGR_0;#;num;#;NaN;:;geneProductAssociation;#;str;#;;:;}//END_GENERIC_DATA',0,0,0,0);
cnap = CNAaddReactionMFN(cnap,'EX_23bdo_e','1 23bdo_c =' ,0,1000,0,nan,0,...
'//START_GENERIC_DATA{;:;deltaGR_0;#;num;#;NaN;:;uncertGR_0;#;num;#;NaN;:;geneProductAssociation;#;str;#;;:;}//END_GENERIC_DATA',0,0,0,0);

%% 2) Define MCS setup
% reaction indices used in Target and Desired regions
r23BDO_ex = find(strcmp(cellstr(cnap.reacID),'EX_23bdo_e'));
rGlc_ex   = find(strcmp(cellstr(cnap.reacID),'EX_glc__D_e'));
rATPM     = find(strcmp(cellstr(cnap.reacID),'ATPM'));
rBM       = find(cnap.objFunc);
% Target region: First compute maximum possible yield, then setup target 
% region accordingly to demand a minimum yield of 30% of the theoretical maximum.
Ymax_23bdo_per_glc = CNAoptimizeYield(cnap,full(sparse(1,r23BDO_ex,1,1,cnap.numr)),full(sparse(1,rGlc_ex,-1,1,cnap.numr)));
Y_thresh = Ymax_23bdo_per_glc * 0.3;
disp(['Minimum product yield threshold set to ' num2str(Y_thresh)]);
T = full(sparse(  [1         1        ], ... % case: single substrate - glucose
                  [r23BDO_ex rGlc_ex  ], ...
                  [1         Y_thresh ],1,cnap.numr));
t = 0;

% Desired regions: 
% D1: Biomass production rate > 0.05 h^-1, 
% D2: ATPM >= 18 mM/gBDW/h
D1 = full(sparse( 1, rBM,   -1,1,cnap.numr));
d1 = -0.05;
D2 = full(sparse( 1, rATPM, -1,1,cnap.numr));
d2 = -18;

T = {T};
t = {t};
D = {D1 D2};
d = {d1 d2};

% knockables: All reactions with gene rules + O2 exchange as a potential knockout
rkoCost = double(cellfun(@(x) ~isempty(x),CNAgetGenericReactionData_as_array(cnap,'geneProductAssociation')));
rkoCost(rkoCost==0) = nan;
rkoCost(strcmp(cellstr(cnap.reacID),'EX_o2_e')) = 1;
% pseudo-gene that marks spontanous reactions is not knockable
[~,~,genes,gpr_rules] = CNAgenerateGPRrules(cnap);
gkoCost = ones(length(genes),1);
gkoCost(ismember(genes,'spontanous')) = nan;

%% 3) MCS Computatior
tic;
[rmcs, gmcs, gcnap, cmp_gmcs, cmp_gcnap, mcs_idx] = CNAgeneMCSEnumerator2(cnap, T, t, D, d,...
                                                    rkoCost,[], ... reaction KO cost, reaction addition cost
                                                    max_solutions,max_num_interv, ...
                                                    gkoCost,[], ... gene KO cost, gene addition cost
                                                    [],options,... gpr_rules, options
                                                    1); % verbose
comp_time = toc;
disp(['Computation time: ' num2str(comp_time) ' s']);

%% 4) validate MCS
if full(~all(all(isnan(gmcs)))) % if mcs have been found
    disp('Verifying mcs');
	valid = verify_mcs(gcnap,gmcs,gcnap.mcs.T,gcnap.mcs.t,gcnap.mcs.D,gcnap.mcs.d);
else
    valid = [];
end
if ~isempty(getenv('SLURM_ARRAY_TASK_ID'))
    filename = ['desired2-' model '-' getenv('SLURM_JOB_ID')];
else
    filename = ['desired2-' model '-' datestr(date,'yyyy-mm-dd')];
end
save([filename '.mat'],'gcnap','gmcs','valid','comp_time');


%% 5) Characterization and ranking of MCS
% Instead of the gene-MCS, their corresponding reaction-representations are analyzed.
% This is preferred, because the reaction-model is smaller and therefore analysis is 
% faster than in the GPR-extended model. Furthermore different gene-MCS can lead to 
% identical 'phenotypes' when translated to the reaction-model and by analyzing rMCS
% only a reduced, non-redundant set of reaction-MCS needs therefore to be considered.
if full(~all(all(isnan(gmcs)))) % if mcs have been found
    disp('Characterizing mcs');
    cnap.local.c_makro = [];
    cnap.local.c_makro = [];
  % 5.1) Lump redundant MCS and create flux bounds for each mutant model
    rmcs(isnan(rmcs)) = -inf; % this step allows to apply 'unique' too remove duplicates
    [rmcs,~,gmcs_rmcs_map] = unique(rmcs','rows');
    rmcs = rmcs';
    rmcs(rmcs == -inf) = nan;
    MCS_mut_lb = repmat({cnap.reacMin},1,size(rmcs,2));
    MCS_mut_ub = repmat({cnap.reacMax},1,size(rmcs,2));
    MCS_mut_lb = arrayfun(@(x) MCS_mut_lb{x}.*(rmcs(:,x)==1 | rmcs(:,x)==0),1:size(rmcs,2),'UniformOutput',0);
    MCS_mut_ub = arrayfun(@(x) MCS_mut_ub{x}.*(rmcs(:,x)==1 | rmcs(:,x)==0),1:size(rmcs,2),'UniformOutput',0);
  % 5.2) Set relevant indices [criterion 2-7] and prepare thermodynamic (MDF) parameters [criterion 9]
    % reaction indices
    [idx,mdfParam] = relev_indc_and_mdf_Param(cnap);
    idx.prod = r23BDO_ex;
    idx.prodYieldFactor = 1;
    idx.subs = rGlc_ex;
    idx.subsYieldFactor = -1;
    idx.bm      = rBM;
  % 5.3) Define core metabolism [criterion 8]
    % Add the new reactions also to the list of reactions that will be
    % considered "core" reactions in the final MCS characterization and ranking
    new_reacs = ismember(cellstr(cnap.reacID),{'ACLDC','BTDD','EX_23bdo_e'});
    reac_in_core_metabolism(new_reacs) = 1;
    lbCore = cnap.reacMin;
    ubCore = cnap.reacMax;
    lbCore(~reac_in_core_metabolism) = 0;
    ubCore(~reac_in_core_metabolism) = 0;
  % 5.4) Costs for genetic interventions  [criterion 10]
    intvCost                  = gcnap.mcs.kiCost;
    intvCost(isnan(intvCost)) = gcnap.mcs.koCost(isnan(intvCost));
    intvCost(gcnap.rType == 'g') = 1;
    gene_and_reac_names = cellstr(gcnap.reacID);
    gene_and_reac_names(gcnap.rType == 'g') = genes; % to avoid the 'GR-' prefix
  % 5.5) Start characterization and ranking
    [MCS_rankingStruct, MCS_rankingTable]...
        = CNAcharacterizeGeneMCS( cnap , MCS_mut_lb, MCS_mut_ub, 1:size(MCS_mut_lb,2),... model, mutants LB,UB, incices ranked mcs
        idx, idx.cytMet, D, d, T, t, mdfParam, ... relevant indices, Desired and Target regions
        lbCore, ubCore, gmcs, intvCost, gene_and_reac_names, gmcs_rmcs_map, ...
        [0:7 9 10], ones(1,9),2); % assessed criteria and weighting factors
    % save ranking and textual gmcs as tab-separated-values
    cell2csv([filename '.tsv'],MCS_rankingTable,char(9));
    text_gmcs = cell(size(gmcs,2),1);
    for i = 1:size(gmcs,2)
        kos = find(~isnan(gmcs(:,i)) & gmcs(:,i) ~= 0);
        for j = 1:length(kos)
            text_gmcs(i,j) = cellstr(gcnap.reacID(kos(j),:));
        end
    end
    cell2csv([filename '-gmcs.tsv'],text_gmcs,char(9));
    save([filename '.mat'],'gmcs_rmcs_map','MCS_rankingStruct','MCS_rankingTable','cnap','rmcs','T','t','D','d','rkoCost','-append');
end

% clear irrelevant variables
a=[setdiff(who,{'cnap','rmcs','D','d','T','t','compression','filename','gcnap',...
                'gmcs','gmcs_rmcs_map','gpr_rules','rmcs','valid','comp_time'});{'a'}];
rmpath(function_path);
clear(a{:});