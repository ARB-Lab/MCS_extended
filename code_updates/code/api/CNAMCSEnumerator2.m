function [mcs, status, cmp_mcs, cmp_cnap, mcs_idx_cmp_full, obj] = ...
    CNAMCSEnumerator2(cnap,T,t,D,d,koCost,kiCost,maxSolutions,maxCost,options,verbose)
%
% ------------------------------------------------
% CellNetAnalyzer API function 'CNAMCSEnumerator2'
% ------------------------------------------------
% --> Minimal cut set enumeration or (iterative) search with multiple target and 
%     desired regions, reaction knockouts, additions and individual intervention costs
%
% Usage: [mcs, status, cmp_cnap, cmp_mcs, mcs_idx_cmp_full, obj] = ...
%           CNAMCSEnumerator2(cnap,T,t,D,d,koCost,kiCost,maxSolutions,maxCost,options,verbose)
% 
% Computes minimal cut set (MCS) based metabolic engineering strain designs.
% The MCS approach can be used to find interventions that turn wild type
% microbes into production hosts. In contrast to bi-level optimization
% techniques, the design specifications are not represented by an optimization 
% goal but denote the desired (e.g. growth) and undesired flux-states (e.g.
% inferior product yield). The MCS algorithm will find interventions that
% shape the solution space of steady state flux-states according to those.
% The design specifications are presented in the way of one or more sets of 
% target constraints (defined by one or multiple matrices T and vectors t) 
% that span the subspace(s) of undesired steady flux states which should be 
% rendered inaccessible by applying a minimal cut set and optionally by a 
% number of sets of desired constraints (defined by one or multiple matrices 
% D and vectors d) that describe the behavior that should be maintained. MCS 
% can consist of gene and/or reaction additions or deletions, each intervention 
% attributed with an individual cost factor.
%
% Example, how desired behavior can be expressed as a region:
%                             D          * r <=  d
%                   ( 0 -1 0 0 0 0 0 0 ) * r <= -0.1
%            <=>                        -r_2 <= -0.1
%            <=>                         r_2 >=  0.1
%    (- More than one row is allowed for the definition of each Desired/Target region)
%    (- Inequalities can be used to express yield thresholds/constraints for Target/Desired region)
%
% The MCS ensure that all target flux vectors r obeying
%
%       cnap.stoichMat * r = 0
%       cnap.reacMin(r) <= r <= cnap.reacMax(r)
%       T * r <= t
%
% will be removed from the solution space and thus the undesired behaviors
% are blocked. 
%       Target{1}: r_prod <= 1 will, for example, lead to an elimination of all
%                              flux states with production rates inferior to 1.
% 
% Whereas (one or more) flux vectors r fulfilling
%
%       cnap.stoichMat * r = 0
%       cnap.reacMin(r) <= r <= cnap.reacMax(r)
%       D * r <= d
%
% must be kept in the solution space, maintaining desired
% functions/behaviors.
%       Desired{1}: r_growth >= 0.1 will, for example, maintain the capability to
%                                   grow with rates higher than 0.1.
%
% ------------------------------------------------
% Input:
% cnap        : a CellNetAnalyzer mass flow project
% T           : <cell>[numT x 1] Matrices specifying (with vector t) the target flux
%                regions as given above. T = {T1,T2, ...}. T1: <double>[rowsT1 x cnap.numr]
% t           : <cell>[numT x 1] Vectors specifying (with matrix T) the target flux
%                regions as given above. t = {t1,t2, ...}. t1: <double>[rowsT1 x 1]
% D           : <cell>[numD x 1] Matrices specifying (with vector d) the desired flux
%                regions as given above. D = {D1,D2, ...}. D1 <double>[rowsD1 x cnap.numr]
% d           : <cell>[numD x 1] Vectors specifying (with matrix D) the desired flux
%                regions as given above. d = {d1,d2, ...}. d1 <double>[rowsD1 x 1]
% koCost      : <double>[cnap.numr x 1] that indicates reactions that 
%                can be knocked out. Notknockable reactions carry NaN. Knockable
%                reactions carry a value that represents their knockout cost. By
%                default all reactions are knockable.
% kiCost      : <double>[cnap.numr x 1] that indicates reactions that 
%                can be added or "knocked in". Notknockable reactions carry NaN. 
%                Addable reactions carry a value that represents their addition cost.
%                If a reaction is set knockable and addable at the same time. It
%                will be assumed that the reaction is addable. By default no
%                reaction is addible.
% maxSolutions: <double> max. number of solutions. Aborts enumeration after 
%                this number of (compressed) solutions have been found. Inf 
%                is allowed (default).
% maxCost     : <double> Upper cost limit for the sum of KO/KI costs in one
%                MCS. If all interventions have a cost of 1 this value 
%                represents the total number of interventions. Inf is
%                allowed (default) but not advised in networks with 100+ 
%                reactions. For full enumerations in genome scale networks 
%                this is typically between 5 and 8. In an iterative search
%                this number is not limited but would still be set to a
%                value that seems feasible in practice (e.g. is the 
%                implementation of 15 gene interventions manageable?).
% options     : <struct> controls different pre- and post-processing steps and. 
%                MILP parameters. Set the according fields if you want to change 
%                the default settings.
%       options.mcs_search_mode              : <double> 
%                                                1: will iteratively search for single solutions. 
%                                                   This is preferred when a full enumeration of intervention
%                                                   strategies is not necessary or too time consuming. Consider
%                                                   this setting also if you chose a high value for maxSumItvCost.
%                                                2: (default)  will do a full bottom-up stepwise enumeration of MCS.
%       options.compression_network_pre_milp : <logical> Network compression (after GPR rule extension and before running the MILP).
%                                                        (default TRUE)
%       options.preproc_check_feas           : <logical> Check if all target and desired regions are feasible in the original model.
%                                                        (default TRUE) Deactivate if you are sure D and T are feasible and you want
%                                                        to reduce the computation time overhead.
%       options.preproc_D_leth               : <logical> Preprocessing step, where minimal knockout combinations are computed
%                                                        that disrupt one or more desired regions (if provided). These MCS/solutions
%                                                        are excluded in the final MCS-MILP. Ideally this speeds up the search for
%                                                        MCS. Depending on the case, it might yet slow down MCS computation.
%                                                        (default FALSE)
%       options.postproc_verify_mcs          : <logical> Verify MCS after the computation to ensure that no false positives are returned.
%                                                        (default TRUE)
%       options.milp_solver                  : <char> 'matlab_cplex' (default), 'native' (supported in the future)
%       options.milp_time_limit              : <double> in seconds. Inf is allowed (default).
%       options.milp_bigM                    : <logical or double> 
%                                                 DOUBLE: The provided value is used as big M.
%                                                 TRUE: M will be set to 110% of the highest reaction boundary. 
%                                                 FALSE or zero: indicator constraints are used (default).
%     == These MILP settings can only be used with the MATLAB-CPLEX API (see options.milp_solver) ==
%       options.milp_split_level             : <logical> Parameter for MILP construction (default FALSE) supported only for matlab_cplex
%       options.milp_reduce_constraints      : <logical> Parameter for MILP construction (default TRUE) supported only for matlab_cplex
%       options.milp_combined_z              : <logical> Parameter for MILP construction (default TRUE) supported only for matlab_cplex
%       options.milp_irrev_geq               : <logical> Parameter for MILP construction (default TRUE) supported only for matlab_cplex
% verbose      : <double> 1 yes (default), 0 no. Select 0 for a silent run.
%
%
% ------------------------------------------------
% Output:
% mcs          : <double>[cnap.numr x numMCS] 0: reaction was left unchanged; 1:
%                 reactions was added; -1: reaction was removed; NaN: reaction
%                 was addition candidate but not added.
% status       : <double> 0 successful, 1 timeout without any solution, 
%                 2 infeasible setup, 3 timeout but (some) solutions were found.
% cmp_mcs      : MCS that were found on the compressed setup.
% cmp_cnap     : compressed CNA mass-flow project used for MCS computation.
%                 The MCS setup is provided in the fields cmp_cnap.mcs.*
% mcs_idx_cmp_full: <double>[1 x numMCS] A vector that maps the decompressed 
%                   mcs to the original compressed ones.
% obj          : Object that contains the MILP setup.
%

%
% This file is part of CellNetAnalyzer. Please visit
% http://www.mpi-magdeburg.mpg.de/projects/cna/cna.html
% for more information and the latest version of CellNetAnalyzer.
%
% Copyright (C) 2000-2020 by Steffen Klamt and Axel von Kamp,
% Max Planck Institute for Dynamics of Complex Technical Systems, Magdeburg, Germany.
%
% Contributors are listed in CONTRIBUTORS.txt.
%
% This software can be used under the terms of our CellNetAnalyzer License.
% A copy of the license agreement is provided in the file named "LICENSE.txt"
% included with this software distribution. The license is also available online at
% http://www2.mpi-magdeburg.mpg.de/projects/cna/license.html
%
% For questions please contact: cellnetanalyzer@mpi-magdeburg.mpg.de
%

if nargin < 11 || isempty(verbose)
    verbose = 1;
end
if isa(T,'double') && isa(t,'double')  % compatibility to old MCSEnumerator
    T = {T};
    t = {t};
end
if nargin < 5 || isempty(D)
    lb_D={};
    ub_D={};
    D = {};
    d = {};
else
    if isa(D,'double') && isa(d,'double') % compatibility to old MCSEnumerator
        D = {D};
        d = {d};
    end
end
if nargin < 6 || isempty(koCost)
    koCost = ones(1,cnap.numr);
elseif length(koCost) ~= cnap.numr % compatibility to previous MCSEnumerator
    notknock_idx = koCost;         % -> notknockables can be defined as reaction indices
    koCost = ones(1,cnap.numr);
    koCost(notknock_idx) = nan;
end
if nargin < 7 || isempty(kiCost)
    kiCost = nan(1,cnap.numr);
end
% Transform ki & ko vector to row vector
kiCost = kiCost(:)';
koCost = koCost(:)';
koCost(~isnan(kiCost)) = nan; % knock-ins 'override' knock-outs

if nargin < 8 || isempty(maxSolutions)
    maxSolutions = inf;
end
if nargin < 9 || isempty(maxCost)
    maxCost = inf;
end
if nargin < 10 || isempty(options)
    options = struct;
end
if isa(options,'struct')
    if ~isfield(options,'mcs_search_mode')
        options.mcs_search_mode = 2;
    end
    if ~isfield(options,'compression_network_pre_milp')
        options.compression_network_pre_milp = true;
    end
    if ~isfield(options,'preproc_check_feas')
        options.preproc_check_feas = true;
    end
    if ~isfield(options,'preproc_D_leth') || isempty(D)
        options.preproc_D_leth = false;
    end
    if ~isfield(options,'postproc_verify_mcs')
        options.postproc_verify_mcs = true;
    end
    if ~isfield(options,'milp_time_limit')
        options.milp_time_limit = inf;
    end
    if ~isfield(options,'milp_bigM') || options.milp_bigM == 0
        options.milp_bigM = false;
    end
    if ~isfield(options,'milp_solver')
        options.milp_solver = 'matlab_cplex';
    end
    if strcmp(options.milp_solver,'matlab_cplex') 
        if ~isfield(options,'milp_split_level')
            options.milp_split_level = false;
        end
        if ~isfield(options,'milp_reduce_constraints')
            options.milp_reduce_constraints = true;
        end
        if ~isfield(options,'milp_combined_z')
            options.milp_combined_z = true;
        end
        if ~isfield(options,'milp_irrev_geq')
            options.milp_irrev_geq = true;
        end
    elseif strcmp(options.milp_solver,'java_cplex')
        if any( [isfield(options,'milp_split_level') isfield(options,'milp_reduce_constraints') ...
                 isfield(options,'milp_combined_z')  isfield(options,'milp_irrev_geq') ] )
            displ('Some MILP settings can only be uses with the MATLAB-Cplex-API and will be ignored here.',verbose)
            try options=rmfield(options,'milp_split_level'); catch, end
            try options=rmfield(options,'milp_reduce_constraints'); catch, end
            try options=rmfield(options,'milp_combined_z'); catch,  end
            try options=rmfield(options,'milp_irrev_geq'); catch,   end
        end
    elseif strcmp(options.milp_solver,'native')
        error('Native MATLAB MILP solver ''linintprog'' is not yet supported');
    else
        error(['Solver ' options.milp_solver ' is not supported']);
    end

else
    error('Function parameter ''options'' must be either a struct or empty.');
end

displ('== MCS Computation ==',verbose);
displ('Settings:',verbose);
displ(options,verbose);

%include stoichiometry for growth for model with explicit biomass composition
if(~isempty(cnap.mue))
	cnap.stoichMat(:,cnap.mue)=cnap.macroComposition*cnap.macroDefault;
end

%% 0. remove external specs and test feasibility of target and desired vectors
if any(cnap.specExternal)
    displ('Removing external species.',verbose);
    cnap = CNAdeleteSpecies(cnap,find(cnap.specExternal),0);
end
if options.preproc_check_feas
    displ('Verifying that D and T regions are feasible in the original model.',verbose);
    feas_D = testRegionFeas(cnap,D,d,2);
    feas_T = testRegionFeas(cnap,T,t,2);
    if any(~feas_T)
        error(['At least one target region (T' num2str(find(~feas_T)) ') is infeasible in the original model.']);
    end
    if any(~feas_D)
        error(['At least one desired region (D' num2str(find(~feas_D)) ') is infeasible in the original model.']);
    end
end

%% 1. Determine boundaries to whole model. (If model is not compressed, this FVA also helps bounding the target regions (5))
displ('FVA to find blocked reactions in uncompressed model.',verbose);
r_notBlocked = find(~(cnap.reacMin == 0 & cnap.reacMax ==0));
reacMin_FVA = zeros(cnap.numr,1);
reacMax_FVA = zeros(cnap.numr,1);
[reacMin_FVA(r_notBlocked), reacMax_FVA(r_notBlocked)] = CNAfluxVariability(cnap,[],cnap.macroDefault,2,r_notBlocked);
displ(['#total_reacs: ' num2str(cnap.numr) ...
      ', #blocked: ' num2str(sum(reacMin_FVA == 0 & reacMax_FVA == 0) + sum(cnap.reacMin == 0 & cnap.reacMax ==0))],verbose);
% All reaction bounds could actually be set here. But to avoid numerical issues
% during compression, only exact zeros are used to identify blocked reactions
% to remove them during compression.
reac_off = find(reacMin_FVA == 0 & reacMax_FVA == 0);
cnap.reacMin(reac_off) = 0;
cnap.reacMax(reac_off) = 0;
%% 2. Flip irreversible negative reactions (lb <= r <= 0) to positive direction (0 <= -r <= -lb)
% saves integer variables and reaction boundary constraints in the Target region later on.
r_irrev_neg = find(reacMax_FVA <= 0 & reacMin_FVA < 0);
[cnap, D, T, flip_map] = flip_negative_reacs(cnap,D,T,r_irrev_neg);
displ(['Flipped ' num2str(length(r_irrev_neg)) ' irreversible reactions.'],verbose);
rmax = reacMax_FVA(r_irrev_neg);
reacMax_FVA(r_irrev_neg) = -reacMin_FVA(r_irrev_neg);
reacMin_FVA(r_irrev_neg) = -rmax;
% flip_map = eye(cnap.numr);

%% 4. compress network
if options.compression_network_pre_milp
    kiCost_uncmp = kiCost;
    koCost_uncmp = koCost;
    % replace the full model variables by the compressed ones
    displ('Compressing network.',verbose);
    [cnap, T, D, koCost, kiCost,  compression_map] = ...
        compress(cnap,T,D,koCost,kiCost,reac_off);
    compression_map = flip_map*compression_map;

    % If desired FVAs were done before compression, boundaries could be
    % remapped with this function:
    %     for i = 1:length(D)
    %         [lb_D{i}, ub_D{i}] = remapBounds(lb_D{i},ub_D{i},compression_map);
    %     end
    % test again feasibility of target and desired vectors
    if options.preproc_check_feas
        displ('Verifying D and T region feasibility in compressed model.',verbose);
        feas_D = testRegionFeas(cnap,D,d,2);
        feas_T = testRegionFeas(cnap,T,t,2);
        if any(~feas_T)
            error(['At least one target region (T' num2str(find(~feas_T)) ') is infeasible in the original model.']);
        end
        if any(~feas_D)
            error(['At least one desired region (D' num2str(find(~feas_D)) ') is infeasible in the original model.']);
        end
    end
    cmp_cnap = cnap;
    cmp_cnap.mcs.T = T;
    cmp_cnap.mcs.t = t;
    cmp_cnap.mcs.D = D;
    cmp_cnap.mcs.d = d;
    cmp_cnap.mcs.kiCost = kiCost;
    cmp_cnap.mcs.koCost = koCost;
    %% (3.1) repeat FVA for model bounds (to be integrated in Target region)
    displ('FVA to determine model bounds.',verbose);
    [reacMin_FVA, reacMax_FVA] = CNAfluxVariability(cnap,[],cnap.macroDefault,2);
else
    compression_map = flip_map;
    cmp_cnap = cnap;
    cmp_cnap.mcs.T = T;
    cmp_cnap.mcs.t = t;
    cmp_cnap.mcs.D = D;
    cmp_cnap.mcs.d = d;
    cmp_cnap.mcs.kiCost = kiCost;
    cmp_cnap.mcs.koCost = koCost;
end

%% Determine boundaries under desired condition and identify further notknockables
displ('FVA to determine model bounds under desired conditions.',verbose && ~isempty(D));
essential = false(cnap.numr,1);
for i = 1:length(D)
    [lb_D{i}, ub_D{i}] = CNAfluxVariability(cnap,[],cnap.macroDefault,2,1:cnap.numr,D{i},d{i});
    % A reaction is essential when its lower bound sign equal.
    % A threshold is used to avoid marking reactions essential that perhaps
    % aren't. (e.g. lb = +1e-14) This precautious measure might lead to
    % more knockables but cannot lead to any false MCS.
    essential = essential | ...
                (sign((abs(lb_D{i})>cnap.epsilon).*lb_D{i}) .* ...
                 sign((abs(ub_D{i})>cnap.epsilon).*ub_D{i}) == 1);
end
displ(['#essential: ' num2str(sum(essential)) ' -> ' num2str(sum(~isnan(koCost(essential)))) ' more reactions are now non-targetable.'],verbose && ~isempty(D));
koCost(essential) = nan; % make essential reactions "notknockable"
knockable = ~isnan(koCost) | ~isnan(kiCost);
% remove redundant desired bounds
for i = 1:length(D)
    % Remove bounds of reactions when simultaneously (1) FVA bound is not identical with the model bound, 
    % because these bounds are implicit and can be derived from the truely limiting bounds. (2) reactions that can
    % be knocked out/in, because they need bounds in the MILP to simulate KOs.
    explicit_D_lb = (abs(lb_D{i}-cnap.reacMin) <= cnap.epsilon)';
    explicit_D_ub = (abs(ub_D{i}-cnap.reacMax) <= cnap.epsilon)';
    % Desired FVA bounds can be replaced with the original model bounds to 
    % avoid numerical issues from the FVA. It works but slows down computation ~10 fold
    lb_D{i}(~explicit_D_lb & ~knockable) = -inf; % unbound all other reactions
    ub_D{i}(~explicit_D_ub & ~knockable) =  inf; % unbound all other reactions
%     ub_D{i}(explicit_D_ub | knockable) = cnap.reacMax(explicit_D_ub | knockable); % replace FBA bounds with model bounds
%     lb_D{i}(explicit_D_lb | knockable) = cnap.reacMin(explicit_D_lb | knockable); % replace FBA bounds with model bounds
    % MILP requires non-inf bounds for knockable reactions. If thouse bounds 
    % are larger than 10^5, numerical issues might occur, so tighten these.
    if any(knockable & (isinf(lb_D{i})' | isinf(ub_D{i})'))
        warning(['The constraint-MCS-MILP requires all knockable reactions to be bound. ' ...
                 'FVA is used to assign upper bounds to these reactions. However, FVA ' ...
                 'found that the flux is not limited in one or more knockable reactions. ' ...
                 'As a workaround, their flux range is limited to +1000 and/or -1000. If you face ' ...
                 'numerical issues or MCS computation does not finish successfully, make sure '...
                 'to bounds to all cycles and pathways that contain knockable reactions '...
                 'before you enter MCS computation.']);
        lb_D{i}(knockable & isinf(lb_D{i})') = -1e3; % Necessary to simulate knockouts.
        ub_D{i}(knockable & isinf(ub_D{i})') =  1e3; % bound knockable reactions even if FVA shows they're unbound. Necessary to simulate knockouts.
    end
end
if ~isempty(D) && any(any(isinf([cell2mat(lb_D) cell2mat(ub_D)]),2) & (~isnan(koCost) | ~isnan(kiCost))')
    error(['Knockable reactions must be bound in the desired system. But the follwing reactions are unbound: '...
        strjoin(cellstr(cnap.reacID(any(isinf([cell2mat(lb_D) cell2mat(ub_D)]) & (~isnan(koCost) | ~isnan(kiCost)),2),:)),', ')])
end

%% 5. integrate (explicit) flux bounds in Target system
% Flux boundaries are relevant for the Target constraints if they are truly flux limiting, 
% then those limits can be reached in an FVA, otherwise these bounds are redundant to the other
% bounds and the network structure.
lb_relevant = abs(reacMin_FVA-cnap.reacMin) <= cnap.epsilon & cnap.reacMin ~= 0;
ub_relevant = abs(reacMax_FVA-cnap.reacMax) <= cnap.epsilon;
r_irrev_pos = reacMin_FVA' >= 0;
% workaround
% r_irrev_pos = zeros(1,cnap.numr);
if any(reacMin_FVA < 0 & reacMin_FVA > -cnap.epsilon)
    warning('Target region: Some reactions might be irreversible but not treated as such. This has probably numerical reasons.');
end
% Bound target regions with model bounds:
% -r <= -lb
%  r <=  ub
T_entries = [-lb_relevant , ub_relevant]';
T_bounds  = [-cnap.reacMin' ; cnap.reacMax']; % [-reacMin_FVA'; reacMax_FVA']; % [-cnap.reacMin' ; cnap.reacMax']; % Using FVA bounds is identical: [-reacMin_FVA'; reacMax_FVA']; . Consider when numerical issues occur.
[~,T_reac] = find(T_entries);
T_i = sparse(1:length(T_reac),T_reac,T_entries(T_entries~=0),length(T_reac),cnap.numr);
t_i = T_bounds(T_entries~=0);
for i = 1:length(T)
    T{i}(size(T{i},1)+(1:length(T_reac)),:) = T_i; % append target matrix
    t{i}(size(t{i},1)+(1:length(T_reac)),:) = t_i; % append target vector
end

%% Initialize CPLEX parameters (memory and CPU)
% Memory allocation
if isunix
    slurm_job_mem = str2double(getenv('SLURM_REQMEM'));
    if ~isnan(slurm_job_mem) % if on SLURM node with allocated memory
        use_mem = round(slurm_job_mem*0.80); % use 80% of allocated memory
%         delete(gcp('nocreate')); % stop parpool to free additional memory
    elseif ~isempty(getenv('SLURM_JOBID')) % if on SLURM node (normal)
        [~,w] = unix('free | grep Mem');
        stats = str2double(regexp(w, '[0-9]*', 'match'));
        use_mem = round(stats(1)/1e3*0.8); % use 80% of total memory
    else
        [~,w] = unix('free | grep Mem');
        stats = str2double(regexp(w, '[0-9]*', 'match'));
        use_mem = round(stats(6)/1e3*0.6); % use 60% of available memory
    end
else
    [~,use_mem] = memory(); % use 90% of the available memory
    use_mem = round(use_mem.PhysicalMemory.Available*1e-6 * 0.9);
end
slurm_job_cpus = str2double(getenv('SLURM_CPUS_ON_NODE'));
if ~isnan(slurm_job_cpus)
    use_threads = slurm_job_cpus;
else
    use_threads = (2*feature('numcores'))-1; % Use one thread less than virtual cores available.
end

%% Exclude desired disrupting knockout-combinations
% Add this by changing the flag 'comp_D_lethalities' if you intend to exclude lethal combinations:
if options.preproc_D_leth
    displ('Searching for knockout combinations that would disrupt desired behavior.',verbose);
    cplex_param.usemem = use_mem;
    cplex_param.usecpus = use_threads;
    desired_disrupt_KO_combin = findDesiredDisrupting(cnap,D,d,lb_relevant,ub_relevant,~isnan(koCost),2,cplex_param,options); % last parameter: Max size of desired-disrupting MCS
end
    
%% 5. Construct MILP
displ('Building MCS MILP.',verbose);
% bigM
if islogical(options.milp_bigM) && (options.milp_bigM == true)
    M = max(max(abs([cnap.reacMin cnap.reacMax]))) * 1.1; % Choose big M 10% bigger than largest bound
elseif options.milp_bigM ~= 0
    M = options.milp_bigM;
else
    M = 0;
end
switch options.milp_solver
    case 'matlab_cplex'
        obj= ConstrainedMinimalCutSetsEnumeratorMILP(  cnap.stoichMat, ...
            r_irrev_pos, [], ...  irreversibles
            T, t, ...             target
            {~isnan(koCost), full(sparse(1,find(~isnan(koCost)),koCost(~isnan(koCost)),1,length(koCost)))}, ... KO
            lb_D, ...        desired bounds
            ub_D, ...
            D, d, ...        desired
            {~isnan(kiCost), full(sparse(1,find(~isnan(kiCost)),kiCost(~isnan(kiCost)),1,length(kiCost)))},... KI
            M,1, ... bigM threshold
            options.milp_split_level, ...
            options.milp_reduce_constraints, ...
            options.milp_combined_z, ...
            options.milp_irrev_geq);
        if M
            obj.cpx.Param.mip.tolerances.integrality.Cur = 1e-8;
        end
        % cplex_param = cplexoptimset('cplex'); % to view all parameters
        obj.cpx.Param.threads.Cur = use_threads;
        obj.cpx.Param.workmem.Cur = use_mem;
        obj.cpx.DisplayFunc = [];
        obj.cpx.Param.output.clonelog.Cur = -1;
        % CPLEX parameters. You can try to increase performance by tuning these.
        % obj.cpx.Param.cplex_param.barrier.algorithm.Cur = 3;
        % obj.cpx.Param.mip.strategy.startalgorithm.Cur = 2; % Activate if CPLEX crashes
        % obj.cpx.Param.emphasis.memory.Cur = 1; % Activate if Out of Memory is an issue
        % obj.cpx.Param.mip.strategy.nodeselect.Cur = 0; % Activate if Out of Memory is still an issue.
        %                                                          Using "Dept first" (0) node selection leads to much longer runtimes, 
        %                                                          but consumes much less memory.
        % obj.cpx.Param.mip.strategy.lbheur.Cur = 1; % Heuristics
        % obj.cpx.Param.mip.display.Cur = 2; % Ignored when setOut = []
        % obj.cpx.Param.emphasis.mip.Cur = 1; % Feasibility emphasis
        % obj.cpx.Param.mip.strategy.probe.Cur = 3; % Agressive probing speeds up binary problems
        numMILPvars = size(obj.cpx.Model.A,2);
        knockable = obj.cpx.Model.ub(obj.z_vars);
    case 'java_cplex'
        obj= ConstrainedMinimalCutSetsEnumerator(  cnap.stoichMat, ...
            r_irrev_pos, [], ...  irreversibles
            T, t, ...             target
            {~isnan(koCost), full(sparse(1,find(~isnan(koCost)),koCost(~isnan(koCost)),1,length(koCost)))}, ... KO
            lb_D, ...        desired bounds
            ub_D, ...
            D, d, ...        desired
            {~isnan(kiCost), full(sparse(1,find(~isnan(kiCost)),kiCost(~isnan(kiCost)),1,length(kiCost)))},... KI
            M);
        cplex_inner = setup_cplex_inner_class_access();
        cplex_param = ilog.cplex.IloCplex().getParameterSet;
        if M
            cplex_param.setParam(cplex_inner.DoubleParam.EpInt, 1e-8);
        end
        % CPLEX parameters. You can try to increase performance by tuning these.
        cplex_param.setParam(cplex_inner.DoubleParam.WorkMem, use_mem); % Allocate memory
        cplex_param.setParam(cplex_inner.IntParam.Threads, use_threads); % Use all allocated SLURM cores.
        % cplex_param.setParam(cplex_inner.IntParam.BarAlg,3);
        % cplex_param.setParam(cplex_inner.IntParam.RootAlg,2); % Activate if CPLEX crashes
        % cplex_param.setParam(cplex_inner.BooleanParam.MemoryEmphasis, 1); % Activate if Out of Memory is an issue
        % cplex_param.setParam(cplex_inner.IntParam.NodeSel, 0); % Activate if Out of Memory is still an issue.
        %                                                          Using "Dept first" (0) node selection leads to much longer runtimes, 
        %                                                          but consumes much less memory.
        % cplex_param.setParam(cplex_inner.IntParam.FPHeur,1); % Heuristics
        % cplex_param.setParam(cplex_inner.IntParam.MIPDisplay, 2); % Ignored when setOut = []
        % cplex_param.setParam(cplex_inner.IntParam.MIPEmphasis, cplex_inner.MIPEmphasis.Feasibility);
        % cplex_param.setParam(cplex_inner.IntParam.Probe, 3); % Agressive probing speeds up binary problems
        obj.cpx.setParameterSet(cplex_param);
        obj.cpx.setOut([]);
        numMILPvars = obj.cpx.getNcols;
        knockable = arrayfun(@(x) x.getUB(),obj.z_vars);
end
if all(isnan(koCost)) && all(isnan(kiCost))
    error('No targetable reactions. Cannot proceed with MCS computation.');
else
    ivCost = kiCost;
    ivCost(~isnan(koCost)) = koCost(~isnan(koCost));
    ivCost(isnan(ivCost))  = 0;
    if isinf(maxCost)
        maxCost = sum(ivCost);
    end
end

% exclude desired desrupting KO-combinations
if options.preproc_D_leth
    obj.add_evs_exclusion_constraints(desired_disrupt_KO_combin, sum(desired_disrupt_KO_combin,1));
end

% Display details on MILP setup size
displ('MILP:',verbose)
milp_info.stoichiometric_matrix_rows = size(cnap.stoichMat,1);
milp_info.stoichiometric_matrix_columns = size(cnap.stoichMat,2);
milp_info.milp_vars = numMILPvars;
milp_info.num_targetable = sum(knockable~=0);
milp_info.time_limit = options.milp_time_limit;
milp_info.milp_memory_mb = use_mem;
milp_info.threads = use_threads;
displ(milp_info,verbose);

%% 5. Compute intervention strategies
starttime = now;
%% 5a iterative search
if options.mcs_search_mode == 1
    disp_text_subspace = @(a,b,c,e) displ([num2str(round(86400*(b - a)*100)/100,'%.2f') ' seconds : ' ...
        'Searching in subspace with ' num2str(sum(c)) ' targetable reactions.'],e);
    disp_text_full_space = @(a,b,c) displ([num2str(round(86400*(b - a)*100)/100,'%.2f') ' seconds : ' ...
        'Searching in full problem space with all targetable reactions.'],c);
    status = 0;
    mcs = nan(cnap.numr,0);
    switch options.milp_solver
        case 'matlab_cplex'
            obj.set_evs_sz_bounds(0,maxCost);
        case 'java_cplex'
            obj.evs_sz.setUB(maxCost);
    end
    endtime = now*86400 + options.milp_time_limit;
    % search for MCS iteratively
    displ('Computing MCS (iteratively).',verbose);
    while size(mcs,2) < maxSolutions && status == 0
        % Search for any (not necessarily minimal) cut set
        disp_text_full_space(starttime,now,verbose);
        obj = set_ub_z_vars(obj,knockable,options);
        [sol, status] = find_mcs(obj, zeros(length(obj.z_vars),1), endtime-now*86400,options);
        if status == 0
            % Minimize found solution by search in subspace
            disp_text_subspace(starttime,now,sol,verbose);
            obj = set_ub_z_vars(obj,sol,options);
            while status == 0 && (size(mcs,2) < maxSolutions)
                [sol, status] = find_mcs(obj, ivCost'.*sol, min(15,endtime-now*86400),options); % spend at most 15 seconds in subspace before refining it
                if  status == 0 && any(sol)
                    % Verify solution
                    [feas_T, feas_D] = verify_mcs(cnap,sol,T,t,D,d,koCost,kiCost,verbose);
                    if all(cell2mat(feas_D)) && all(~cell2mat(feas_T)) && ivCost*sol <= maxCost
                        mcs = [mcs, sol];
                        if verbose % Console output MCS
                            displ(mcs2text(sol,cnap.reacID,koCost,kiCost),verbose);
                            text = [num2str(round(86400*(now - starttime)*100)/100,'%.2f') ' seconds : '];
                            weighted_mcs = round(ivCost*mcs*1000)/1000;
                            for mcs_size = unique(weighted_mcs)
                                text = [text, num2str(sum(weighted_mcs==mcs_size)), ' of cost ', num2str(mcs_size) ,';  '];
                            end
                            displ(text,verbose);
                        end
                    elseif any(~cell2mat(feas_D))
                        displ('Invalid solution found (Desired constraint(s) violated). Solution excluded from pool.',verbose);
                    elseif any(cell2mat(feas_T)) && all(cell2mat(feas_D)) % exclude case (needs to be implemented)
                        warning('Invalid solution found (Target not eliminated). Solution excluded from pool. Correct solutions might be overlooked from here on!');
                    elseif any(cell2mat(feas_T)) && any(~cell2mat(feas_D)) % exclude case (needs to be implemented) 
                        warning('Invalid solution found (Target not eliminated). Solution excluded from pool. Correct solutions might be overlooked from here on!');
                    else % if solution was too large
                        displ('.',verbose);
                    end
                    obj.add_evs_exclusion_constraints(sol, sum(sol));
                elseif status == 3 && any(sol) % Refine subspace again
                    disp_text_subspace(starttime,now,sol,verbose);
                    obj.cpx.Model.ub(obj.z_vars) = sol;
                    status = 0;
                end
            end
            status = 0;
        end
    end
    fprintf(char(10));
%% 5b full stepwise enumeration
elseif options.mcs_search_mode == 2
    [obj, mcs, ~]= obj.shortest_minimal_cut_sets(maxCost, maxSolutions, options.milp_time_limit, 2);
    status = getMILPstatus(obj,options);
end
displ(['Pure MILP computation time: ' num2str(86400*(now - starttime)) ' s'],verbose);
global time_new_milp;
time_new_milp = 86400*(now - starttime);

if size(mcs,2) == maxSolutions
    displ('Solution pool limit reached.',verbose);
end

%% Validate MCS and remove redundant
if options.postproc_verify_mcs
    displ('Verifying compressed MCS.',verbose);
    [feas_T, feas_D] = verify_mcs(cnap,mcs,T,t,D,d,koCost,kiCost,verbose);
    valid = all(~cell2mat(feas_T),2) & all(cell2mat(feas_D),2);
    if any(~valid)
        displ(['Removing ' num2str(sum(~valid)) ' invalid mcs.'],verbose)
        mcs = mcs(:,valid);
    elseif size(mcs,2)>0
        displ('All MCS valid.',verbose);
    end
end
% redundant MCS can only occur if costs for certain interventions are zero
if ~all(all(isnan(mcs))) && any(koCost==0) && any(kiCost==0)
    redund = find_redundant_mcs(mcs);
    if any(redund)  
        mcs = mcs(:,~redund);
        displ(['Removed ' num2str(sum(redund)) ' redundant MCS.'],verbose);
    elseif ~all(all(isnan(mcs)))
        displ('All MCS unique.',verbose);
    end
end


%% 6. decompress intervention strategies
if ~all(all(isnan(mcs)))
    % preparing output of the compressed MCS
    cmp_mcs = mcs;
    % set knock-in/out candidates to [-1 : knocked out] or 
    %                                [1  : knocked in] or 
    %                                [nan: not knocked in]
    cmp_mcs = double(cmp_mcs);
    cmp_mcs(~isnan(koCost),:) = -cmp_mcs(~isnan(koCost),:);
    for i = find(~isnan(kiCost(:)))'
        cmp_mcs(i, cmp_mcs(i,:)==0) = nan;
    end
    
    if options.compression_network_pre_milp
        % decompression
        compSols = size(mcs,2);
        displ('Expand MCS.',verbose);
        kiCost = kiCost_uncmp;
        koCost = koCost_uncmp;
        ivCost = kiCost;
        ivCost(~isnan(koCost)) = koCost(~isnan(koCost));
        [mcs, dum, mcs_idx_cmp_full] = expand_mcs(mcs, compression_map');
        % check if knock-out and knock-in costs are still met
        ivCost(isnan(ivCost))  = maxCost+1; % higher number, so that mcs containing this intervention are filtered
        mcs_affordable = ivCost*mcs <= maxCost;
        mcs = mcs(:,mcs_affordable);
        mcs_idx_cmp_full = mcs_idx_cmp_full(mcs_affordable);
    else
        mcs_idx_cmp_full = 1:size(mcs,2);
    end
    % set knock-in/out candidates to [-1 : knocked out] or 
    %                                [1  : knocked in] or 
    %                                [nan: not knocked in]
    mcs = double(mcs);
    mcs(~isnan(koCost),:) = -mcs(~isnan(koCost),:);
    for i = find(~isnan(kiCost(:)))'
        mcs(i, mcs(i,:)==0) = nan;
    end
    
    %% 6. Prepare Model output and status report
    if size(mcs,2) > 0 && (status == 0 || status == 1 || status == 2)
        if status == 2 % In this case the MILP was not infeasible but returned infeasible because all solutions were already found
            status = 0;
        end
        if options.compression_network_pre_milp
            text =  ['. (' num2str(compSols) ' before network decompression / MCS expansion).']; 
        else
            text = [];
        end
        displ(['MCS found: ' num2str(size(mcs,2)) text],verbose);

    end
    displ('The time limit was reached.', verbose && status == 1);
else
    cmp_mcs = [];
    mcs_idx_cmp_full = [];
    switch status
        case 1
            displ('No MCS found: Timeout.',verbose);
        case 2
            displ('No MCS found: Infeasible.',verbose);
    end
end
end

%% Functions:
%% 1. Flip reactions
function [cnap, D, T, flip_map] = flip_negative_reacs(cnap,D,T,flipped_reacs)
% identify reactions to be flipped and create map
flip_map = eye(cnap.numr);
for i = flipped_reacs'
    flip_map(i,i) = -1;
end
% flipping reactions in Model
rmax = cnap.reacMax(flipped_reacs);
cnap.reacMax(flipped_reacs) = -cnap.reacMin(flipped_reacs);
cnap.reacMin(flipped_reacs) = -rmax;
cnap.stoichMat(:,flipped_reacs) = -cnap.stoichMat(:,flipped_reacs);
for i = flipped_reacs'
    rname = [strtrim(cnap.reacID(i,:)) '_flipped'];
    cnap.reacID(i,1:length(rname)) = rname;
end
% flipping Target and desired
for i = 1:length(T)
    T{i}(:,flipped_reacs) = -T{i}(:,flipped_reacs);
end
for i = 1:length(D)
    D{i}(:,flipped_reacs) = -D{i}(:,flipped_reacs);
end
end

%% 2. compress
function [cmp_cnap, cmp_T, cmp_D, cmp_koCost, cmp_kiCost, cmp_mapReac] = compress(cnap,T,D,koCost,kiCost,r_off)
    % non_compress_reacs = []; % compress all reactions
    % Add this to array if reactions from target and desired regions or
    % knock-in reactions should be exempt from compression: 
    non_compress_reacs = unique(find(any([cell2mat(D') ; cell2mat(T') ; ~isnan(kiCost)],1)))';
    %
    % Add this if you get any Java Errors/Warnings:
    %     javastderr= java.lang.System.err;
    %     java.lang.System.setErr(java.io.PrintStream('cplex_stderr.log'));
    %% Compress
    % remove conservation relations: on
    % remove choke points: off
    % use rational numbers to avoid numeric issues: on
    % suppress console output: on
    [dum1,dum2,cmp_mapReac,dum3,cmp_cnap]=CNAcompressMFNetwork(cnap,non_compress_reacs,[],1,0,1,r_off,0);
    %     java.lang.System.setErr(javastderr);
    %
    %% Remap MCS specifiers
    cmp_T = cellfun(@(x) x*cmp_mapReac,T,'UniformOutput',0);
    cmp_D = cellfun(@(x) x*cmp_mapReac,D,'UniformOutput',0);
    lumpedReacs = double(cmp_mapReac ~= 0);
    lumpedReacs(lumpedReacs == 0) = nan;
    is_kiable = any(~isnan(lumpedReacs.*repmat(kiCost',1,cmp_cnap.numr)),1);
    cmp_kiCost(~is_kiable) = nan; % all reactions where none of the lumped reactions have KI costs are not addable
    cmp_kiCost( is_kiable) = arrayfun(@(x) sum(kiCost(~isnan(kiCost) & ~isnan(lumpedReacs(:,x))')),find(is_kiable)); % for KI, all lumped reactions need to be knocked in
    cmp_koCost = min(lumpedReacs.*repmat(koCost',1,cmp_cnap.numr),[],1); % for KO, only cheapest of the lumped reactions needs to be knocked out
    cmp_koCost(~isnan(cmp_kiCost)) = nan; % KIs overwrite KOs because the lumped reaction can be blocked without cost by not performing a KI
end

%% 3. Find Desired disrupting knockouts
% search for lethal combinations up to size 2
% (don't consider knock-ins / act if all reactions were active)
function leth_combin = findDesiredDisrupting(cnap,D,d,lb_relevant,ub_relevant,knockable,numKOs,cpxparam,options)
% Add this to generate more desired-lethals
    D_entries = [-lb_relevant , ub_relevant]';
    D_bounds  = [-cnap.reacMin' ; cnap.reacMax'];
    [~,D_reac] = find(D_entries);
    D_i = sparse(1:length(D_reac),D_reac,D_entries(D_entries~=0),length(D_reac),cnap.numr);
    d_i = D_bounds(D_entries~=0);
    for i = 1:length(D)
        D{i}(size(D{i},1)+(1:length(D_reac)),:) = D_i; % append target matrix
        d{i}(size(d{i},1)+(1:length(D_reac)),:) = d_i; % append target vector
    end
    leth_combin = nan(size(cnap.stoichMat,2),0);
    for i = 1:length(D)
        switch options.milp_solver
            case 'matlab_cplex'
                obj= ConstrainedMinimalCutSetsEnumeratorMILP(  cnap.stoichMat, ...
                cnap.reacMin'>=0, [], ...  irreversibles
                D{i}, d{i}, ...             target
                knockable, ... KO
                {}, {}, {}, {}, ...        desired
                [],... %  KI
                0,1,options.milp_split_level,options.milp_reduce_constraints,options.milp_combined_z);
                obj.cpx.Param.threads.Cur = cpxparam.usecpus;
                obj.cpx.Param.workmem.Cur = cpxparam.usemem;
                obj.cpx.Param.output.clonelog.Cur = -1;
                obj.cpx.DisplayFunc = [];
            case 'java_cplex'
                obj= ConstrainedMinimalCutSetsEnumerator(  cnap.stoichMat, ...
                    cnap.reacMin'>=0, [], ...  Irr
                    D{i}, d{i}, ...  Target - Here Desired region because we compute lethals
                    knockable);
                cplex_inner = setup_cplex_inner_class_access();
                cplex_param = ilog.cplex.IloCplex().getParameterSet;
                cplex_param.setParam(cplex_inner.IntParam.Threads, cpxparam.usecpus); % Use all allocated SLURM cores.
                cplex_param.setParam(cplex_inner.DoubleParam.WorkMem, cpxparam.usemem); % Allocate memory
                obj.cpx.setParameterSet(cplex_param);
                obj.cpx.setOut([]);
        end
        [obj, leth, status]= obj.shortest_minimal_cut_sets(numKOs, inf, inf, 2);
        leth_combin = [leth_combin, leth];
    end
end

%% 4. optimization
function [solution, status] = find_mcs(obj, obj_fun, timeLim,options)
    solution = nan;
    if timeLim<=0
        status = 1;
        return;
    elseif isinf(timeLim)
        timeLim = 1e75;
    end

    switch options.milp_solver
        case 'matlab_cplex'
            % set timelimit and remove MIP starts
            obj.obj_expr(obj.z_vars) = obj_fun;
            % Set new objective function
            obj.set_objective_function();
            obj.cpx.Param.timelimit.Cur = timeLim;
            obj.delete_all_mip_starts();
        case 'java_cplex'
            cp_inner = setup_cplex_inner_class_access();
            % remove old objective function
            if ~isempty(obj.cpx.getObjective()) , obj.cpx.remove(obj.cpx.getObjective()); 
            end
            % Set new objective function if provided
            if any(obj_fun)
                  obj.cpx.addMinimize().setExpr(obj.cpx.scalProd(obj_fun, obj.z_vars));
            end
            if ~isinf(timeLim),                   obj.cpx.setParam(cp_inner.DoubleParam.TiLim,timeLim);
            end
            for j = 0:(obj.cpx.getNMIPStarts-1),  obj.cpx.deleteMIPStarts(0);
            end
    end
    obj.solve();
    % If time limit is reached
    status = getMILPstatus(obj,options);
    switch status
        case 0 % If optimal
            solution = retrieve_optimal_solution(obj);
        case 1 % If time limit reached without solution
        case 2 % If infeasible
        case 3 % If feasible and time limit reached
            solution = retrieve_optimal_solution(obj);
        otherwise
            warning('Status unknown');
    end
    solution = round(solution);
    switch options.milp_solver
        case 'matlab_cplex'
            obj.clear_objective_function();
        case 'java_cplex'
            if ~isempty(obj.cpx.getObjective()) , obj.cpx.remove(obj.cpx.getObjective()); 
            end
    end
end

%% 5. fix z-vars to search in subspace of found cut set
function obj = set_ub_z_vars(obj,cs,options)
    switch options.milp_solver
        case 'matlab_cplex'
            obj.cpx.Model.ub(obj.z_vars) = cs;
        case 'java_cplex'
            arrayfun(@(x,y) x.setUB(y),obj.z_vars,cs);
    end
end

%% 6. return status depending on 
function status = getMILPstatus(obj,options)
    switch options.milp_solver
        case 'matlab_cplex'
            if obj.isOptimal() && ~obj.timeLimitReached() % completed
                status = 0;
            elseif obj.timeLimitReached() && ~obj.timeLimitReachedFeasible()
                status = 1;
            elseif obj.isInfeasible()
                status = 2;
            elseif obj.timeLimitReachedFeasible()
                status = 3;
            else
                status = -1;
            end
        case 'java_cplex'
            cp_inner = setup_cplex_inner_class_access();
            jStatus = obj.getStatus();
            CpxStatus = obj.getCplexStatus();
            if jStatus.equals(cp_inner.Status.Optimal) && ~CpxStatus.equals(cp_inner.CplexStatus.AbortTimeLim) % completed
                status = 0;
            elseif ~isinf(options.milp_time_limit) && CpxStatus.equals(cp_inner.CplexStatus.AbortTimeLim) % timelimit reached with or without solutions
                status = 1;
            elseif ~jStatus.equals(cp_inner.Status.Feasible) && ~CpxStatus.equals(cp_inner.CplexStatus.AbortTimeLim) % infeasible
                status = 2;
            elseif jStatus.equals(cp_inner.Status.Feasible) && CpxStatus.equals(cp_inner.CplexStatus.AbortTimeLim) % time limit reached, feasible
                status = 3;
            else
                status = -1;
            end
    end
end

%% 7. verify MCS
function [feas_T, feas_D] = verify_mcs(cnap,mcs,T,t,D,d,koCost,kiCost,verbose)
% this function is a shorter version of an externally available one.
    if license('test','Distrib_Computing_Toolbox') && isempty(getCurrentTask()) && (~isempty(ver('parallel'))  || ~isempty(ver('distcomp'))) 
        pool = gcp('nocreate');
        if ~isempty(pool) && pool.NumWorkers > 1
            parforarg = pool.NumWorkers;
        else
            parforarg = 0;
        end
    else
        parforarg = 0;
    end
    mcs = mcs.*~isnan(koCost') | (1-mcs).*~isnan(kiCost');
    feas_T = cell(size(mcs,2),1);
    feas_D = cell(size(mcs,2),1);
    parfor(i = 1:size(mcs,2),parforarg)
        cnap_valid = cnap;
        cnap_valid.reacMin(mcs(:,i)) = 0;
        cnap_valid.reacMax(mcs(:,i)) = 0;
        feas_D{i} = testRegionFeas(cnap_valid,D,d,2);
        feas_T{i} = testRegionFeas(cnap_valid,T,t,2);
        if any(feas_T{i})
            displ(['At least one target region (T' strrep(num2str(find(feas_T{i})'),'  ',',') ') is feasible in mutant model of mcs ' num2str(i)],verbose);
        end
        if any(~feas_D{i})
            displ(['At least one desired region (D' strrep(num2str(find(~feas_D{i})'),'  ',',') ') is infeasible in mutant model of mcs ' num2str(i)],verbose);
        end
    end
end

%% 8. text output of MCS
function ko_ki_text = mcs2text(mcs,reacID,koCost,kiCost)
    kos = find(mcs & ~isnan(koCost)');
    kis = find(mcs & ~isnan(kiCost)');
    if ~isempty(kis) 
        kis = [strjoin(strcat('+',cellstr(reacID(kis,:))),', ') ', '];
    else
        kis = [];
    end
    if ~isempty(kos) 
        kos = strjoin(strcat('-',cellstr(reacID(kos,:))),', ');
    else
        kos = [];
    end
    ko_ki_text = [kis kos];
end

%% 9. find redundant mcs
function redund = find_redundant_mcs(mcs)
% are there redundancies in a set of mcs
% return value is an integer vector that indicates whether the mcs is
% redundant in the set
% 0: true mcs
% 1: another mcs contains a subset of this mcs (smallest subset indicated)
% or the same mcs already occurred earlier
mcs = mcs ~= 0 & ~isnan(mcs);
mcs = sparse(mcs(sum(mcs,2)~=0,:));
redund = nan(1,size(mcs,2));
next = false(1,size(mcs,2));
next(1) = true;
while true
    i = next;
    unclear = find(isnan(redund) & ~i);
    s = relation(mcs(:,i),mcs(:,unclear));
    redund(unclear( s>=2 )) = 1; % superset or identical set
    if any(s == 1) % subset exists
        redund(i) = 1;
    else
        redund(i) = 0;
    end
    next(i) = false;
    next(find(isnan(redund),1)) = true;
    if ~any(next), break; end
end
function s = relation(mcs_ref,mcs_compare)
    contnset = all(mcs_compare( mcs_ref,:),1); % contained in the set
    diffset  = any(mcs_compare(~mcs_ref,:),1); % have nonzero values outside the set
    s = zeros(size(mcs_compare,2),1);
    s(~contnset & ~diffset) = 1; % subset
    s( contnset &  diffset) = 2; % superset
    s( contnset & ~diffset) = 3; % identical set
end
end

%% === utilities ===

% Suppress code warnings
%#ok<*ASGLU>
%#ok<*AGROW>

%% Spare functions
%% 1. Remap bounds with matrix (e.g. after compression)
function [splt_lb, splt_ub] = remapBounds(lb,ub,map)
    [splt_lb, splt_ub] = deal(nan(size(map,2),1));
    for k=1:size(map,2)
        forw=find(map(:,k)>0);
        revs=find(map(:,k)<0);
        if(~isempty(forw))
            splt_lb(k) = max(lb(forw)./map(forw,k));
            splt_ub(k) = min(ub(forw)./map(forw,k));
        end
        if(~isempty(revs))
            splt_lb(k) = max([splt_lb(k);ub(revs)./map(revs,k)]);
            splt_ub(k) = min([splt_ub(k);lb(revs)./map(revs,k)]);
        end
    end
end
