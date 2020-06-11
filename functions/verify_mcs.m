function [feas_T, feas_D] = verify_mcs(cnap,mcs,T,t,D,d,koCost,kiCost,verbose)
% verify MCS
% When MCS are applied to a model, all sets of TARGET constraints should
% become INVALID and all DESIRED constraints should still be VALID.
%
% ---
% Input:
% ---
% cnap: CNA mass flow project
% mcs: <double>[numreacs x numMCS] vector or matrix that contains the
%      minimal cut sets. 
%      3 formats are possible:
%           (1) Only knockouts, KOs are indicated with "1" in the Matrix
%           (2) KOs and additions. KOs "-1", KIs "1" and non-added
%               candidates "nan"
%           (3) KOs and KIs are marked with "1". This requires a koCost and
%               kiCost vector to indicate which intervention is a KI and
%               which is a KO. If reaction is targetable, one of these
%               vectors carry a non-nan value at the reaction index.
% T,t,D,d: (cell arrays that contain double matrices: see MCSEnumerator2) 
%          Target and desired constraints
% koCost, kiCost: Only needed when mcs format is type 3. Non-nan values
%                 mark if reaction was KO- or KI-candidate
% verbose: <0,1> Flag. Text output when an MCS is invalid
%
% ---
% Output:
% ---
% feas_T: <cell>[numMCS x 1]{<double>[1 x numT]} Indicates feasible target
%         constrains for each MCS. 1: feasible, 0: infeasible.
%   >> special case << If only one output is queried, the function returns
%                      whether the MCSs are valid (1) (all targets infeasible 
%                      and all desired feasible) or invalid (0).
% feas_D: <cell>[numMCS x 1]{<double>[1 x numD]} Indicates feasible desired
%         constrains for each MCS. 1: feasible, 0: infeasible
% 
    if nargin<9
        verbose = 0;
    end
    if nargin>6 && ~isempty(koCost) && ~isempty(kiCost) % Type 3 MCS
        mcs = mcs.*~isnan(koCost') | (1-mcs).*~isnan(kiCost');
    elseif any(any(mcs == -1)) % Type 2 MCS
        mcs = isnan(mcs) | mcs==-1;
    end
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
    mcs = logical(mcs);
    feas_T = cell(size(mcs,2),1);
    feas_D = cell(size(mcs,2),1);
    parfor(i = 1:size(mcs,2),parforarg)
        cnap_valid = cnap;
        cnap_valid.reacMin(mcs(:,i)) = 0;
        cnap_valid.reacMax(mcs(:,i)) = 0;
		feas_D{i} = testRegionFeas(cnap,D,d,2);
		feas_T{i} = testRegionFeas(cnap,T,t,2);
        if verbose
            if any(feas_T{i})
                disp(['At least one target region (T' strrep(num2str(find(feas_T{i})'),'  ',',') ') is feasible in mutant model of mcs ' num2str(i)]);
            end
            if any(~feas_D{i})
                disp(['At least one desired region (D' strrep(num2str(find(~feas_D{i})'),'  ',',') ') is infeasible in mutant model of mcs ' num2str(i)]);
            end
        end
    end
    if nargout == 1
        feas_T = all(~cell2mat(feas_T),2) & all(cell2mat(feas_D),2);
    end
end