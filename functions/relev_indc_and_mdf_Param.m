function [idx,mdfParam] = relev_indc_and_mdf_Param(cnap)
% function is used to find reaction and species indices that are used for
% the characterization and ranking of MCS
    % relevant reaction indices
    idx.o2      = find(~cellfun(@isempty ,regexp(cellstr(cnap.reacID),'.*EX_o2_e.*','match')));
    idx.atpm    = find(~cellfun(@isempty ,regexp(cellstr(cnap.reacID),'.*ATPM','match')));
    % relevant species indices
    idx.cytMet = find(~cellfun(@isempty ,regexp(cellstr(cnap.specID),'_c$','match')));
    % other important species
%     idx.pi    = find(strcmp(cellstr(cnap.specID), 'pi_c')); 
% (adding idx.pi activates the output of coupling mechanism analysis (not yet functional))
    idx.h     = find(strcmp(cellstr(cnap.specID), 'h_c'));
    idx.h2o   = find(strcmp(cellstr(cnap.specID), 'h2o_c'));
    idx.atp   = find(strcmp(cellstr(cnap.specID), 'atp_c'));
    idx.adp   = find(strcmp(cellstr(cnap.specID), 'adp_c'));
    idx.amp   = find(strcmp(cellstr(cnap.specID), 'amp_c'));
    idx.nad   = find(strcmp(cellstr(cnap.specID), 'nad_c'));
    idx.nadh  = find(strcmp(cellstr(cnap.specID), 'nadh_c'));
    idx.nadp  = find(strcmp(cellstr(cnap.specID), 'nadp_c'));
    idx.nadph = find(strcmp(cellstr(cnap.specID), 'nadph_c'));
    idx.co2_e = find(strcmp(cellstr(cnap.specID), 'co2_c'));
    idx.glc_e = find(strcmp(cellstr(cnap.specID), 'glc__D_e'));
    % MDF setup (thermodynamic benchmark)
    mdfParam.Cmin    = 1e-6*ones(cnap.nums,1);
    mdfParam.Cmin(idx.glc_e) = 1e-6;
    mdfParam.Cmax = 0.02*ones(cnap.nums,1);
    mdfParam.Cmax(idx.co2_e) = 1e-4;
    mdfParam.Cmax(idx.glc_e) = 0.055557;
    mdfParam.fixed_ratios(1,1:3) = [idx.atp   idx.adp   10];
    mdfParam.fixed_ratios(2,1:3) = [idx.adp   idx.amp    1];
    mdfParam.fixed_ratios(3,1:3) = [idx.nad   idx.nadh  10];
    mdfParam.fixed_ratios(4,1:3) = [idx.nadph idx.nadp  10];
    mdfParam.RT = 8.31446*300/1000; % Computation of MDF in kJ
    mdfParam.bottlenecks = 0; % change to 1 to compute thermodynamic bottlenecks
    mdfParam.G0 = cell2mat(CNAgetGenericReactionData_as_array(cnap,'deltaGR_0'));
    mdfParam.uncert = cell2mat(CNAgetGenericReactionData_as_array(cnap,'uncertGR_0'));
end