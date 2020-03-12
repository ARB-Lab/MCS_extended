% This script reproduces the results from Figure 4 and 5. 
%
% In the example, a reaction-based setup is given together with a set of 
% different GPR rules for each reaction. CNAgeneMCSEnumerator2 is used to 
% compute gene-based MCS for this setup. The computation is supported by a
% GPR-rule compression and a network compression. After solving the MCS-MILP, 
% the MCS are expanded from the compressed to the full network. Although
% not used by the CNAgeneMCSEnumerator2, the function returns the full
% metabolic network with integrated GPR rules, to enable further analysis of 
% the computed MCS.
%
% Correspondence: cellnetanalyzer@mpi-magdeburg.mpg.de
% -Mar 2020

% Model Setup
cnap = struct();
cnap.specID = {'S';'A';'C';'F';'H';'D';'P';'Z';'BM'};
cnap.stoichMat = zeros(length(cnap.specID),0);
cnap = CNAgenerateMFNetwork(cnap);
cnap = CNAaddReactionMFN(cnap,'rs_up','= S' ,0,10);
cnap = CNAaddReactionMFN(cnap,'rd_ex','D =' ,0,100);
cnap = CNAaddReactionMFN(cnap,'rp_ex','P =' ,0,100);
cnap = CNAaddReactionMFN(cnap,'r_bm' ,'BM =',0,100,  -1); % objective
cnap = CNAaddReactionMFN(cnap,'r1','S = A'     ,0   ,100);
cnap = CNAaddReactionMFN(cnap,'r2','A + Z = BM',0   ,100);
cnap = CNAaddReactionMFN(cnap,'r3','S = P + Z' ,0   ,100);
cnap = CNAaddReactionMFN(cnap,'r4','S = C'     ,0   ,100);
cnap = CNAaddReactionMFN(cnap,'r5','C = H + Z' ,0   ,100);
cnap = CNAaddReactionMFN(cnap,'r6','H = D'     ,0   ,100);
cnap = CNAaddReactionMFN(cnap,'r7','A = F'     ,0   ,100);

dict_gpr_rules = {'r1' '(g1 and g2) or (g1 and g3)';
                  'r2' 'g2 or g4';
                  'r3' '(g3 and g6) or g8';
                  'r4' '(g4 and g8) or (g1 and g8)';
                  'r5' '(g1 and g5 and g9) or (g1 and g7)';
                  'r6' 'g4 and g1';
                  'r7' 'g6 or g7'};
for i = 1:size(dict_gpr_rules)
    idx = find(strcmp(cellstr(cnap.reacID),dict_gpr_rules{i,1}));
    cnap = CNAsetGenericReactionData(cnap,idx,'geneProductAssociation',dict_gpr_rules{i,2});
end

rd_ex = find(strcmp(cellstr(cnap.reacID),'rd_ex'));
r_bm  = find(strcmp(cellstr(cnap.reacID),'r_bm'));

% Target regions
T = full(sparse(  1, rd_ex , -1,1,cnap.numr));
t = -1;
% Desired region
D = full(sparse( 1,r_bm,-1,1,cnap.numr));
d = -1;

% Definition of knockable reactions
koable = find(~cellfun(@isempty,(regexp(cellstr(cnap.reacID),'r\d')))); % all numbered reactions are knockable
koCost = nan(cnap.numr,1);
koCost(koable) = 1;
% note that if only reaction KOs are defined and GPR rules are given, the algorithm
% will trace back targetable reactions to targetable genes. It is, however, possible
% to define up to 4 different vectors for addable and deletable genes and reactions
% rKoCost double<1,len(reactions)>: Costs for reaction  knock outs
% rKicost double<1,len(reactions)>: Costs for reaction "knock ins"
% gKoCost double<1,len(genes)>: Costs for gene knock outs
% gKiCost double<1,len(genes)>: Costs for gene knock ins
% A detailed explanation is given in the main funcion.

% MCS computation
[rmcs, gmcs, gcnap, cmp_gmcs, cmp_gcnap, mcs_idx] = CNAgeneMCSEnumerator2(cnap, {T} , {t} , {D} , {d} ,...
                                                    koCost,[], ... reackoCost,reackiCost
                                                    inf,cnap.numr,inf,... max_solutions,max_num_interv,time_limit
                                                    0, 2, [],[], ... use_bigM, enum_method, gkoCost, gkiCost
                                                    [],[0 1 1],1,1); % gpr_rules,use_compression,
% Additions are denoted with 1, deletions with -1 (here no solution contains deletions)
disp([cellstr(gcnap.reacID) num2cell(gmcs)]);
disp('The MCSs refer to the GPR-rule extended metabolic network. Knockouts are marked by -1.');