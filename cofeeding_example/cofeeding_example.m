% This script reproduces the results from Figure 2. The example demonstrates 
% how different strain designs that require a single substrate or co-feeding 
% can be found with a single MCS computation setup, using reaction additions 
% and multiple target regions.
%
% Correspondence: cellnetanalyzer@mpi-magdeburg.mpg.de
% -Jun 2020

% Model Setup
cnap = struct();
cnap.specID = {'S';'A';'B';'M1';'M2';'M3';'BM';'P';'Q';'R';'U';'ATP';'Z'};
cnap.stoichMat = zeros(length(cnap.specID),0);
cnap = CNAgenerateMFNetwork(cnap);
cnap = CNAaddReactionMFN(cnap,'rs_up','= S' ,0,10);
cnap = CNAaddReactionMFN(cnap,'ru_up','= U' ,0,10);
cnap = CNAaddReactionMFN(cnap,'ru_ex','U =' ,0,100);
cnap = CNAaddReactionMFN(cnap,'rp_ex','P =' ,0,100);
cnap = CNAaddReactionMFN(cnap,'rq_ex','Q =' ,0,100);
cnap = CNAaddReactionMFN(cnap,'rr_ex','R =' ,0,100);
cnap = CNAaddReactionMFN(cnap,'r_bm' ,'BM =',0,100,  -1); % objective
cnap = CNAaddReactionMFN(cnap,'r01','S = A + ATP'   ,0   ,100);
cnap = CNAaddReactionMFN(cnap,'r02','A = B + Z'     ,0   ,100);
cnap = CNAaddReactionMFN(cnap,'r03','A = B'         ,0   ,100);
cnap = CNAaddReactionMFN(cnap,'r04','S = B'         ,0   ,100);
cnap = CNAaddReactionMFN(cnap,'r05','Z + B = P'     ,0   ,100);
cnap = CNAaddReactionMFN(cnap,'r06','B = M1'        ,0   ,100);
cnap = CNAaddReactionMFN(cnap,'r07','Z + M1 = Q'    ,0   ,100);
cnap = CNAaddReactionMFN(cnap,'r08','M1 = M2'       ,0   ,100);
cnap = CNAaddReactionMFN(cnap,'r09','Z + M2 = R'    ,0   ,100);
cnap = CNAaddReactionMFN(cnap,'r10','M2 = M3'       ,0   ,100);
cnap = CNAaddReactionMFN(cnap,'r11','ATP + M3 = BM' ,0   ,100);
cnap = CNAaddReactionMFN(cnap,'r12','M3 = U'        ,-100,100);

rs_up = find(strcmp(cellstr(cnap.reacID),'rs_up'));
ru_up = find(strcmp(cellstr(cnap.reacID),'ru_up'));
ru_ex = find(strcmp(cellstr(cnap.reacID),'ru_ex'));
rp_ex = find(strcmp(cellstr(cnap.reacID),'rp_ex'));
r_bm  = find(strcmp(cellstr(cnap.reacID),'r_bm'));

% Target regions
Y = 0.4;
T1 = full(sparse(  [1      1      2     ], ... % case: single substrate
                   [rp_ex  rs_up  rs_up ], ...
                   [1      -Y     -1    ],2,cnap.numr));
t1 =  [  0 ; -1 ];
T2 = full(sparse(  [1     1      1      2      3     ], ... % case: co-feeding
                   [rp_ex rs_up  ru_up  ru_ex  rs_up ], ...
                   [1     -Y     -Y     1      -1    ],3,cnap.numr));
t2 =  [  0 ; 0 ; -1];
% Desired region
D = full(sparse( 1,r_bm,-1,1,cnap.numr));
d = -1;

% Definition of deletable, addable and non-targetable
koable = find(~cellfun(@isempty,(regexp(cellstr(cnap.reacID),'r\d\d')))); % all numbered reactions are knockable
kiable = ru_up;
koCost = nan(cnap.numr,1);
kiCost = nan(cnap.numr,1);
koCost(koable) = 1;
kiCost(kiable) = 1;

% MCS computation
mcs  = CNAMCSEnumerator2(cnap, {T1,T2}, {t1,t2}, D, d,koCost,kiCost);
disp([ {'MCS No.  ->'} num2cell(1:size(mcs,2)) ;[cellstr(cnap.reacID) num2cell(mcs)]]);
disp('Additions are marked with 1, deletions with -1, NaN marks addition candidates that were not added.');