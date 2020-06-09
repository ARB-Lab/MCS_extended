function [optFlux, success, status, optval]= CNAoptimizeFlux(cnap, fixedFluxes, c_macro, solver, verbose, pFBA, A_ieq, b_ieq)
%
% CellNetAnalyzer API function 'CNAoptimizeFlux'
% ---------------------------------------------
% --> flux optimization (flux balance analysis): the ojective function cnap.objFunc
%     is minimized by linear programming (LP)
%
% Usage: [optFlux, success, status, optval]= CNAoptimizeFlux(cnap, fixedFluxes, c_macro, solver, dispval,pFBA,A_ieq,b_ieq)
%
% cnap is a CellNetAnalyzer (mass-flow) project variable and mandatory argument.
%
% The other arguments are optional:
%
%   fixedFluxes: a (cnap.numr x 1) or empty vetor specifying predefined (fixed) fluxes.
%     If fixedFluxes(i) is a number then the flux through reaction i is fixed at this
%     given rate; if fixedFluxes(i) is NaN or if fixedFluxes is empty then the flux through
%     reaction i is free within the bounds given by cnap.reacMin(i) and cnap.reacMax(i).
%     Default:[]
%
%  c_macro: vector containing the concentrations (g/gDW) of the macromolecules
%     if a variable biomass composition has been defined (cnap.mue not empty).
%     Can be empty when cnap.mue or cnap.macroComposition is empty. If it
%     is empty and cnap.mue is not empty then cnap.macroDefault is used.
%     (default: cnap.macroDefault)
%
%   solver: selects the LP solver
%     0: GLPK (glpklp)
%     1: Matlab Optimization Toolbox (linprog)
%     2: CPLEX (cplexlp)
%     (default: 0)
%
%   verbose: controls the output printed to the console
%     -1: no output, even no warnings
%     0: no solver output, but warnings and information on final result will be shown
%     1: as option '0' but with additional solver output
%     (default: 0)
%
%   pFBA: flag that indicates whether a parsimonious FBA is to be done
%     0: no, a classical FBA is done
%     1: yes, the optimal flux vector with minimum sum of
%        absolute fluxes (secondary objective) is delivered
%
%   A_ieq: matrix of size (s x cnap.numr) for additional inequality constraints. May be empty.
%     default: []
%
%   b_ieq: vector of size (s x 1) for additional inequality constraints. May be empty.
%     default: []
%     A_ieq and b_ieq must either be both defined or both undefined (empty).
%     They define additional inequalities for the flux vector: A_ieq * r <= b_ieq
%     where r is the (optimal) flux vecor to be found.
%
%
% The following results are returned:
%
%   optFlux: if success == true: flux vector representing the optimal solution;
%            if succes == false: a flux vector with NaN entries indicating that
%            the result returned by the LP solver is not be meaningful
%
%   success: flag indicating whether the optimization was successful
%
%   status: solver status; for interpretation check the documentation of
%     the selected LP solver
%
%   optval: the optimum found for the objective function
%     (corresponds to cnap.objFunc'*optFlux)

% This file is part of CellNetAnalyzer. Please visit
% http://www.mpi-magdeburg.mpg.de/projects/cna/cna.html
% for more information and the latest version of CellNetAnalyzer.
%
% Copyright (C) 2000-2019 by Steffen Klamt and Axel von Kamp,
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

if(verLessThan('matlab','8.0'))
    error(nargchk(1, 8, nargin));
else
    narginchk(1, 8);
end

optFlux= nan(cnap.numr,1);
success= false;
status= NaN;
optval=NaN;

if(nargin<7)
    A_ieq=[];
end
if(nargin<8)
    b_ieq=[];
end
if ((size(A_ieq,1)~=size(b_ieq,1)) || (~isempty(A_ieq) && size(A_ieq,2) ~= cnap.numr))
    warning('CNAoptimizeFlux: Dimensions of A_ieq (s x cnap.numr) and/or of b_ieq (s x 1) are not consistent!')
    return;
end
if(nargin<6)
    pFBA=0;
end
if(nargin<5)
    verbose=0;
end
if(nargin<4)
    solver=0;
else
    if(~ismember(solver,[0,1,2]))
        error('Solver must be 0, 1, or 2.');
    end
    cnap.local.takelp= solver;
end
if(nargin>2 && ~isempty(c_macro))
    cnap.local.c_makro= c_macro;
else
    c_macro = cnap.macroDefault;
end

if nargin > 1 && ~isempty(fixedFluxes)
    fixedFluxes= reshape(fixedFluxes, length(fixedFluxes), 1);
else
    fixedFluxes=nan(cnap.numr,1);
end

LPavail=LP_solver_availability(true);
if(LPavail(solver+1)==false)
    solvers={'GLPK (glpk)','MATLAB (linprog)','CPLEX (cplexlp)'};
    error(['Solver ',solvers{solver+1},' not found. Please check whether you have porperly installed the toolbox and added the path!']);
end


%%%%%%%%%% New

truelyFixed=~isnan(fixedFluxes);
if ~all((fixedFluxes(truelyFixed) >= cnap.reacMin(truelyFixed)) & ... are the fixed fluxes in the
        (fixedFluxes(truelyFixed) <= cnap.reacMax(truelyFixed)))%      flux ranges of the model?
    displ('Infeasibility: the fixed fluxes do not lie within the boundaries defined in the model. Exit.',verbose>=0);
    return;
end

%% Prepare LP
N = initsmat(cnap.stoichMat,cnap.mue,cnap.macroComposition,c_macro,cnap.specInternal);
b = zeros(cnap.numis,1);
lb = cnap.reacMin;
ub = cnap.reacMax;
lb(truelyFixed) = fixedFluxes(truelyFixed);
ub(truelyFixed) = fixedFluxes(truelyFixed);
obj=cnap.objFunc;

switch solver
    case 0 % GLPK
        solvtext='glpk';
        ctype = repmat('S',1,cnap.numis);
        vtype= repmat('C', 1, cnap.numr);
        Nex=[N;A_ieq];
        bex=[b;b_ieq];
        ctype=[ctype,repmat('U', 1, length(b_ieq))];
        [optFlux, optval, status] = glpk(obj, Nex, bex, lb, ub, ctype,vtype);

    case 1 % linprog
        % opts = optimoptions('linprog','Display','off');
        % [optFlux,optval,status] = linprog(obj, [], [], N, b, lb, ub, linprog_options);
        solvtext='linprog';
        if verbose>=1
            try
                opts = optimoptions('linprog','Display','final');
            catch
                opts= optimset('Display', 'final');
            end
        else
            try
                opts = optimoptions('linprog','Display','off');
            catch
                opts= optimset('Display', 'off');
            end
        end

        [optFlux, optval, status] = linprog(obj, A_ieq, b_ieq, N, b, lb, ub,[],opts);

    case 2 % CPLEX
        solvtext='cplexlp';

        opts= cplexoptimset('cplex');
        opts.simplex.tolerances.feasibility = 1e-9;
        if verbose>=1
            if isfield(opts,'Display')
                opts.Display = 'on';
            elseif isfield(opts,'display')
                opts.display = 'on';
            elseif isfield(opts,'paramdisplay')
                opts.paramdisplay = 1;
            else
                opts= cplexoptimset('diagnostics', 'on');
            end
        else
            if isfield(opts,'Display')
                opts.Display = 'off';
            elseif isfield(opts,'display')
                opts.display = 'off';
            elseif isfield(opts,'paramdisplay')
                opts.paramdisplay = 0;
            else
                opts= cplexoptimset('diagnostics', 'off');
            end
        end

        [optFlux,optval,status] = cplexlp(obj, A_ieq, b_ieq, N, b, lb, ub, [],opts);

        % cgp= Cplex();
        % cgp.Param.emphasis.numerical.Cur= 1;
        % cgp.Param.simplex.tolerances.optimality.Cur= cgp1.Param.simplex.tolerances.optimality.Min;
        % cgp.Model.A= N;
        % cgp.Model.ub= ub;
        % cgp.Model.lb= lb;
        % cgp.Model.lhs= b;
        % cgp.Model.rhs= b;
        % cgp.Model.obj=obj';
        % cgp.DisplayFunc=[];
        % cgp.Model.sense='minimize';
        % x= cgp.solve();
        % status=x.status;
        % optval=x.objval;
end

if((solver == 0 && status == 5) || (solver == 1 && status == 1) || (solver == 2 && status == 1))  %optimal
    success = 1;
elseif ((solver == 0 && (status == 6  || status == 111)) || (solver == 1 &&  status == -3) || (solver == 2 && (status == 2  || status == 4  || status == -8 || status == -3)))  %% unbounded
    success = 0;
    displ(['Warning: solution seems to be unbounded. Status of solver ',solvtext,' : ',num2str(status)],verbose>=0);
elseif ((solver == 0 && (status == 3  || status == 110)) || (solver == 1 &&  status == -2) || (solver == 2 && (status == 3  || status == -2)))  %infeasible
    success = 0;
    displ(['Warning: problem seems to be infeasible. Status of solver ',solvtext,' : ',num2str(status)],verbose>=0);
else
    success = 0;
    displ(['Warning: optimal solution NOT found. Status of solver ',solvtext,' : ',num2str(status)],verbose>=0);
end

if ~success
    optFlux= nan(cnap.numr,1);
    optval=NaN;
    return;
else
    if(verbose>=0)
        disp('Optimization terminated.');
        resi=sum(abs(N*optFlux));
        if(resi>cnap.epsilon)
            disp(['Warning: Sum of residuals: ',num2str(resi)]);
        else
            disp(['Sum of residuals: ',num2str(resi)]);
        end
    end
end


if pFBA
    displ('Parsimonious FBA: Optimum found! Now minimizing sum of fluxes.',verbose>=0);

    % split network into forward and backward reactios
    cnapnew.stoichMat=[N -N];
    cnapnew.reacMax=[max(0,ub);max(0,-lb)];
    cnapnew.reacMin= [max(0,lb);max(0,-ub)];
    % minimize the total flux
    cnapnew.objFunc=ones(cnap.numr*2,1);
    cnapnew=CNAgenerateMFNetwork(cnapnew);

    %integrate original (extra) inequality constraints + constraint for optimal value
    Aobj=[A_ieq -A_ieq; obj' -obj'; -obj' obj'];
    bobj=[b_ieq; optval;-optval];

    [optpFBA, success, status, optvalpFBA]= CNAoptimizeFlux(cnapnew,[], [], solver,-1,0,Aobj,bobj);

    if ~success
        optFlux= nan(cnap.numr,1);
        success= false;
        status= NaN;
        optval= NaN;
        return;
    else
        optFlux=optpFBA(1:cnap.numr)-optpFBA(cnap.numr+1:end);
        if(verbose>=0)
            disp('Optimization terminated.');
            disp(['Result of parsimonious FBA: Sum of all fluxes: ',num2str(optvalpFBA)]);
            resi=sum(abs(N*optFlux));
            if(resi>cnap.epsilon)
                disp(['Warning: Sum of residuals: ',num2str(resi)]);
            else
                disp(['Sum of residuals: ',num2str(resi)]);
            end
        end
    end
end
