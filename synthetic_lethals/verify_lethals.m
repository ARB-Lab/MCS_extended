function [valid,max_mue] = verify_lethals(cnap,mcs)
    cnap_bu = cnap;
    valid   = nan(size(mcs,2),1);
    max_mue = nan(size(mcs,2),1);
    idx_mue = logical(cnap.objFunc);
    T = cnap.mcs.T{1};
    t = cnap.mcs.t{1};
    parfor i = 1:size(mcs,2)
        cnap = cnap_bu;
        cnap.reacMin(logical(mcs(:,i))) = 0;
        cnap.reacMax(logical(mcs(:,i))) = 0;
        fv = CNAoptimizeFlux(cnap,[],[],2,-1,0,T,t);
        if any(isnan(fv))
            valid(i) = 1;
            max_mue(i) = fv(idx_mue);
        else
            valid(i) = 0;
            fv = CNAoptimizeFlux(cnap,[],[],2,-1);
            max_mue(i) = fv(idx_mue);
        end
    end
end

