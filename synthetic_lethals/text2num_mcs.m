function num_gmcs = text2num_mcs(gmcs,cnap)
    reac_names = cellstr(cnap.reacID);
    gmcs = cellfun(@(x) strcat('GP-',x),gmcs,'UniformOutput',0);
    num_gmcs = zeros(cnap.numr,size(gmcs,1));
    for i = 1:size(gmcs,1)
        num_gmcs(:,i) = -ismember(reac_names,gmcs{i});
    end
end

