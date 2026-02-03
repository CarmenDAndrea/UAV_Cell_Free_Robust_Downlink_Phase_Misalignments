function D = MF_APClusterControl_Cmax(D,nbrOfUEs,gainOverNoise_dB, C_max)

for indexUE = 1:nbrOfUEs

    servingAPs = find(D(:,indexUE) == 1);

    if C_max < length(servingAPs)
        [~,indexAPs] = sort(gainOverNoise_dB(servingAPs,indexUE),'descend');
        indexAPs_keep_id = indexAPs(1:C_max);
        indexAPs_drop_id = indexAPs(C_max+1:end);
        
        D(servingAPs(indexAPs_keep_id),indexUE) = 1;
        D(servingAPs(indexAPs_drop_id),indexUE) = 0;

    else

    end

end
