function [D,pilotIndex,masterAPs] = MF_APselectionBook(pilotIndex,nbrOfUEs,nbrOfAPs,masterAPs,mastersOfUEs,U_max,tau_p,gainOverNoise_dB)
nbrOfmasterAPs = ones(nbrOfUEs,1);
D = zeros(nbrOfAPs,nbrOfUEs);

% Master AP assignment
for indexUE = 1:nbrOfUEs
    % Coordination (master) AP
    [masterAPs,mastersOfUEs] = MF_mastersAPs(masterAPs,mastersOfUEs,nbrOfmasterAPs,U_max,indexUE,gainOverNoise_dB);
    D = masterAPs;

    % Pilot assignment
    pilotIndex = PilotAssignment(pilotIndex,indexUE,tau_p,gainOverNoise_dB,mastersOfUEs(indexUE));
end

% Non-master AP selection
for indexAP = 1:nbrOfAPs
    tmax = min(tau_p,nbrOfUEs);
    for t = 1:tmax
        % Find which UEs are using the pilot "t" in AP L.
        pilotUEs = find(t==pilotIndex(:,1));
        % If at least one UE is served by the AP in the pilot "t", this condition will not be satisfied
        if (sum(D(indexAP,pilotUEs)) == 0) && (sum(D(indexAP,:)) < U_max)
            % Assigns the pilot to best UE, the UEindex
            [~,UEindex] = max(gainOverNoise_dB(indexAP,pilotUEs));
            %Serve this UE if the channel is at most "threshold" weaker than the master AP's channel
            D(indexAP,pilotUEs(UEindex)) = 1;
        end
    end
end