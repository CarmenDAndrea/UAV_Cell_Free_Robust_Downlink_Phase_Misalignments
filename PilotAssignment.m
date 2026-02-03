%Pilot assignment
function pilotIndex = PilotAssignment(pilotIndex,indexUE,tau_p,gainOverNoise_dB,indexSortedAPs)
gainOverNoise_Lin = db2pow(gainOverNoise_dB);

if indexUE <= tau_p
    pilotIndex(indexUE) = indexUE;
else %Assign pilot for remaining UEs
    %Compute received power from to the master AP from each pilot
    pilotinterference = zeros(tau_p,1);    
    for t = 1:tau_p
        pilotinterference(t) = sum(gainOverNoise_Lin(indexSortedAPs(end), pilotIndex(1:indexUE-1)==t));
    end
    %Find the pilot with the least receiver power
    [~,bestpilot] = min(pilotinterference);
    pilotIndex(indexUE) = bestpilot;    
end