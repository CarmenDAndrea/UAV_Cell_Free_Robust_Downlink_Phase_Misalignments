% This function performs the Largest-large-scale-fading-based (LSFB) AP
% selection scheme
% Author: Marx Freitas

% INPUT:
% pilotIndex        = the pilot assigned to each UE
%                     Dim: nbrOfUEs x 1
% nbrOfUEs          = the number of UEs (K) in the newtork
% nbrOfAPs          = the number of APs (L) in the newtork
% gainOverNoise_dB  = the channel gain normalized by noise
%                     Dim: nbrOfAPs x nbrOfUEs
% tau_p             = the number of orthogonal pilot signals
% masterAPs         = the master AP of each UE
%                     Dim: nbrOfAPs x nbrOfUEs

% OUTPUT:
% D          = the AP cluster of each UE
%              Dim: nbrOfAPs x nbrOfUEs
% pilotIndex = the pilot assigned to each UE
%              Dim: nbrOfUEs x 1
% masterAPs  = the master AP of each UE
%              Dim: nbrOfAPs x nbrOfUEs

function [D, pilotIndex, masterAPs] = StrongestAPs(pilotIndex, nbrOfUEs, nbrOfAPs, gainOverNoise_dB, tau_p, C_max)

% Convert the channel gain from dB to linear
gainOverNoise_lin = db2pow(gainOverNoise_dB);

% Prepare to store the AP cluster
D = zeros(nbrOfAPs,nbrOfUEs);

%=========================================================================%
% ========================== LSFB AP selection ========================== %
%=========================================================================%
for indexUE = 1:nbrOfUEs
    
    % Sort the APs in desceding order
    [~, indexAP_sorted] = sort(gainOverNoise_lin(:,indexUE), 'descend');
        
    % Set the AP cluster
    D(indexAP_sorted(1:C_max),indexUE) = 1;

end

%=========================================================================%
% ========================== PILOT ASSIGNMENT =========================== %
%=========================================================================%
gainOverNoise_lin(D == 0) = 0; % It consider only the APs serving the UE

for indexUE = 1:nbrOfUEs

    % Identify the AP serving the UE presenting the strongest channel gain
    [~,strongestAP_id] = max(gainOverNoise_lin(:,indexUE));

    % Assign orthogonal pilots to the first tau_p UEs
    pilotIndex = PilotAssignment(pilotIndex,indexUE,tau_p,gainOverNoise_dB,strongestAP_id);

    % The AP with strongest channel gain is considered the "master" AP of
    % the UE. Note that there is no master AP in this scheme. This line is
    % used only to allow the code to run in different AP selection schemes
    masterAPs(strongestAP_id, indexUE) = 1;

end

% REFERENCES
% [3] H. Q. Ngo, L. Tran, T. Q. Duong, M. Matthaiou, and E. G. Larsson,
% "On the total energy efficiency of cell-free massive MIMO," IEEE Trans.
% Green Commun. Netw, vol. 2, no. 1, pp. 25–39, 2018.