% ----------------------------------------------------------------------- %
% AP SELECTION
% ----------------------------------------------------------------------- %
[D, pilotIndex, ~] = MF_APselectionSchemes(APselectionMethod, ...
    nbrOfAPs, nbrOfUEs, tau_p, U_max, gainOverNoise_dB, L_k_DSTBC);

% ----------------------------------------------------------------------- %
% AP CLUSTER CONTROL
% ----------------------------------------------------------------------- %
% Limit the maximum number of APs serving the user equipment (UE)
D = MF_APClusterControl_Cmax(D, nbrOfUEs, gainOverNoise_dB, L_k_DSTBC);

% ----------------------------------------------------------------------- %
% EXTEND THE AP CLUSTER UNTIL L_k_DSTBC
% ----------------------------------------------------------------------- %
% Ensure that the UE will be connected to L_k APs
D = extendAPclusterToCmin(D, L_k_DSTBC, gainOverNoise_dB, U_max);

% ----------------------------------------------------------------------- %
%  CHANNEL ESTIMATION
% ----------------------------------------------------------------------- %
% Set channel estimation method
channelestimationMethod = 'LMMSE';% 'phaseAwareMMSE';'LMMSE';'LS';

[Hhat, C] = functionChannelEstimatesWithRician(R, Hmean ,HmeanFase, HnoPM, ...
    nbrOfRealizations, nbrOfAPs, nbrOfUEs, N, tau_p, pilotIndex, ...
    poweUL_mW, channelestimationMethod, covariance_matrix, R_imp);

function D = extendAPclusterToCmin(D, C_max, gainOverNoise_dB, U_max)
% This function ensures that each UE is connected to C_max APs. That is, if
% the previous AP selection (or AP cluster adjustments) did not guarantee
% C_max connections to the APs, this function is activated. First, the
% function attempts to assign C_max APs to the UEs while respecting the
% scalability condition, i.e., U_max. If this is not possible, the function
% will override the scalability condition to provide more APs to the UE.
% The scalability condition will only be ignored if the number of UEs K,
% exceeds L * U_max, where L is the total number of APs in the network

% Author: Marx M. M. Freitas
% INPUT:
% D                 = the AP cluster of each UE
%                     Dim: nbrOfAPs x nbrOfUEs
% C_max             = the maximum number of APs serving the UE
% gainOverNoise_dB  = the channel gain normalized by noise
%                     Dim: nbrOfAPs x nbrOfUEs
% U_max             = the maximum number of UEs that each AP can serve

% OUTPUT
% D          = the AP cluster of each UE
%              Dim: nbrOfAPs x nbrOfUEs

% Compute the number of APs serving each UE, i.e., L_k
APsPerUE = sum(D, 1);
UEsWithSmallCmax = find(APsPerUE < C_max);

for id = 1:length(UEsWithSmallCmax)
    gainOverNoise_lin = db2pow(gainOverNoise_dB);

    % Identify which APs are already serving this UE
    servingAPs = find(D(:, UEsWithSmallCmax(id)) == 1);

    % Sort the APs by gain over noise (from strongest to weakest)
    [~, APs] = sort(gainOverNoise_lin(:, UEsWithSmallCmax(id)), 'descend');

    % Remove the APs that are already serving this UE
    remainingAPs = APs(~ismember(APs, servingAPs));

    % Number of APs that still need to be assigned
    Ek = C_max - length(servingAPs);

    % Available APs that have not yet reached U_max
    availableAPs = remainingAPs(sum(D(remainingAPs,:), 2) < U_max);

    % If there are available APs, we select the best ones among them
    if length(availableAPs) >= Ek
        strongestAPs = availableAPs(1:Ek);
    else
        % If there are not enough available APs, we select all the available ones
        strongestAPs = availableAPs;
        missing = Ek - length(strongestAPs);

        if missing > 0
            % We select already overloaded APs (ignoring U_max) to complete C_max
            overloadedAPs = remainingAPs(~ismember(remainingAPs, strongestAPs));
            strongestAPs = [strongestAPs; overloadedAPs(1:min(missing, length(overloadedAPs)))];
        end
    end

    % Update the D matrix to ensure L_k_DSTBC connections
    D(strongestAPs, UEsWithSmallCmax(id)) = 1;
end

end