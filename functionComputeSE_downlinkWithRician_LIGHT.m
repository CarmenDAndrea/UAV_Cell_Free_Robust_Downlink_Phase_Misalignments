% function [SE_P_MMSE, SE_P_RZF, SE_LP_MMSE, SE_MR, w_MR, w_LP_MMSE, w_P_MMSE, w_PRZF, w_MMSE]...
%     = functionComputeSE_downlinkWithRician_LIGHT(Hhat,H,D,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,...
%     rho_dist_mW, gainOverNoisedB, PowerDL_mW)

% function [w_MR, w_LP_MMSE, w_P_MMSE, w_PRZF, w_MMSE, b_VLSFB_MR, b_VLSFB_LP_MMSE]...
%     = functionComputeSE_downlinkWithRician_LIGHT(Hhat,H,D,C,tau_c,tau_p,nbrOfRealizations,N,K,L,p,...
%     rho_dist_mW, gainOverNoisedB, PowerDL_mW)

function [w_MR, w_LP_MMSE, w_P_MMSE, w_P_RZF, b_VLSFB_MR, b_VLSFB_LP_MMSE]...
    = functionComputeSE_downlinkWithRician_LIGHT(Hhat, H, D, C, tau_c, tau_p, nbrOfRealizations, N, K, L, p,...
    rho_dist_mW, gainOverNoisedB, PowerDL_mW)

% Compute downlink SE for different transmit precoding schemes using the
% capacity bound in Theorem 6.1 for the centralized schemes and the
% capacity bound in Corollary 6.3 for the distributed schemes. Compute the
% genie-aided downlink SE from Corollary 6.6 for the centralized and the
% distributed operations.

%
% INPUT:
% Hhat              = Matrix with dimension L*N  x nbrOfRealizations x K
%                    where (:,n,k) is the estimated collective channel to
%                    UE k in channel realization n.
% H                 = Matrix with dimension L*N  x nbrOfRealizations x K
%                     with the true channel realizations. The matrix is
%                     organized in the same way as Hhat.
% D                 = DCC matrix for cell-free setup with dimension L x K
%                     where (l,k) is one if AP l serves UE k and zero otherwise
% B                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                     the spatial correlation matrix of the channel estimate
%                     between AP l and UE k, normalized by noise variance
% C                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                     the spatial correlation matrix of the channel
%                     estimation error between AP l and UE k,
%                    normalized by noise variance
% tau_c             = Length of coherence block
% tau_p             = Length of pilot sequences
% nbrOfRealizations = Number of channel realizations
% N                 = Number of antennas per AP
% K                 = Number of UEs
% L                 = Number of APs
% p                 = Uplink transmit power per UE (same for everyone)
% R                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                     the spatial correlation matrix between AP l and UE k,
%                     normalized by noise
% pilotIndex        = Vector containing the pilot assigned to each UE
% rho_dist          = Matrix with dimension L x K where (l,k) is the power
%                     allocated to UE k by AP l in the distributed downlink
%                     operation
% gainOverNoisedB   = Matrix with dimension L x K where (l,k) is the channel
%                     gain (normalized by the noise variance) between AP l
%                     and UE k
% rho_tot           = Maximum allowed transmit power for each AP
%
% OUTPUT:
% SE_P_MMSE         = SEs achieved with P-MMSE precoding in (6.17)
% SE_P_RZF          = SEs achieved with P-RZF precoding in (6.18)
% SE_LP_MMSE        = SEs achieved with LP-MMSE precoding in (6.33)
% SE_MR             = SEs achieved with MR precoding in (6.26)

% w_MR              = The precoding matrix of MR scheme.
%                     Dim: L x K x nbrOfRealizations
% w_LP_MMSE         = The precoding matrix of LP-MMSE scheme.
%                     Dim: L x K x nbrOfRealizations
% w_P_MMSE          = The precoding matrix of P-MMSE scheme.
%                     Dim: L x K x nbrOfRealizations
% w_PRZF            = The precoding matrix of P-RZF scheme.
%                     Dim: L x K x nbrOfRealizations


% This Matlab function was developed to generate simulation results to:
%
% Ozlem Tugfe Demir, Emil Bjornson and Luca Sanguinetti (2021),
% "Foundations of User-Centric Cell-Free Massive MIMO",
% Foundations and Trends in Signal Processing: Vol. 14: No. 3-4,
% pp 162-472. DOI: 10.1561/2000000109
%
% This is version 1.0 (Last edited: 2021-01-31)
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% monograph as described above.

% Some comments: Marx Freitas

% Store the N x N identity matrix
eyeN = eye(N);

% Compute the prelog factor assuming only downlink data transmission
% prelogFactor = (1-tau_p/tau_c);

% Prepare to store the terms that appear in SEs
scaling_MR = zeros(L,K);
% interUserGains_MR = zeros(K,K,nbrOfRealizations);

% signal_MR1 = zeros(K,1);
% interf_MR1 = zeros(K,1);

% signal_P_MMSE = zeros(K,1);
% interf_P_MMSE = zeros(K,1);
scaling_P_MMSE = zeros(K,1);
portionScaling_P_MMSE = zeros(L,K);
% interUserGains_P_MMSE = zeros(K,K,nbrOfRealizations);

% signal_MMSE = zeros(K,1);
% interf_MMSE = zeros(K,1);
% scaling_MMSE = zeros(K,1);
% portionScaling_MMSE = zeros(L,K);
% interUserGains_MMSE = zeros(K,K,nbrOfRealizations);

% signal_P_RZF = zeros(K,1);
% interf_P_RZF = zeros(K,1);
scaling_P_RZF = zeros(K,1);
portionScaling_PRZF = zeros(L,K);
% interUserGains_P_RZF = zeros(K,K,nbrOfRealizations);

% signal_LP_MMSE = zeros(K,1);
% interf_LP_MMSE = zeros(K,1);
scaling_LP_MMSE = zeros(L,K);
% interUserGains_LP_MMSE = zeros(K,K,nbrOfRealizations);

% MF editing
% Prepare to store the power consumed in the APs in each precoding scheme
w_MR      = zeros(L*N,K,nbrOfRealizations);
w_LP_MMSE = zeros(L*N,K,nbrOfRealizations);
w_P_MMSE  = zeros(L*N,K,nbrOfRealizations);
w_P_RZF    = zeros(L*N,K,nbrOfRealizations);
% w_MMSE    = zeros(L*N,K,nbrOfRealizations);

%%%% w_MR_hat = zeros(L*N,K,nbrOfRealizations);
%%%% w_LP_MMSE_hat = zeros(L*N,K,nbrOfRealizations);

% Prepare to store VLSFP results
signal_LSFP_LP_MMSE = zeros(K,K,L);
signal2_LSFP_LP_MMSE = zeros(K,K,L);
scaling_LSFP_LP_MMSE = zeros(L,K);

signal_LSFP_MR  = zeros(K,K,L);
signal2_LSFP_MR = zeros(K,K,L);
scaling_LSFP_MR = zeros(L,K);


%%
% ----------------------------------------------------------------------- %
% COMBINING AND SCALING FACTORS %
% ----------------------------------------------------------------------- %
% First generate the combining, then the scaling factors
%Go through all channel realizations
for n = 1:nbrOfRealizations
    % ------------------------------------------------------------------- %
    % MR and LP-MMSE schemes
    % ------------------------------------------------------------------- %
    % Go through all APs
    for l = 1:L
        %Extract channel realizations from all UEs to AP l
        Hallj = reshape(H(1+(l-1)*N:l*N,n,:),[N K]);

        % Extract which UEs are served by the AP
        servedUEs = find(D(l,:)==1);

        % Compute sum of estimation error covariance matrices of the UEs served by AP l
        Cserved = sum(C(:,:,l,servedUEs),4);

        % Compute MR and LP-MMSE combining
        V_MR = reshape(Hhat((l-1)*N+1:l*N,n,:),[N K]);
        V_LP_MMSE = p*((p*(V_MR(:,servedUEs)*V_MR(:,servedUEs)')+p*Cserved+eyeN)\V_MR(:,servedUEs));

        % Compute the precoding scaling factor by Monte Carlo methods
        scaling_MR(l,servedUEs) = scaling_MR(l,servedUEs) + sum(abs(V_MR(:,servedUEs)).^2,1)/nbrOfRealizations;
        scaling_LP_MMSE(l,servedUEs) = scaling_LP_MMSE(l,servedUEs) + sum(abs(V_LP_MMSE).^2,1)/nbrOfRealizations;

        for ind = 1:length(servedUEs)
            %Extract UE index
            k = servedUEs(ind);

            % Signal and scaling factor for VLSFP, MR
            signal_LSFP_MR(:, k, l)  = signal_LSFP_MR(:,k,l)  + Hallj'*V_MR(:, ind)/nbrOfRealizations;
            signal2_LSFP_MR(:, k, l) = signal2_LSFP_MR(:,k,l) + abs(Hallj'*V_MR(:, ind)).^2/nbrOfRealizations;
            scaling_LSFP_MR(l,k)     = scaling_LSFP_MR(l,k)   + sum(abs(V_MR(:, ind)).^2,1)/nbrOfRealizations;

            % Signal and scaling factor for VLSFP, LP_MMSE
            signal_LSFP_LP_MMSE(:, k, l)  = signal_LSFP_LP_MMSE(:, k, l)  + Hallj'*V_LP_MMSE(:, ind)/nbrOfRealizations;
            signal2_LSFP_LP_MMSE(:, k, l) = signal2_LSFP_LP_MMSE(:, k, l) + abs(Hallj'*V_LP_MMSE(:, ind)).^2/nbrOfRealizations;
            scaling_LSFP_LP_MMSE(l, k)    = scaling_LSFP_LP_MMSE(l, k)    + sum(abs(V_LP_MMSE(:, ind)).^2,1)/nbrOfRealizations;
        end
    end

    % Final do loop dos APs (L)
    clear V_MR V_LP_MMSE Cserved;

    % ------------------------------------------------------------------- %
    % P-MMSE and P-RZF schemes
    % ------------------------------------------------------------------- %
    % Go through all UEs
    for k = 1:K
        % Determine the set of serving APs
        servingAPs = find(D(:,k)==1);
        La = length(servingAPs);

        % Determine which UEs that are served by partially the same set
        % of APs as UE k, i.e., the set in (5.15)
        servedUEs = sum(D(servingAPs,:),1)>=1;

        % Extract channel realizations and estimation error correlation
        % matrices for the APs that involved in the service of UE k
        Hhatallj_active = zeros(N*La,K);
        C_tot_blk = zeros(N*La,N*La);
        C_tot_blk_partial = zeros(N*La,N*La);

        for l = 1:La
            Hhatallj_active((l-1)*N+1:l*N,:) = reshape(Hhat((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            C_tot_blk((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),:),4);
            C_tot_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),servedUEs),4);
        end

        % Compute MMSE, P-MMSE and P-RZF combining
        V_P_MMSE = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+p*C_tot_blk_partial+eye(La*N))\Hhatallj_active(:,k));
        % V_MMSE = p*((p*(Hhatallj_active*Hhatallj_active')+p*C_tot_blk+eye(La*N))\Hhatallj_active(:,k));
        V_P_RZF = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+eye(La*N))\Hhatallj_active(:,k));

        % Compute the precoding scaling factor by Monte Carlo methods
        scaling_P_MMSE(k) = scaling_P_MMSE(k) + sum(abs(V_P_MMSE).^2,1)/nbrOfRealizations;
        % scaling_MMSE(k) = scaling_MMSE(k) + sum(abs(V_MMSE).^2,1)/nbrOfRealizations;
        scaling_P_RZF(k) = scaling_P_RZF(k) + sum(abs(V_P_RZF).^2,1)/nbrOfRealizations;

        % Go through all the serving APs
        for l=1:La
            % Auxiliar variables for the combining vectors
            V_P_MMSE2 = V_P_MMSE((l-1)*N+1:l*N,:);
            % V_MMSE2 = V_MMSE((l-1)*N+1:l*N,:);
            V_P_RZF2 = V_P_RZF((l-1)*N+1:l*N,:);

            % Portiong scaling factors
            portionScaling_P_MMSE(servingAPs(l),k) = portionScaling_P_MMSE(servingAPs(l),k) + sum(abs(V_P_MMSE2).^2,1)/nbrOfRealizations;
            % portionScaling_MMSE(servingAPs(l),k)   = portionScaling_MMSE(servingAPs(l),k) + sum(abs(V_MMSE2).^2,1)/nbrOfRealizations;
            portionScaling_PRZF(servingAPs(l),k)   = portionScaling_PRZF(servingAPs(l),k) + sum(abs(V_P_RZF2).^2,1)/nbrOfRealizations;
        end
    end

    % Final do loop dos UEs (K)
    clear Hhatallj_active C_tot_blk C_tot_blk_partial;
    clear V_P_MMSE V_MMSE V_P_RZF;
end

% Normalize the norm squares of the portions for the normalized
% centralized precoders
portionScaling_P_MMSE = portionScaling_P_MMSE./repmat(scaling_P_MMSE.',[L 1]);
% portionScaling_MMSE   = portionScaling_MMSE./repmat(scaling_MMSE.',[L 1]);
portionScaling_PRZF   = portionScaling_PRZF./repmat(scaling_P_RZF.',[L 1]);

% ----------------------------------------------------------------------- %
% CENTRALIZED POWER ALLOCATION
% ----------------------------------------------------------------------- %
% The parameters for the centralized downlink power allocation in (7.43)
upsilon = -0.5; kappa = 0.5; theta = 0.2;

% Compute the power allocation coefficients for centralized precoding
% according to (7.43), fractional power allocation
rho_P_MMSE = functionCentralizedPowerAllocation(K, gainOverNoisedB, D, PowerDL_mW, portionScaling_P_MMSE, upsilon, kappa, theta);
% rho_MMSE   = functionCentralizedPowerAllocation(K, gainOverNoisedB, D, PowerDL_mW, portionScaling_MMSE, upsilon, kappa, theta);
rho_PRZF   = functionCentralizedPowerAllocation(K, gainOverNoisedB, D, PowerDL_mW,  portionScaling_PRZF, upsilon, kappa, theta);

% Compute power allocation coefficients for VLSFP
b_VLSFB_MR = vlsfp(L, K, D, signal_LSFP_MR, signal2_LSFP_MR, scaling_LSFP_MR);
b_VLSFB_LP_MMSE = vlsfp(L, K, D, signal_LSFP_LP_MMSE, signal2_LSFP_LP_MMSE, scaling_LSFP_LP_MMSE);

% [b_MR, ~] = functionCentralizedPowerAllocationVLSFP(L, K, gainOverNoisedB, D, PowerDL_mW, b_VLSFB_MR, upsilon, kappa, theta);
% [b_LP_MMSE, ~] = functionCentralizedPowerAllocationVLSFP(L, K, gainOverNoisedB, D, PowerDL_mW, b_VLSFB_LP_MMSE, upsilon, kappa, theta);

clear gainOverNoisedB PowerDL_mW upsilon kappa
clear portionScaling_P_MMSE portionScaling_MMSE;

%%
% ----------------------------------------------------------------------- %
%  PRECODING VECTORS
% ----------------------------------------------------------------------- %
% Generates the precoding from combining vectors

% Go through all channel realizations
for n = 1:nbrOfRealizations
    % Matrix to store Monte-Carlo results for this realization
    % interf_MR_n = zeros(K,K);
    % interf_P_MMSE_n = zeros(K,K);
    % interf_MMSE_n = zeros(K,K);
    % interf_P_RZF_n = zeros(K,K);
    % interf_LP_MMSE_n = zeros(K,K);

    % ------------------------------------------------------------------- %
    % MR and LP-MMSE schemes
    % ------------------------------------------------------------------- %
    % Go through all APs
    for l = 1:L

        % Extract channel realizations from all UEs to AP l
        Hallj = reshape(H((l-1)*N+1:l*N,n,:),[N K]);

        % Extract channel estimates from all UEs to AP l
        Hhatallj = reshape(Hhat((l-1)*N+1:l*N,n,:),[N K]);

        % Extract which UEs are served by AP l
        servedUEs = find(D(l,:)==1);

        % Compute sum of estimation error covariance matrices of the UEs
        % served by AP l
        Cserved = sum(C(:,:,l,servedUEs),4);

        % Compute MR combining
        V_MR = Hhatallj(:,servedUEs);

        % Compute LP-MMSE combining
        V_LP_MMSE = p*((p*(V_MR*V_MR')+p*Cserved+eyeN)\V_MR);

        % Go through all UEs served by the AP
        for ind = 1:length(servedUEs)

            % Extract UE index
            k = servedUEs(ind);

            % Normalize MR precoding
            w = V_MR(:,ind)*sqrt(rho_dist_mW(l,k)/scaling_MR(l,k));
            w_MR((l-1)*N+1:l*N,k,n) = w;
            %%%%% w_MR_hat((l-1)*N+1:l*N,k,n) = V_MR(:,ind)*sqrt(1/scaling_MR(l,k));

            % Compute gain of the signal from UE that arrives at other UEs
            % interUserGains_MR(:,k,n) = interUserGains_MR(:,k,n) + Hallj'*w;

            % signal_MR1(k) = signal_MR1(k) + (Hallj(:,k)'*w)/nbrOfRealizations;
            % interf_MR_n(:,k) = interf_MR_n(:,k) + Hallj'*w;

            % Normalize LP-MMSE precoding
            w = V_LP_MMSE(:,ind)*sqrt(rho_dist_mW(l,k)/scaling_LP_MMSE(l,k));
            w_LP_MMSE((l-1)*N+1:l*N,k,n) = w;
            %%%%% w_LP_MMSE_hat((l-1)*N+1:l*N,k,n) = V_LP_MMSE(:,ind)*sqrt(1/scaling_LP_MMSE(l,k));

            % Compute realizations of the terms inside the expectations
            % of the signal and interference terms in Corollary 6.3

            % signal_LP_MMSE(k) = signal_LP_MMSE(k) + (Hallj(:,k)'*w)/nbrOfRealizations;
            % interf_LP_MMSE_n(:,k) = interf_LP_MMSE_n(:,k) + Hallj'*w;

            % Compute gain of the signal from UE that arrives at other UEs
            % interUserGains_LP_MMSE(:,k,n) = interUserGains_LP_MMSE(:,k,n) + Hallj'*w;
        end
    end

    % Final do loop dos APs
    clear Hallj Hhatallj V_MR V_LP_MMSE;

    % ------------------------------------------------------------------- %
    % P-RZF and P-MMSE schemes
    % ------------------------------------------------------------------- %
    % Consider the centralized schemes
    % Go through all UEs
    for k = 1:K

        % Determine the set of serving APs
        servingAPs = find(D(:,k)==1);

        La = length(servingAPs);

        % Determine which UEs that are served by partially the same set
        % of APs as UE k, i.e., the set in (5.15)
        servedUEs = sum(D(servingAPs,:),1)>=1;

        % Extract channel realizations and estimation error correlation
        % matrices for the APs that involved in the service of UE k
        Hallj_active = zeros(N*La,K);

        Hhatallj_active = zeros(N*La,K);
        C_tot_blk = zeros(N*La,N*La);
        C_tot_blk_partial = zeros(N*La,N*La);

        for l = 1:La
            Hallj_active((l-1)*N+1:l*N,:) = reshape(H((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            Hhatallj_active((l-1)*N+1:l*N,:) = reshape(Hhat((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            C_tot_blk((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),:),4);
            C_tot_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),servedUEs),4);
        end

        % Compute P-MMSE precoding
        w = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+p*C_tot_blk_partial+eye(La*N))\Hhatallj_active(:,k));

        % Apply power allocation
        w = w*sqrt(rho_P_MMSE(k)/scaling_P_MMSE(k));

        % MF Modification
        for index = 1:length(servingAPs)
            l = servingAPs(index);
            w_P_MMSE(((l-1)*N+1):(l*N),k,n) = w((index-1)*N+1:index*N);
        end

        % Compute realizations of the terms inside the expectations
        % of the signal and interference terms in Theorem 6.1
        % signal_P_MMSE(k) = signal_P_MMSE(k) + (Hallj_active(:,k)'*w)/nbrOfRealizations;
        % interf_P_MMSE_n(:,k) = interf_P_MMSE_n(:,k) + Hallj_active'*w;

        % Compute gain of the signal from UE that arrives at other UEs
        % interUserGains_P_MMSE(:,k,n) = interUserGains_P_MMSE(:,k,n) + Hallj_active'*w;

        %%Compute MMSE combining
        % w = p*((p*(Hhatallj_active*Hhatallj_active')+p*C_tot_blk+eye(La*N))\Hhatallj_active(:,k));

        %% Apply power allocation
        % w = w*sqrt(rho_MMSE(k)/scaling_MMSE(k));

        % %MF Modification
        % for index = 1:length(servingAPs)
        %     l = servingAPs(index);
        %     w_MMSE(((l-1)*N+1):(l*N),k,n) = w((index-1)*N+1:index*N);
        % end

        %Compute realizations of the terms inside the expectations
        %of the signal and interference terms in Theorem 6.1
        % signal_MMSE(k) = signal_MMSE(k) + (Hallj_active(:,k)'*w)/nbrOfRealizations;
        % interf_MMSE_n(:,k) = interf_MMSE_n(:,k) + Hallj_active'*w;

        %Compute gain of the signal from UE that arrives at other UEs
        % interUserGains_MMSE(:,k,n) = interUserGains_MMSE(:,k,n) + Hallj_active'*w;

        % Compute P-RZF combining
        w = p*((p*(Hhatallj_active(:,servedUEs)*Hhatallj_active(:,servedUEs)')+eye(La*N))\Hhatallj_active(:,k));

        % Apply power allocation
        w = w*sqrt(rho_PRZF(k)/scaling_P_RZF(k));

        % MF Modification
        for index = 1:length(servingAPs)
            l = servingAPs(index);
            w_P_RZF(((l-1)*N+1):(l*N),k,n) = w((index-1)*N+1:index*N);
        end

        % Compute realizations of the terms inside the expectations
        % of the signal and interference terms in Theorem 6.1
        % signal_P_RZF(k) = signal_P_RZF(k) + (Hallj_active(:,k)'*w)/nbrOfRealizations;
        % interf_P_RZF_n(:,k) = interf_P_RZF_n(:,k) + Hallj_active'*w;

        % % Compute gain of the signal from UE that arrives at other UEs
        % interUserGains_P_RZF(:,k,n) = interUserGains_P_RZF(:,k,n) + Hallj_active'*w;

    end

    % Final do loop dos UEs
    clear Hallj_active Hhatallj_active C_tot_blk C_tot_blk_partial;
    clear w servedUEs servingAPs;


    % Compute interference power in one realization
    % interf_P_MMSE = interf_P_MMSE + sum(abs(interf_P_MMSE_n).^2,2)/nbrOfRealizations;
    % interf_MMSE = interf_MMSE + sum(abs(interf_MMSE_n).^2,2)/nbrOfRealizations;
    % interf_P_RZF = interf_P_RZF + sum(abs(interf_P_RZF_n).^2,2)/nbrOfRealizations;

    % interf_MR1 = interf_MR1 + sum(abs(interf_MR_n).^2,2)/nbrOfRealizations;
    % interf_LP_MMSE = interf_LP_MMSE + sum(abs(interf_LP_MMSE_n).^2,2)/nbrOfRealizations;

end


clear interUserGains_MR interUserGains_LP_MMSE interUserGains_P_MMSE interUserGains_P_RZF interUserGains_MMSE;
clear rho_P_MMSE rho_MMSE rho_PRZF H_hat H rho_dist_mW;
clear scaling_MR scaling_LP_MMSE scaling_P_MMSE scaling_MMSE scaling_P_RZF;

end


% LOCAL FUNCTIONS
function b_vlsfp = vlsfp(L, K, D, signal, signal2, scaling)

%Prepare to save simulation results
b_vlsfp = zeros(L,K);

%Prepare to store the interference matrix in (7.26)
Ck = zeros(L,L,K,K);

%Go through all UEs
for k = 1:K
    %Find the APs that serve UE k
    servingAPs = find(D(:,k)==1);
    %The number of APs that serve UE k
    La = length(servingAPs);
    %Compute the vector in (7.25) for UE k (only the non-zero indices correspondig to
    %serving APs are considered)
    fkk = conj(vec(signal(k,k,servingAPs)))./sqrt(scaling(servingAPs,k));

    CCk = eye(La);
    %Go through all UEs
    for i = 1:K
        %Compute the matrices
        if i==k
            Cik = fkk*fkk';
        else
            Cik = diag(1./sqrt(scaling(servingAPs,k)))...
                *(conj(vec(signal(i,k,servingAPs)))...
                *conj(vec(signal(i,k,servingAPs)))')...
                *diag(1./sqrt(scaling(servingAPs,k)));
        end

        for j = 1:La
            Cik(j,j) = signal2(i,k,servingAPs(j))/scaling(servingAPs(j),k);
        end

        Ck(1:La,1:La,i,k) = Cik;

        CCk = CCk + Cik;
    end

    %Compute normalized virtual LSFD vectors for UE k
    b_vlsfp(servingAPs,k) =  CCk\fkk;
    b_vlsfp(:,k) = b_vlsfp(:,k)/norm(b_vlsfp(:,k));
end

end



% % ======================================================================= %
% % =================== DOWNLINK SPECTRAL EFFICIENCY ====================== %
% % ======================================================================= %
%
% % Compute SE in Corollary 6.3 with MR  using the closed-form expressions in Corollary 6.4
% SE_MR = prelogFactor*real(log2(1+(abs(signal_MR1).^2) ./ (interf_MR1 - abs(signal_MR1).^2 + 1)));
%
% % Compute SE  in Corollary 6.3 with LP-MMSE
% SE_LP_MMSE = prelogFactor*real(log2(1+(abs(signal_LP_MMSE).^2) ./ (interf_LP_MMSE - abs(signal_LP_MMSE).^2 + 1)));
%
% % Compute SE in Theorem 6.1 with P-MMSE
% SE_P_MMSE = prelogFactor*real(log2(1+(abs(signal_P_MMSE).^2) ./ (interf_P_MMSE - abs(signal_P_MMSE).^2 + 1)));
%
% % Compute SE in Theorem 6.1 with P-RZF
% SE_P_RZF = prelogFactor*real(log2(1+(abs(signal_P_RZF).^2) ./ (interf_P_RZF - abs(signal_P_RZF).^2 + 1)));

function [b_k, rho_vlsfp] = functionCentralizedPowerAllocationVLSFP(L, K, gainOverNoisedB, D, rho_tot, b_vlsfp, upsilon, kappa, theta)

%% Centralized fractional power allocation
pow = zeros(K,1);
maxPow = zeros(K,1);
rho_vlsfp = zeros(L,K);
b_k = zeros(L,K);

for k = 1:K
    %Extract which UEs are served by AP l
    servingAPs = find(D(:,k)==1);
    pow(k) =  sum(gainOverNoisedB(servingAPs,k).^theta).^upsilon;
    pow(k) = pow(k)/max(abs(b_vlsfp(servingAPs,k)).^2)^kappa;
    maxPow(k) = max(abs(b_vlsfp(servingAPs,k)).^2);
end
normalizationFactor = 0;

for ell = 1:L
    servedUEs = find(D(ell,:)==1);
    temporScalar = maxPow(servedUEs)'*pow(servedUEs)/rho_tot;
    normalizationFactor = max(normalizationFactor,temporScalar);
end

rho_cent = pow/normalizationFactor;

for k = 1:K

    for l = 1:L
        rho_vlsfp(l,k) = rho_cent(k)*abs(b_vlsfp(l,k))^2;
    end

    b_k(:,k) = sqrt(rho_cent(k))*b_vlsfp(:,k);
end

end