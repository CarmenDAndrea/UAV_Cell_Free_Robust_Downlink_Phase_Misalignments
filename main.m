clear; clc; tic
%=========================================================================%
%% Input parametes
%=========================================================================%
% Choose the access point (AP) selection strategy utilized to associate the
% user-equipment (UE) with a subset of APs.
APselectionMethod = 'StrongestAPs'; % Options: 'CanonicalCF', 'UCC', 'LSFB', 'ScalableCF', 'matchedDecision', 'MD_LSFB', 'APselectionBook', 'StrongestAPs'
phaseMisalignments = 'Yes'; % 'Yes', 'No';
alpha_pi = pi; % Degree of misalignment % pi/4
nbrOfUEs = 20; % Number of UEs
nbrOfAPs = 40; % Number of APs
N = 4; % Number of antenna elements per AP (N)
U_max = 10; % Maximum number of UEs per AP
L_k_DSTBC = 2; % Number of APs per user
M = 8; % Set the modulation order
nbrOfSetups = 400; % number of monte-Carlo setups
nbrOfRealizations = 400; % number of channel realizations

% Set the number of samples per coherence block and for pilot signaling
tau_c = 200; % Number of samples per coherence block
tau_p = 10; % Number of sample for pilot signaling

% Transmission power of uplink (UL) and downlink (DL)
poweUL_mW = 100;
PowerDL_mW = 200;

% Propagation model parameters
BW_Hz = 20e6; % Bandwidth in Hz
fc_GHz = 3.5; % Frequency center in GHz
noiseFigure_dB = 8; % Noise figure in dB

%=========================================================================%
% Processing
%=========================================================================%

% Maximum number of UEs per AP
U_max = min(U_max, tau_p);

% Bits
modulationType = 'psk'; % Options: 'psk', 'qam'
bitMap = generateBitMap(M, modulationType);

% Generate constellation and amicable matrices
GenerateConstellation

BER_tot = cell(nbrOfSetups, 1);

% Go through all setups
for nSetup = 1:nbrOfSetups

    % Display simulation progress
    disp(['Setup ' num2str(nSetup) ' out of ' num2str(nbrOfSetups)]);

    % Places the APs and UEs in the coverage area, generates the channel
    % and compute the covariance matrix
    MF_APs_and_UEs_Positioning_and_ChannelandCovarianceEstimation

    % Performs AP clustering and channel estimation
    MF_APselection_and_ChannelEstimation

    % Power allocation distributed processing
    rho_dist_mW = DLpowerAllocation(nbrOfUEs, nbrOfAPs, D, PowerDL_mW, gainOverNoise_dB);
    
    % Generating the precoding vectors
    [w_MR, w_LP_MMSE, w_P_MMSE, w_P_RZF, ~, ~] = functionComputeSE_downlinkWithRician_LIGHT(Hhat, HnoPM, D, C,...
    tau_c, tau_p, nbrOfRealizations, N, nbrOfUEs, nbrOfAPs, poweUL_mW, rho_dist_mW, gainOverNoise_dB, PowerDL_mW);

    % Data the Detection, SER and BER
    Differential_and_Traditional_CellFree

    BER_tot{nSetup} = BER;

    clear SER BER

end

toc

