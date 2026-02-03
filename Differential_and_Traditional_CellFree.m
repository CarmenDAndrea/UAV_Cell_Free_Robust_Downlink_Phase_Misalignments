% This code computes the bit-error-rate (BER) and spectral efficiency (SE)
% of conventional cell-free massive MIMO systems and differential schemes,
% i.e., DSTBC and DPSK

% ----------------------------------------------------------------------- %
% Prepare to store the received signal
% ----------------------------------------------------------------------- %
noiseVariance_mW = db2pow(noiseVariance_dBm);

% Prepare to store the bit error rate (BER)
BER_DSTBC_MR      = zeros(nbrOfUEs, nbrOfRealizations);
BER_DSTBC_LP_MMSE = zeros(nbrOfUEs, nbrOfRealizations);
BER_DSTBC_P_MMSE  = zeros(nbrOfUEs, nbrOfRealizations);
BER_DSTBC_P_RZF   = zeros(nbrOfUEs, nbrOfRealizations);

BER_DPSK_MR      = zeros(nbrOfUEs, nbrOfRealizations);
BER_DPSK_LP_MMSE = zeros(nbrOfUEs, nbrOfRealizations);
BER_DPSK_P_MMSE  = zeros(nbrOfUEs, nbrOfRealizations);
BER_DPSK_P_RZF   = zeros(nbrOfUEs, nbrOfRealizations);

BER_CF_PSK_MR      = zeros(nbrOfUEs, nbrOfRealizations);
BER_CF_PSK_LP_MMSE = zeros(nbrOfUEs, nbrOfRealizations);
BER_CF_PSK_P_MMSE  = zeros(nbrOfUEs, nbrOfRealizations);
BER_CF_PSK_P_RZF   = zeros(nbrOfUEs, nbrOfRealizations);

BER_CF_PSK_MR_noPM      = zeros(nbrOfUEs, nbrOfRealizations);
BER_CF_PSK_LP_MMSE_noPM = zeros(nbrOfUEs, nbrOfRealizations);
BER_CF_PSK_P_MMSE_noPM  = zeros(nbrOfUEs, nbrOfRealizations);
BER_CF_PSK_P_RZF_noPM   = zeros(nbrOfUEs, nbrOfRealizations);

% ----------------------------------------------------------------------- %
% PROCESSING
% ----------------------------------------------------------------------- %
for indexRealization = 1:nbrOfRealizations

    % ------------------------------------------------------------------- %
    % DIFFERENTIAL OSTBC
    % ------------------------------------------------------------------- %
    indicesTx_DSTBC               = ceil (M*rand(nbrOfUEs, nbrOfSymbols_DSTBC));
    transmittedSymbols_DSTBC      = constellation(indicesTx_DSTBC);

    % ------------------------------------------------------------------- %
    % Settings for code matrices
    % ------------------------------------------------------------------- %
    % Compute the matrix W (relates the previous and current code matrices
    C_0 = eye(N_period, N_period); % Set C_0 to an identity
    C_t = zeros(nbrOfUEs, N_period, N_period, nbrOfBlocks);
    cil = zeros(nbrOfUEs, nbrOfAPs, nbrOfBlocks, N_period);

    % ------------------------------------------------------------------- %
    % THE CPU GENERATES THE CODE MATRICES (BLOCKS) FOR EACH UE
    % ------------------------------------------------------------------- %
    for indexUE =1:nbrOfUEs

        % Prepare to store the code matrix of UE k
        C_t_aux = zeros(N_period, N_period, nbrOfBlocks);

        for t = 1:nbrOfBlocks

            if t == 1
                C_t_aux(:,:, t)  = C_0;
                C_t(indexUE, :, :, t) = C_t_aux(:,:, t);
            else

                % Set the symbols to be transmitted in each code matrix (blocl)
                firstSymbol_id = (t-2)*n_s + 1;
                lastSymbol_id  = (t-1)*n_s;
                symbols_per_block = transmittedSymbols_DSTBC(indexUE, firstSymbol_id:lastSymbol_id);

                % Generates an orthogonal code matrix to transmit the symbols
                X_block_t = OSTBC_matrices(N_period, symbols_per_block);
                X_block_t = X_block_t/sqrt(n_s);

                % Relashionship between code matrices (block) in two time instants
                C_t_aux(:,:, t)  = C_t_aux(:,:, t-1)*X_block_t;
                C_t(indexUE, :, :, t) = C_t_aux(:,:, t);
            end
        end
    end

    % ------------------------------------------------------------------- %
    % THE CPU SPLITS C_k AMONG THE APs SERVING THE UE
    % ------------------------------------------------------------------- %
    for t=1:nbrOfBlocks

        for indexUE = 1:nbrOfUEs
            servingAPs = find(D(:, indexUE) == 1);

            % Compute the number of APs serving the UE
            L_k = length(servingAPs);

            if (L_k < L_k_DSTBC) || (L_k > L_k_DSTBC)
                error('L_k  as to be equal to L_k_DSTBC');
            end

            for index = 1:length(servingAPs)
                cil(indexUE, servingAPs(index), t, :)  = C_t(indexUE, index, :, t);
            end
        end
    end

    % ----------------------------------------------------------------------- %
    % DATA DETECTION: DIFFERENTIAL OSTBC
    % ----------------------------------------------------------------------- %

    receivedSymbols_DSTBC_MR = receivedSymbols_DSTBC(indexRealization, nbrOfUEs, nbrOfBlocks, N_period,...
        nbrOfAPs, N, nbrOfSymbols_DSTBC, n_s, HnoPM, noiseVariance_dBm, D, w_MR, cil, constellation, An, Bn, alpha_pi);

    receivedSymbols_DSTBC_LP_MMSE = receivedSymbols_DSTBC(indexRealization, nbrOfUEs, nbrOfBlocks, N_period,...
        nbrOfAPs, N, nbrOfSymbols_DSTBC, n_s, HnoPM, noiseVariance_dBm, D, w_LP_MMSE, cil, constellation, An, Bn, alpha_pi);

    receivedSymbols_DSTBC_P_MMSE = receivedSymbols_DSTBC(indexRealization, nbrOfUEs, nbrOfBlocks, N_period,...
        nbrOfAPs, N, nbrOfSymbols_DSTBC, n_s, HnoPM, noiseVariance_dBm, D, w_P_MMSE, cil, constellation, An, Bn, alpha_pi);

    receivedSymbols_DSTBC_P_RZF = receivedSymbols_DSTBC(indexRealization, nbrOfUEs, nbrOfBlocks, N_period,...
        nbrOfAPs, N, nbrOfSymbols_DSTBC, n_s, HnoPM, noiseVariance_dBm, D, w_P_RZF, cil, constellation, An, Bn, alpha_pi);

    % ------------------------------------------------------------------- %
    % DATA DETECTION: DIFFERENTIAL PSK AND TRADITIONAL CELL-FREE
    % ------------------------------------------------------------------- %
    tau_d = tau_c-tau_p; % Number of data symbols

    % Transmitted symbols
    indicesTx_DPSK   = ceil (M*rand(nbrOfUEs, tau_d));
    indicesTx_CF_PSK = ceil (M*rand(nbrOfUEs, tau_d));

    transmittedSymbols_DPSK   = constellation(indicesTx_DPSK);
    transmittedSymbols_CF_PSK = constellation(indicesTx_CF_PSK);

    % Reference symbol for differential PSK
    transmittedSymbols_DPSK(:, 1)   = 1;

    [receivedSymbols_CF_PSK_MR, receivedSymbols_CF_PSK_MR_noPM, receivedSymbols_DPSK_MR] = receivedSymbolsCF_DPSK_opt(indexRealization, nbrOfUEs, ...
        nbrOfAPs, N, tau_d, HnoPM, noiseVariance_dBm, transmittedSymbols_DPSK, transmittedSymbols_CF_PSK, constellation, D, w_MR, alpha_pi);

    [receivedSymbols_CF_PSK_LP_MMSE, receivedSymbols_CF_PSK_LP_MMSE_noPM, receivedSymbols_DPSK_LP_MMSE] = receivedSymbolsCF_DPSK_opt(indexRealization, nbrOfUEs, ...
        nbrOfAPs, N, tau_d, HnoPM, noiseVariance_dBm, transmittedSymbols_DPSK, transmittedSymbols_CF_PSK, constellation, D, w_LP_MMSE, alpha_pi);

    [receivedSymbols_CF_PSK_P_RZF, receivedSymbols_CF_PSK_P_RZF_noPM, receivedSymbols_DPSK_P_RZF] = receivedSymbolsCF_DPSK_opt(indexRealization, nbrOfUEs, ...
        nbrOfAPs, N, tau_d, HnoPM, noiseVariance_dBm, transmittedSymbols_DPSK, transmittedSymbols_CF_PSK, constellation, D, w_P_RZF, alpha_pi);

    [receivedSymbols_CF_PSK_P_MMSE, receivedSymbols_CF_PSK_P_MMSE_noPM, receivedSymbols_DPSK_P_MMSE] = receivedSymbolsCF_DPSK_opt(indexRealization, nbrOfUEs, ...
        nbrOfAPs, N, tau_d, HnoPM, noiseVariance_dBm, transmittedSymbols_DPSK, transmittedSymbols_CF_PSK, constellation, D, w_P_MMSE, alpha_pi);

    % ------------------------------------------------------------------- %
    % COMPUTE BER, SER
    % ------------------------------------------------------------------- %
    % DSTBC
    BER_DSTBC_MR(:, indexRealization) =...
        SER_BER(nbrOfUEs, transmittedSymbols_DSTBC, receivedSymbols_DSTBC_MR, constellation, bitMap);

    BER_DSTBC_LP_MMSE(:, indexRealization) =...
        SER_BER(nbrOfUEs, transmittedSymbols_DSTBC, receivedSymbols_DSTBC_LP_MMSE, constellation, bitMap);

    BER_DSTBC_P_MMSE(:, indexRealization) =...
        SER_BER(nbrOfUEs, transmittedSymbols_DSTBC, receivedSymbols_DSTBC_P_MMSE, constellation, bitMap);

    BER_DSTBC_P_RZF(:, indexRealization) =...
        SER_BER(nbrOfUEs, transmittedSymbols_DSTBC, receivedSymbols_DSTBC_P_RZF, constellation, bitMap);

    % DPSK
    BER_DPSK_MR(:, indexRealization) =...
        SER_BER(nbrOfUEs, transmittedSymbols_DPSK(:,2:end), receivedSymbols_DPSK_MR(:,2:end), constellation, bitMap);

    BER_DPSK_LP_MMSE(:, indexRealization) =...
        SER_BER(nbrOfUEs, transmittedSymbols_DPSK(:,2:end), receivedSymbols_DPSK_LP_MMSE(:,2:end), constellation, bitMap);

    BER_DPSK_P_MMSE(:, indexRealization) =...
        SER_BER(nbrOfUEs, transmittedSymbols_DPSK(:,2:end), receivedSymbols_DPSK_P_MMSE(:,2:end), constellation, bitMap);

    BER_DPSK_P_RZF(:, indexRealization) =...
        SER_BER(nbrOfUEs, transmittedSymbols_DPSK(:,2:end), receivedSymbols_DPSK_P_RZF(:,2:end), constellation, bitMap);


    % Conventional Cell-Free
    BER_CF_PSK_MR(:, indexRealization) =...
        SER_BER(nbrOfUEs, transmittedSymbols_CF_PSK, receivedSymbols_CF_PSK_MR, constellation, bitMap);

    BER_CF_PSK_LP_MMSE(:, indexRealization) =...
        SER_BER(nbrOfUEs, transmittedSymbols_CF_PSK, receivedSymbols_CF_PSK_LP_MMSE, constellation, bitMap);

    BER_CF_PSK_P_MMSE(:, indexRealization) =...
        SER_BER(nbrOfUEs, transmittedSymbols_CF_PSK, receivedSymbols_CF_PSK_P_MMSE, constellation, bitMap);

    BER_CF_PSK_P_RZF(:, indexRealization) =...
        SER_BER(nbrOfUEs, transmittedSymbols_CF_PSK, receivedSymbols_CF_PSK_P_RZF, constellation, bitMap);


    % Conventional Cell-Free with no phase misalignment
    BER_CF_PSK_MR_noPM(:, indexRealization) =...
        SER_BER(nbrOfUEs, transmittedSymbols_CF_PSK, receivedSymbols_CF_PSK_MR_noPM, constellation, bitMap);

    BER_CF_PSK_LP_MMSE_noPM(:, indexRealization) =...
        SER_BER(nbrOfUEs, transmittedSymbols_CF_PSK, receivedSymbols_CF_PSK_LP_MMSE_noPM, constellation, bitMap);

    BER_CF_PSK_P_MMSE_noPM(:, indexRealization) =...
        SER_BER(nbrOfUEs, transmittedSymbols_CF_PSK, receivedSymbols_CF_PSK_P_MMSE_noPM, constellation, bitMap);

    BER_CF_PSK_P_RZF_noPM(:, indexRealization) =...
        SER_BER(nbrOfUEs, transmittedSymbols_CF_PSK, receivedSymbols_CF_PSK_P_RZF_noPM, constellation, bitMap);

end

% ----------------------------------------------------------------------- %
% BER
% ----------------------------------------------------------------------- %
BER = struct;

BER.BER_DSTBC_MR  = BER_DSTBC_MR;
BER.BER_DPSK_MR   = BER_DPSK_MR;
BER.BER_CF_PSK_MR = BER_CF_PSK_MR;
BER.BER_CF_PSK_MR_noPM = BER_CF_PSK_MR_noPM;

% BER LP-MMSE
BER.BER_DSTBC_LP_MMSE  = BER_DSTBC_LP_MMSE;
BER.BER_DPSK_LP_MMSE   = BER_DPSK_LP_MMSE;
BER.BER_CF_PSK_LP_MMSE = BER_CF_PSK_LP_MMSE;
BER.BER_CF_PSK_LP_MMSE_noPM = BER_CF_PSK_LP_MMSE_noPM;

% BER P-MMSE
BER.BER_DSTBC_P_MMSE  = BER_DSTBC_P_MMSE;
BER.BER_DPSK_P_MMSE   = BER_DPSK_P_MMSE;
BER.BER_CF_PSK_P_MMSE = BER_CF_PSK_P_MMSE;
BER.BER_CF_PSK_P_MMSE_noPM = BER_CF_PSK_P_MMSE_noPM;

% BER P-RZF
BER.BER_DSTBC_P_RZF  = BER_DSTBC_P_RZF;
BER.BER_DPSK_P_RZF   = BER_DPSK_P_RZF;
BER.BER_CF_PSK_P_RZF = BER_CF_PSK_P_RZF;
BER.BER_CF_PSK_P_RZF_noPM = BER_CF_PSK_P_RZF_noPM;

function BER = computeBER(transmittedSymbols, receivedSymbols, constellation, bitMap)

[~, indicesTx] = ismember(transmittedSymbols, constellation);
transmittedbits = bitMap(indicesTx, :);
transmittedbits = transmittedbits(:);

[~, indicesRx] = ismember(receivedSymbols, constellation);
receivedbits = bitMap(indicesRx, :);
receivedbits = receivedbits(:);

BER = sum(receivedbits ~= transmittedbits)/length(transmittedbits);

end

function receivedSymbols = receivedSymbols_DSTBC(indexRealization, nbrOfUEs, nbrOfBlocks, N_period, nbrOfAPs, ...
    N, nbrOfSymbols, n_s, HnoPM, noiseVariance_dBm, D, w_precoding, cil, constellation, An, Bn, alpha_pi)

% Prepare to store the signal transmitted by AP l
s_l_DSTBC = zeros(nbrOfAPs, N, N_period*nbrOfBlocks); % OK

% Prepare to store the received signal
y_i_t_DSTBC = zeros(nbrOfUEs, N_period*nbrOfBlocks);

% Converting the noise variance to mW
noiseVariance_mW = db2pow(noiseVariance_dBm);

% Prepare to store the phase misalignments
exp_phi_DL = zeros(1, nbrOfAPs);

% Signal transmitted by the AP l
for indexAP = 1:nbrOfAPs

    % Accessing the channel between the UE l and AP l
    startIdx_h = (indexAP-1)*N  + 1; endIdx_h   = startIdx_h + N - 1; % OK

    % Finding the UEs served by AP l
    servedUEs = find(D(indexAP, :) == 1); % OK

    if isempty(servedUEs)
        % Do nothing
    else
        % Compute the signal transmitted by AP l
        for id = 1:length(servedUEs) % OK
            % Precodings
            wil = w_precoding(startIdx_h:endIdx_h, servedUEs(id), indexRealization);
            cil_aux = cil(servedUEs(id), indexAP, :, :); % Extração dos elementos desejados
            cil_aux = reshape(cil_aux, [], N_period);   % Converte para (G × N_p)
            cil_aux = reshape(cil_aux.', 1, []);        % Transpõe e converte para (1 × (G × N_p))

            % Signal transmited by AP l
            s_l_DSTBC(indexAP, :, :) = s_l_DSTBC(indexAP, :, :) + reshape(wil * cil_aux, [1, N, N_period*nbrOfBlocks]);
        end
    end
    % Generating phase misalignments
    phi_shift_DL_aux = -alpha_pi + 2*alpha_pi*randn;
    exp_phi_DL_aux = exp(1i*phi_shift_DL_aux);
    exp_phi_DL(1, indexAP) = diag(exp_phi_DL_aux);
end

% Received signal
for indexUE = 1:nbrOfUEs
    % Prepare to store the received signal at UE i, at block t
    y_i_DSTBC = 0;

    for indexAP = 1:nbrOfAPs
        % Accessing the channel between the UE l and AP l
        startIdx_h = (indexAP-1)*N  + 1; endIdx_h   = startIdx_h + N - 1; % OK
        hli_hermitian = exp_phi_DL(1, indexAP)*transpose(conj(HnoPM(startIdx_h:endIdx_h, indexRealization, indexUE))); % OK

        % Received signal at UE i due to AP l
        y_i_DSTBC = y_i_DSTBC + hli_hermitian*reshape(s_l_DSTBC(indexAP, :, :), [N, N_period*nbrOfBlocks]); % OK
    end

    % Computing the noise
    noise = sqrt(0.5)*sqrt(noiseVariance_mW)*(randn(1, N_period*nbrOfBlocks) + 1i*randn(1, N_period*nbrOfBlocks));
    noise = noise./sqrt(noiseVariance_mW);
    y_i_t_DSTBC(indexUE, :) = y_i_DSTBC + noise; % OK
end

% Prepare to store the received symbols
receivedSymbols = zeros(nbrOfUEs, nbrOfSymbols);

% ----------------------------------------------------------------------- %
% ML DETECTION: DIFFERENTIAL OSTBC
% ----------------------------------------------------------------------- %

for indexUE = 1:nbrOfUEs
    Y  = zeros(1, N_period, nbrOfBlocks);

    for t = 1:nbrOfBlocks
        startIdx_t = (t-1)*N_period  + 1;
        endIdx_t   = startIdx_t + N_period - 1; % OK
        y_k_t_aux = y_i_t_DSTBC(indexUE, startIdx_t:endIdx_t);
        Y(1, :, t) = y_k_t_aux;

        if t == 1
            % Do nothing
        else

            % Firs step of differential detection
            argument_DSTBC = transpose(conj(Y(1, :, t)))*Y(1, :, t-1);

            % Detecting using amicable orthogonal designs
            firstSymbol_id = (t-2)*n_s + 1;
            lastSymbol_id  = (t-1)*n_s;
            symbols_id_t   = firstSymbol_id:lastSymbol_id;

            for indexSymbol_aux = 1:n_s
                real_computation = real(trace(An(:,:,indexSymbol_aux)*argument_DSTBC))*real(constellation);
                imag_computation = imag(trace(Bn(:,:,indexSymbol_aux)*argument_DSTBC))*imag(constellation);
                ML = -real_computation + imag_computation;

                % Computing the received symbols
                min_ML = min(ML(1, :));
                id = find(ML(1, :) == min_ML);
                receivedSymbols(indexUE, symbols_id_t(indexSymbol_aux)) = constellation(id);
            end
        end
    end
end
end


% ----------------------------------------------------------------------- %
% ADDITIONAL FUNCTIONS
% ----------------------------------------------------------------------- %

function [receivedSymbols_CF_PSK, receivedSymbols_CF_PSK_noPM, receivedSymbols_DPSK] = receivedSymbolsCF_DPSK_opt(indexRealization, nbrOfUEs,...
    nbrOfAPs, N, tau_d, HnoPM, noiseVariance_dBm, transmittedSymbols_DPSK, transmittedSymbols_CF_PSK, constellation, D, w_precoding, alpha_pi)
% ----------------------------------------------------------------------- %
% RECEIVED SIGNAL: DIFFERENTIAL PSK AND TRADITIONAL CELL-FREE
% ----------------------------------------------------------------------- %
% Prepare to store the received signals
y_i_t_DPSK         = zeros(nbrOfUEs, tau_d);
y_i_t_CF_PSK       = zeros(nbrOfUEs, tau_d);
y_i_t_CF_PSK_noPM  = zeros(nbrOfUEs, tau_d);
cil_aux_DPSK  = zeros(nbrOfUEs, tau_d);

% Prepare to store the received symbols
receivedSymbols_CF_PSK       = zeros(nbrOfUEs, tau_d);
receivedSymbols_CF_PSK_noPM  = zeros(nbrOfUEs, tau_d);
receivedSymbols_DPSK         = zeros(nbrOfUEs, tau_d);

% Generating the transmitted signal for DPSK
cil_aux_DPSK(:,1) = transmittedSymbols_DPSK(:,1); % Inicializa o primeiro símbolo corretamente
cil_aux_DPSK(:, 2:tau_d) = cumprod(transmittedSymbols_DPSK(:, 2:tau_d), 2);

% Prepare to store the signal transmitted by AP l
s_l_DPSK = zeros(nbrOfAPs, N, tau_d); % OK
s_l_CF_PSK = zeros(nbrOfAPs, N, tau_d); % OK

% Converte the noiseVariance to mW
noiseVariance_mW = db2pow(noiseVariance_dBm);

% Prepare to store the phase misalignments
exp_phi_DL = zeros(1, nbrOfAPs);

% Signal transmitted by the AP l
for indexAP = 1:nbrOfAPs

    % Accessing the channel between the UE l and AP l
    startIdx_h = (indexAP-1)*N  + 1; endIdx_h   = startIdx_h + N - 1; % OK

    % Finding the UEs served by AP l
    servedUEs = find(D(indexAP, :) == 1); % OK

    if isempty(servedUEs)
        % Do nothing
    else
        % Compute the signal transmitted by AP l
        for id = 1:length(servedUEs) % OK
            % Precodings
            wil = w_precoding(startIdx_h:endIdx_h, servedUEs(id), indexRealization);

            % Signal transmited by the AP l
            s_l_DPSK(indexAP, :, :)   = s_l_DPSK(indexAP, :, :)   + reshape(wil * cil_aux_DPSK(servedUEs(id), :), [1, N, tau_d]);
            s_l_CF_PSK(indexAP, :, :) = s_l_CF_PSK(indexAP, :, :) + reshape(wil*transmittedSymbols_CF_PSK(servedUEs(id), :), [1, N, tau_d]); % OK
        end
    end
    % Generating phase misalignments
    phi_shift_DL_aux = -alpha_pi + 2*alpha_pi*randn;
    exp_phi_DL_aux = exp(1i*phi_shift_DL_aux);
    exp_phi_DL(1, indexAP) = diag(exp_phi_DL_aux);
end

% Received signal
for indexUE = 1:nbrOfUEs

    % Prepare to store the received signal at UE i, time instant t
    y_i_DPSK = 0; y_i_CF_PSK = 0; y_i_CF_PSK_noPM = 0;

    for indexAP = 1:nbrOfAPs
        % Accessing the channel between the UE l and AP l
        startIdx_h = (indexAP-1)*N  + 1; endIdx_h   = startIdx_h + N - 1; % OK
        hli_hermitian = exp_phi_DL(1, indexAP)*transpose(conj(HnoPM(startIdx_h:endIdx_h, indexRealization, indexUE))); % OK
        hli_hermitian_noPM = transpose(conj(HnoPM(startIdx_h:endIdx_h, indexRealization, indexUE))); % OK

        % Received signal at UE i due to AP l
        y_i_DPSK = y_i_DPSK + hli_hermitian*reshape(s_l_DPSK(indexAP, :, :), [N, tau_d]); % OK
        y_i_CF_PSK = y_i_CF_PSK + hli_hermitian*reshape(s_l_CF_PSK(indexAP, :, :), [N, tau_d]); % OK
        y_i_CF_PSK_noPM = y_i_CF_PSK_noPM + hli_hermitian_noPM*reshape(s_l_CF_PSK(indexAP, :, :), [N, tau_d]); % OK
    end

    % Computing the noise
    noise = sqrt(0.5)*sqrt(noiseVariance_mW)*(randn(1, tau_d) + 1i*randn(1, tau_d));
    noise = noise./sqrt(noiseVariance_mW);

    % Store the received signal at UE i, time instant t
    y_i_t_DPSK(indexUE, :)        = y_i_DPSK + noise; % OK
    y_i_t_CF_PSK(indexUE, :)      = y_i_CF_PSK + noise; % OK
    y_i_t_CF_PSK_noPM(indexUE, :) = y_i_CF_PSK_noPM + noise; % OK
end

% ----------------------------------------------------------------------- %
% ML DETECTION: DIFFERENTIAL PSK AND TRADITIONAL CELL-FREE
% ----------------------------------------------------------------------- %
for indexUE = 1:nbrOfUEs
    % ----------------------------------------- %
    % Computing
    % ----------------------------------------- %
    for t = 1:tau_d
        % Compute the ML detection for the traditional cell-free
        ML_aux_CF_PSK = abs(mod(angle(y_i_t_CF_PSK(indexUE, t)), 2*pi) - mod(angle(constellation), 2*pi)).^2; % OK
        [~, id_CF_PSK] = min(ML_aux_CF_PSK); % OK
        receivedSymbols_CF_PSK(indexUE, t) = constellation(1, id_CF_PSK); % OK

        ML_aux_CF_PSK_noPM = abs(mod(angle(y_i_t_CF_PSK_noPM(indexUE, t)), 2*pi) - mod(angle(constellation), 2*pi)).^2; % OK
        [~, id_CF_PSK_noPM] = min(ML_aux_CF_PSK_noPM); % OK
        receivedSymbols_CF_PSK_noPM(indexUE, t) = constellation(1, id_CF_PSK_noPM); % OK

        if t == 1
            receivedSymbols_DPSK(indexUE, t) = 1;
        else
            argument_DPSK = conj(y_i_t_DPSK(indexUE, t))*y_i_t_DPSK(indexUE, t-1); % OK
            ML_aux_DPSK   = real(argument_DPSK.*constellation);
            [~, id_DPSK]  = max(ML_aux_DPSK); % OK
            receivedSymbols_DPSK(indexUE, t) = constellation(1, id_DPSK); % OK
        end
    end
end

end

function BER = SER_BER(nbrOfUEs, transmittedSymbols, receivedSymbols, constellation, bitMap)

% Initialize outputs
BER = zeros(nbrOfUEs, 1);

% Loop through each user
for indexUE = 1:nbrOfUEs

    % Get transmitted and received symbols for this user
    txSymbols = transmittedSymbols(indexUE, :);
    rxSymbols = receivedSymbols(indexUE, :);

    % Compute Bit Error Rate (BER)
    BER(indexUE, 1) = computeBER(txSymbols, rxSymbols, constellation, bitMap);

end
end