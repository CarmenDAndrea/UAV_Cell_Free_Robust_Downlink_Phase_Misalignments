% Set the stream pattern per UE
streamPerUE_CF_DPSK = 'Same';

% Number of data symbols
tau_d = tau_c-tau_p;

% ----------------------------------------------------------------------- %
% Generating the transmitted symbols
% ----------------------------------------------------------------------- %

% ------------------------- %
% Same stream per UE
% ------------------------- %
if strcmp(streamPerUE_CF_DPSK, 'Same')
    indicesTx_DPSK(:, 1, :)   = ceil (M*rand(nbrOfUEs, 1, tau_d));
    indicesTx_CF_PSK(:, 1, :) = ceil (M*rand(nbrOfUEs, 1, tau_d));

    for indexStreamUE = 2:N_streamsUE
        indicesTx_DPSK(:, indexStreamUE, :)   = indicesTx_DPSK(:, 1, :);
        indicesTx_CF_PSK(:, indexStreamUE, :) = indicesTx_CF_PSK(:, 1, :);
    end
else
    % ------------------------- %
    % Different streams per UE
    % ------------------------- %
    indicesTx_DPSK   = ceil (M*rand(nbrOfUEs, N_streamsUE, tau_d));
    indicesTx_CF_PSK = ceil (M*rand(nbrOfUEs, N_streamsUE, tau_d));

    % Different streams per UE
    transmittedSymbols_DPSK   = constellation(indicesTx_DPSK);
    transmittedSymbols_CF_PSK = constellation(indicesTx_CF_PSK);

    % Reference symbol for differential PSK
    transmittedSymbols_DPSK(:, :, 1)   = 1;
end

% Different streams per UE
transmittedSymbols_DPSK   = constellation(indicesTx_DPSK);
transmittedSymbols_CF_PSK = constellation(indicesTx_CF_PSK);

% Reference symbol for differential PSK
transmittedSymbols_DPSK(:, :, 1)   = 1;

clear indicesTx_DPSK indicesTx_CF_PSK cil