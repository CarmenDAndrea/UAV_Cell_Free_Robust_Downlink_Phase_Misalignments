% ----------------------------------------------------------------------- %
% Generate the constellation and symbols
% ----------------------------------------------------------------------- %
% Generating QPSK symbols
constellation = MF_QPSK_unitaryConstellation(M);

% Number of code matrices (blocks)
nbrOfBlocks = floor((tau_c-tau_p)/L_k_DSTBC);

% Configuring the symbol period
N_period = configureSymbolsPeriods_N(L_k_DSTBC);

% Number of symbols in each code matrix (block)
n_s = numberOfsymbolsOSTBC(N_period);

% Number of useful symbols symbols to be transmitted
nbrOfSymbols_DSTBC = n_s*(nbrOfBlocks-1);

% Using Equation 9.6.30
[An, Bn] = amicable_OSTBC_matrices(N_period, n_s);