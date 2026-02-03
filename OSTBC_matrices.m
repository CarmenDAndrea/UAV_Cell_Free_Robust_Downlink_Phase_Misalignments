% This function computes the orthogonal space-time block coding (OSTBC) matrices
% for a MIMO system. OSTBC matrices are deterministic matrices used to perform
% space-time block transmissions, aiming to improve diversity by transmitting
% symbols redundantly over time and space.
% Author: Marx Freitas

% INPUT
% N_period = Number of symbol periods per code matrix (bloc)
% s  = DIM: 1 x Nt - symbols to be transmitted at time block "t" (denoted
%       as block "k" in some references)

% OUTPUT
% X   = DIM: Nt x N_period - the "code matrix" to be transmitted.

% Observations
% "N_period" represents the number of times that a symbol "s" was
% transmitted over the time, i.e., it is transmitted in N epochs
% "N_period" is also called "Number of time intervals"

% Please, note that "N" does not represent the time of a symbol, but the
% time of a block (code matrix), which contains many symbols

function X = OSTBC_matrices (N_period, s)

switch N_period


    case 2
        % Example 7.2: The Alamouti Code is an OSTBC
        % Consider the Alamouti code in Section 6.3.1:
        % Here N = ns = 2 and hence the code rate is equal to R = 1
        X = [s(1),  conj(s(2));...
            s(2), -conj(s(1))];

    case 4
        % Example 7.4: OSTBC for nt = 4.
        % For nt = 4, N = 4, ns = 3 the following code is an orthogonal STBC:
        % This code has also a rate of R = 3/4
        X = [s(1), 0, s(2), -s(3);...
            0, s(1), conj(s(3)), conj(s(2));...
            -conj(s(2)), -s(3), conj(s(1)), 0;...
            conj(s(3)), -s(2), 0, conj(s(1))];

    case 8
        % Example 7.5: OSTBC for nt = 8.
        % For nt = 8, N = 8, ns = 4 the following code is an orthogonal
        % STBC (with rate R = 1/2):

        X = [s(1), 0, 0, 0, s(4), 0, s(2), -s(3);...
            0, s(1), 0, 0, 0, s(4), conj(s(3)), conj(s(2));...
            0, 0, s(1), 0, -conj(s(2)), -s(3), conj(s(4)), 0;...
            0, 0, 0, s(1), conj(s(3)), -s(2), 0, conj(s(4));...
            -conj(s(4)), 0, s(2), -s(3), conj(s(1)), 0, 0, 0;...
            0, -conj(s(4)), conj(s(3)), conj(s(2)), 0, conj(s(1)), 0, 0;...
            -conj(s(2)), -s(3), -s(4), 0, 0, 0, conj(s(1)), 0;...
            conj(s(3)), -s(2), 0, -s(4), 0, 0, 0, conj(s(1))];

end