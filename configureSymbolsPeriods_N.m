% This function sets the number of symbol periods used by each orthogonal
% code matrix. The symbol period refers to the number of transmission time
% that symbols are transmitted.

% Example: Nt x N_period (3 transmitting antennas and 4 symbol periods
%      [s1, s2, s3,    0,
% X =  -s2, s3,  0, -s1*,
%       0 , s1, s3, -s2,];

% INPUT
% Nt = The number of transmitting antennas

% OUTPUT
% N_period = The symbol period

function N_period = configureSymbolsPeriods_N(Nt)

switch Nt

    case 2
        N_period = 2;

    case 3
        N_period = 4;

    case 4
        N_period = 4;

    case 5
        N_period = 8;

    case 6
        N_period = 8;

    case 7
        N_period = 8;

    case 8
        N_period = 8;

end