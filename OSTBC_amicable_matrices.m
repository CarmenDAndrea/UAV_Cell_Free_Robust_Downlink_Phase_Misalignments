% Author: Marx Freitas

function X = OSTBC_amicable_matrices (An, Bn, s, n_s)

X = 0;

for indexSymbol = 1: n_s

    s_n_k_bar = (s(indexSymbol) + conj(s(indexSymbol)))/2;
    s_n_k_til = (s(indexSymbol) - conj(s(indexSymbol)))/2;
    X = X + (s_n_k_bar*An(:,:, indexSymbol) + s_n_k_til*Bn(:, :, indexSymbol));

end

end