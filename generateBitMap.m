function bitMap = generateBitMap(M, modulationType)
    % Generate bit mapping using Gray coding for QAM and PSK
    % Input:
    %   M - Modulation order (must be a power of 2)
    %   modulationType - 'qam' or 'psk'
    % Output:
    %   bitMap - Bit mapping for the constellation indices

    if log2(M) ~= round(log2(M))
        error('M must be a power of 2.');
    end

    k = log2(M); % Bits per symbol

    % Generate symbol indices (0 to M-1)
    symbols = (0:M-1)';

    % Apply Gray-coded modulation
    if strcmp(modulationType, 'qam')
        modulatedSymbols = qammod(symbols, M, 'gray', 'UnitAveragePower', true);
    elseif strcmp(modulationType, 'psk')
        modulatedSymbols = pskmod(symbols, M, pi/M, 'gray'); % pi/M phase shift for proper Gray coding
    else
        error('Unsupported modulation type. Use "qam" or "psk".');
    end

    % Obtain Gray-coded indices by sorting the modulated symbols
    [~, grayIndices] = sort(modulatedSymbols); 

    % Ensure the Gray indices are correctly mapped from 0 to M-1
    grayIndices = uint32(grayIndices - 1); % Convert to zero-based indexing

    % Convert Gray indices to binary bit mapping
    bitMap = de2bi(grayIndices, k, 'left-msb');
end
