function pamConst = MF_pamSquareConstellation(M)
% function qamConst=ak_qamSquareConstellation(M)
%Generates a QAM square constellation with M symbols.

if M<4
    error('Number M of symbols must be >= 4');
end

pamConst = -(M-1):2:M-1; %use PAM constellation

% % MF modification: Unity-power symbols
% % Calculando a potência média dos símbolos
% P_media = mean(abs(qamConst).^2);
% 
% % Normalizando para garantir que a potência média seja 1
% qamConst = qamConst / sqrt(P_media);