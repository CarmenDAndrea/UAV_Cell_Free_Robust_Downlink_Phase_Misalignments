function qpskConst= MF_QPSK_unitaryConstellation(M)


% % Geração da constelação QPSK no intervalo [M-1, M+1]
% scaleFactor = (M+1) - (M-1); % Amplitude do intervalo
% offset = M-1; % Deslocamento inicial
% qpskConst = offset + scaleFactor * exp(1i * (0:M-1) * 2 * pi / M);
% 
% % Ajuste para média zero
% meanValue = mean(qpskConst); % Calcula a média
% qpskConst = qpskConst - meanValue; % Subtrai a média para centralizar
% 
% % Exibição da constelação
% disp('Constelação QPSK ajustada para média zero:');
% disp(qpskConst);
% 
% % Verificar a nova média
% disp(['Média ajustada da constelação: ', num2str(mean(qpskConst))]);



% Geração da constelação QPSK no círculo unitário
qpskConst = exp(1i * (0:M-1) * 2 * pi / M);


% % Exibição da constelação no plano complexo
% disp('Constelação QPSK:');
% disp(qpskConst);
% 
% % Visualização no plano complexo
% figure;
% plot(real(qpskConst), imag(qpskConst), 'o', 'MarkerSize', 10, 'LineWidth', 2);
% grid on;
% axis equal;
% xlabel('Parte Real');
% ylabel('Parte Imaginária');
% title('Constelação QPSK');
