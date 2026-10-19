%% =====================================================================
%  Transformada de Fourier - Funcoes Basicas: RECT  <->  SINC
%  Exercicio "hands on" - Sinais e Sistemas / PDS
%
%  Par de transformadas:
%       x(t) = A * rect(t/tau)   <-->   X(f) = A * tau * sinc(f*tau)
%  com sinc(u) = sin(pi*u)/(pi*u)   (sinc normalizada)
%
%  O que o script faz:
%   1) Varre varios valores de tau (largura da porta) e plota rect e |X(f)|
%   2) Compara a TF numerica (FFT) com a TF analitica
%   3) Mostra o 1o nulo em f = 1/tau e o produto tempo x banda constante
%   4) Varia a amplitude A (2o parametro) e o atraso t0 (3o parametro)
%   5) Abre um modo interativo com sliders para tau, A e t0
%
%  Nao precisa de toolbox: a funcao sinc e definida no proprio arquivo.
%  Autor: Prof. Newton Licciardi Jr.
% =====================================================================
clear; close all; clc;

%% ---------------- Parametros de simulacao ---------------------------
fs   = 1000;              % taxa de amostragem [amostras/s]
T    = 40;                % janela total de observacao [s]  (resolucao df = 1/T)
t    = -T/2 : 1/fs : T/2 - 1/fs;   % eixo do tempo
N    = numel(t);
f    = (-N/2 : N/2-1) * (fs/N);    % eixo da frequencia (centrado)

taus = [0.5 1 2 4];       % larguras da porta a comparar [s]
A    = 1;                 % amplitude

cores = lines(numel(taus));

%% ---------------- Figura 1: efeito da largura tau --------------------
figure('Name','Rect x Sinc - variando tau','Color','w','Position',[80 80 1100 700]);

subplot(2,1,1); hold on; grid on;
for k = 1:numel(taus)
    x = A * rectp(t / taus(k));
    plot(t, x, 'LineWidth', 1.8, 'Color', cores(k,:), ...
        'DisplayName', sprintf('\\tau = %.1f s', taus(k)));
end
xlim([-3 3]); ylim([-0.1 1.3]);
xlabel('t [s]'); ylabel('x(t)');
title('Dominio do tempo: x(t) = A\cdot rect(t/\tau)');
legend('Location','northeast');

subplot(2,1,2); hold on; grid on;
for k = 1:numel(taus)
    Xa = A * taus(k) * sincn(f * taus(k));          % TF analitica
    plot(f, Xa, 'LineWidth', 1.8, 'Color', cores(k,:), ...
        'DisplayName', sprintf('\\tau = %.1f s  (1o nulo em %.2f Hz)', taus(k), 1/taus(k)));
end
xlim([-5 5]);
xlabel('f [Hz]'); ylabel('X(f)');
title('Dominio da frequencia: X(f) = A\cdot\tau\cdot sinc(f\tau)');
legend('Location','northeast');

%% ---------------- Figura 2: FFT numerica x TF analitica --------------
tau_teste = 1;
x   = A * rectp(t / tau_teste);
Xn  = fftshift(fft(ifftshift(x))) / fs;      % aproxima a integral de Fourier
Xa  = A * tau_teste * sincn(f * tau_teste);

figure('Name','FFT x analitica','Color','w','Position',[120 120 1000 420]);
plot(f, real(Xn), 'b', 'LineWidth', 2, 'DisplayName','FFT (numerica)'); hold on; grid on;
plot(f, Xa, 'r--', 'LineWidth', 1.5, 'DisplayName','\tau\cdot sinc(f\tau) (analitica)');
xlim([-6 6]);
xlabel('f [Hz]'); ylabel('X(f)');
title(sprintf('Validacao numerica, \\tau = %.1f s  |  erro maximo = %.2e', ...
      tau_teste, max(abs(real(Xn) - Xa))));
legend;

%% ---------------- Figura 3: espectro em dB (lobulos laterais) --------
figure('Name','Espectro em dB','Color','w','Position',[160 160 1000 420]);
hold on; grid on;
for k = 1:numel(taus)
    Xa = abs(taus(k) * sincn(f * taus(k)));
    plot(f, 20*log10(Xa / taus(k) + eps), 'LineWidth', 1.5, 'Color', cores(k,:), ...
        'DisplayName', sprintf('\\tau = %.1f s', taus(k)));
end
yline(-13.26, 'k--', '1o lobulo lateral: -13,3 dB', 'LabelHorizontalAlignment','left');
xlim([-5 5]); ylim([-60 3]);
xlabel('f [Hz]'); ylabel('|X(f)|/X(0)  [dB]');
title('Espectro normalizado em dB: a forma nao muda, so a escala');
legend('Location','northeast');

%% ---------------- Figura 4: relacao inversa tempo x banda ------------
tau_v  = linspace(0.1, 5, 200);
banda  = 1 ./ tau_v;                 % largura do lobulo principal (meia) = 1/tau
figure('Name','Tempo x Banda','Color','w','Position',[200 200 1000 420]);
subplot(1,2,1);
plot(tau_v, banda, 'LineWidth', 2); grid on;
xlabel('\tau [s]'); ylabel('1o nulo  f_0 = 1/\tau  [Hz]');
title('Porta mais larga  =>  espectro mais estreito');
subplot(1,2,2);
plot(tau_v, tau_v .* banda, 'LineWidth', 2); grid on; ylim([0 2]);
xlabel('\tau [s]'); ylabel('\tau \cdot f_0');
title('Produto tempo x banda = 1 (constante)');

%% ---------------- Figura 5: 2o parametro - amplitude A ---------------
As = [0.5 1 2];
figure('Name','Variando A','Color','w','Position',[240 240 1100 600]);
for k = 1:numel(As)
    subplot(2,1,1); hold on; grid on;
    plot(t, As(k)*rectp(t), 'LineWidth', 1.8, 'DisplayName', sprintf('A = %.1f', As(k)));
    subplot(2,1,2); hold on; grid on;
    plot(f, As(k)*sincn(f), 'LineWidth', 1.8, 'DisplayName', sprintf('A = %.1f', As(k)));
end
subplot(2,1,1); xlim([-2 2]); ylim([-0.1 2.3]); xlabel('t [s]'); ylabel('x(t)');
title('Variando a amplitude A (\tau = 1 s fixo)'); legend;
subplot(2,1,2); xlim([-5 5]); xlabel('f [Hz]'); ylabel('X(f)');
title('A escala o espectro inteiro, mas os nulos NAO se movem (linearidade)'); legend;

%% ---------------- Figura 6: 3o parametro - atraso t0 -----------------
t0s = [0 0.5 1.5];
figure('Name','Variando t0','Color','w','Position',[280 280 1100 750]);
for k = 1:numel(t0s)
    x  = rectp(t - t0s(k));
    Xn = fftshift(fft(ifftshift(x))) / fs;
    subplot(3,1,1); hold on; grid on;
    plot(t, x, 'LineWidth', 1.8, 'DisplayName', sprintf('t_0 = %.1f s', t0s(k)));
    subplot(3,1,2); hold on; grid on;
    plot(f, abs(Xn), 'LineWidth', 1.8, 'DisplayName', sprintf('t_0 = %.1f s', t0s(k)));
    subplot(3,1,3); hold on; grid on;
    plot(f, angle(Xn), 'LineWidth', 1.2, 'DisplayName', sprintf('t_0 = %.1f s', t0s(k)));
end
subplot(3,1,1); xlim([-2 3]); ylim([-0.1 1.3]); xlabel('t [s]'); ylabel('x(t)');
title('Deslocando a porta: x(t - t_0)'); legend;
subplot(3,1,2); xlim([-5 5]); xlabel('f [Hz]'); ylabel('|X(f)|');
title('Modulo identico para todos os t_0'); legend;
subplot(3,1,3); xlim([-3 3]); xlabel('f [Hz]'); ylabel('\angle X(f) [rad]');
title('So a fase muda: e^{-j2\pi f t_0} (fase linear)'); legend;

%% ---------------- Verificacao de Parseval ----------------------------
fprintf('\n===== Verificacao de Parseval (energia) =====\n');
fprintf('  tau [s] | E tempo (A^2*tau) | E frequencia (integral |X|^2)\n');
for k = 1:numel(taus)
    x  = A * rectp(t / taus(k));
    Xn = fftshift(fft(ifftshift(x))) / fs;
    Et = sum(abs(x).^2) / fs;
    Ef = sum(abs(Xn).^2) * (fs/N);
    fprintf('   %4.1f   |      %7.4f       |      %7.4f\n', taus(k), Et, Ef);
end

%% ---------------- Modo interativo (sliders) --------------------------
rect_sinc_interativo();

%% =====================================================================
%  Funcoes locais
%% =====================================================================
function y = rectp(u)
%RECTP  Porta unitaria: 1 para |u|<1/2, 1/2 nas bordas, 0 fora.
    y = double(abs(u) < 0.5);
    y(abs(abs(u) - 0.5) < 1e-12) = 0.5;
end

function y = sincn(u)
%SINCN  sinc normalizada sin(pi*u)/(pi*u), sem depender de toolbox.
    y = ones(size(u));
    nz = (u ~= 0);
    y(nz) = sin(pi*u(nz)) ./ (pi*u(nz));
end

function rect_sinc_interativo()
%RECT_SINC_INTERATIVO  Janela com sliders para tau, A e t0.
    fs = 1000; T = 40;
    t  = -T/2 : 1/fs : T/2 - 1/fs;
    N  = numel(t);
    f  = (-N/2 : N/2-1) * (fs/N);

    fig = figure('Name','Rect x Sinc - INTERATIVO','Color','w', ...
                 'Position',[320 120 1100 760]);
    ax1 = axes('Parent',fig,'Position',[0.08 0.62 0.86 0.32]);
    ax2 = axes('Parent',fig,'Position',[0.08 0.20 0.86 0.32]);

    sTau = mkSlider(fig, 'tau [s]', 0.1, 5,  1,   0.10);
    sA   = mkSlider(fig, 'A',       0.1, 3,  1,   0.06);
    sT0  = mkSlider(fig, 't0 [s]', -3,   3,  0,   0.02);

    atualizar();

    function s = mkSlider(parent, rotulo, vmin, vmax, v0, ypos)
        uicontrol(parent,'Style','text','Units','normalized', ...
            'Position',[0.08 ypos 0.08 0.03],'String',rotulo, ...
            'BackgroundColor','w','HorizontalAlignment','right');
        s = uicontrol(parent,'Style','slider','Units','normalized', ...
            'Position',[0.18 ypos 0.62 0.03],'Min',vmin,'Max',vmax,'Value',v0, ...
            'Callback',@(~,~) atualizar());
        addlistener(s,'ContinuousValueChange',@(~,~) atualizar());
    end

    function atualizar()
        tau = sTau.Value; A = sA.Value; t0 = sT0.Value;
        x   = A * rectp((t - t0) / tau);
        Xn  = fftshift(fft(ifftshift(x))) / fs;
        Xa  = A * tau * sincn(f * tau);

        plot(ax1, t, x, 'b', 'LineWidth', 2); grid(ax1,'on');
        xlim(ax1,[-6 6]); ylim(ax1,[-0.1 3.3]);
        xlabel(ax1,'t [s]'); ylabel(ax1,'x(t)');
        title(ax1, sprintf('x(t) = %.2f \\cdot rect((t - %.2f)/%.2f)   |   energia = %.3f', ...
              A, t0, tau, A^2*tau));

        plot(ax2, f, abs(Xn), 'b', 'LineWidth', 2); hold(ax2,'on');
        plot(ax2, f, abs(Xa), 'r--', 'LineWidth', 1.2);
        xline(ax2,  1/tau, 'k:'); xline(ax2, -1/tau, 'k:');
        hold(ax2,'off'); grid(ax2,'on');
        ylim(ax2,[0 A*tau*1.1]);
        xlabel(ax2,'f [Hz]'); ylabel(ax2,'|X(f)|');
        legend(ax2,{'|FFT|','|A\tau sinc(f\tau)|','nulos \pm1/\tau'},'Location','northeast');
        title(ax2, sprintf('|X(f)|: pico = A\\tau = %.2f  |  1o nulo = 1/\\tau = %.2f Hz', ...
              A*tau, 1/tau));
    end
end
