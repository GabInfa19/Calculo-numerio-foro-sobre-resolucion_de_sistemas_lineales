% =====================================================================
% ESTIMACION ESPECTRAL CLASICA: EFECTO DE LA VENTANA (Octave / MATLAB)
% =====================================================================
% Idea general de la practica:
%   - Una "ventana" es una curva que se multiplica por la senal antes de
%     calcular su espectro. Cambia dos cosas:
%       * la RESOLUCION (que tan cerca pueden estar dos frecuencias para
%         verse separadas), que depende del ancho del lobulo principal;
%       * la FUGA espectral (energia que se "derrama" a frecuencias
%         vecinas), que depende del nivel de los lobulos secundarios.
%   - El PERIODOGRAMA estima la densidad espectral de potencia (DEP) de
%     una senal con la FFT. Es muy ruidoso.
%   - El PERIODOGRAMA PROMEDIO (Welch/Bartlett) divide la senal en
%     trozos, calcula el periodograma de cada uno y los promedia, lo que
%     reduce el ruido de la estimacion.
%
% Version corregida del codigo del grupo. Cambios marcados con [CAMBIO].
% =====================================================================

% Se limpian graficas, variables y Consola de Comando
clear; clc; close all;

% Paquete de senales para OCTAVE (Si usas MatLab comentarlo).
% En MATLAB, rectwin/hann/hamming/blackman/periodogram/pwelch vienen del
% Signal Processing Toolbox.
pkg load signal;

% [CAMBIO] Semilla fija: asi el "ruido aleatorio" es el mismo en cada
% corrida y las figuras se pueden repetir (en MATLAB tambien sirve rng(1)).
randn('state', 1);

% ---------------------------------------------------------------------
% PARAMETROS
% ---------------------------------------------------------------------
N = 1024;        % Numero de muestras (longitud de las ventanas y de la senal).
                 % Con N grande los lobulos de las ventanas se ven bien
                 % definidos; con N pequeno (p. ej. 128) el lobulo principal
                 % es mas ancho y la resolucion es peor.
Nfft = 4096;     % Puntos de la FFT. Si Nfft > N, la FFT rellena con ceros
                 % ("zero-padding"): solo INTERPOLA la curva (se ve mas
                 % suave), no mejora la resolucion real.

% ---------------------------------------------------------------------
% PASO 1: VENTANAS Y SU ESPECTRO
% ---------------------------------------------------------------------
% Se generan las cuatro ventanas como vectores columna de N puntos.
wRect = rectwin(N);          % Rectangular: todo 1. Lobulo principal mas
                             %   estrecho, pero lobulos laterales altos (-13 dB).
wHann = hann(N);             % Hanning: lobulo lateral ~ -31 dB.
wHamm = hamming(N);          % Hamming: lobulo lateral ~ -42 dB.
% [CAMBIO] En Octave blackman(N) devuelve -1.4e-17 en los extremos por
% redondeo, y pwelch rechaza ventanas con valores negativos. max(.,0) lo
% corrige sin cambiar la ventana en la practica.
wBlck = max(blackman(N), 0); % Blackman: lobulo lateral ~ -58 dB, pero
                             %   lobulo principal mas ancho.

% Espectro de cada ventana, en tres pasos:
%  1. fft(w, Nfft): Transformada Rapida de Fourier con Nfft puntos.
%  2. fftshift(): reordena el resultado para que la frecuencia 0 quede en
%     el centro de la grafica (la FFT la entrega al inicio).
%  3. abs(): magnitud, porque la FFT entrega numeros complejos.
WfRect = abs(fftshift(fft(wRect, Nfft)));
WfHann = abs(fftshift(fft(wHann, Nfft)));
WfHamm = abs(fftshift(fft(wHamm, Nfft)));
WfBlck = abs(fftshift(fft(wBlck, Nfft)));

% Eje de frecuencia normalizada de -pi a pi (rad/muestra), con un punto por
% cada punto de la FFT.
fVentana = linspace(-pi, pi, Nfft);

figure('Name', 'Paso 1: Espectro de las Ventanas');
plot(fVentana, WfRect); hold on;
plot(fVentana, WfHann);
plot(fVentana, WfHamm);
plot(fVentana, WfBlck);
xlim([-0.05 0.05]);              % Zoom: solo se ve el lobulo principal.
                                 % Se nota que Rectangular es la mas
                                 % estrecha y Blackman la mas ancha.
title('Caracteristicas Espectrales de las Ventanas');
xlabel('Frecuencia (rad/muestra)');   % [CAMBIO] antes decia (\pi rad/muestra),
                                      % pero el eje esta en radianes.
ylabel('Magnitud');
legend('Rectangular', 'Hanning', 'Hamming', 'Blackman');
grid on;
hold off;

% [CAMBIO, opcional] Mismo espectro pero en dB y normalizado a su maximo.
% La escala lineal no deja ver los lobulos secundarios (son muy pequenos);
% en dB si se lee su nivel respecto al lobulo principal, que es lo que pide
% el enunciado. El 1e-12 evita calcular log10(0).
figure('Name', 'Paso 1b: Espectro de las Ventanas en dB');
Wf = {WfRect, WfHann, WfHamm, WfBlck};
for i = 1:4
    plot(fVentana, 20*log10(Wf{i}/max(Wf{i}) + 1e-12)); hold on;
end
xlim([-0.2 0.2]); ylim([-120 5]);
title('Espectro de las Ventanas (dB, normalizado)');
xlabel('Frecuencia (rad/muestra)'); ylabel('Magnitud (dB)');
legend('Rectangular', 'Hanning', 'Hamming', 'Blackman');
grid on; hold off;

% ---------------------------------------------------------------------
% EXPERIMENTO 1.a: DOS SINUSOIDES IGUALES + RUIDO BLANCO
% ---------------------------------------------------------------------
% Pregunta: cada ventana logra ver las dos frecuencias por separado?
n = 0:N-1;               % indices de las muestras

% Frecuencias angulares dadas en el enunciado (rad/muestra). Estan muy
% cerca entre si: w1 = 0.5236 y w2 = 0.4488.
w1 = 2*pi/12;
w2 = 2*pi/14;

% Ruido blanco: media 0 y varianza 1. randn genera numeros con
% distribucion normal estandar, uno por muestra.
ruido = randn(1, N);

% Senal: dos sinusoides de amplitud 1 mas el ruido. La MISMA realizacion de
% ruido se usa para las cuatro ventanas, asi las diferencias entre graficas
% se deben solo a la ventana.
A1a = 1;
A2a = 1;
xa = A1a*sin(w1*n) + A2a*sin(w2*n) + ruido;

figure('Name', 'Experimento 1.a');
ventanas = {wRect, wHann, wHamm, wBlck};   % se agrupan para usar un ciclo for
titulos_a = {'1.a Rectangular', '1.a Hanning', '1.a Hamming', '1.a Blackman'};

for i = 1:4
    subplot(2, 2, i);
    % periodogram(x, ventana, Nfft): multiplica x por la ventana, calcula
    % la FFT y normaliza. Devuelve la DEP (Pxx) y la frecuencia w en
    % rad/muestra de 0 a pi.
    [Pxx, w] = periodogram(xa, ventanas{i}, Nfft);
    plot(w, Pxx);
    xlim([0 0.3*pi]);            % zoom a la zona donde estan w1 y w2
    title(titulos_a{i});
    xlabel('Frecuencia (rad/muestra)');
    ylabel('Potencia');
    grid on;
end

% ---------------------------------------------------------------------
% EXPERIMENTO 1.b: UNA SINUSOIDE FUERTE Y OTRA MUY DEBIL + RUIDO
% ---------------------------------------------------------------------
% Pregunta: que ventana deja ver la sinusoide debil (amplitud 0.01, es
% decir 40 dB por debajo de la fuerte)? Ojo: con ruido de varianza 1 la
% sinusoide debil queda por debajo del ruido y no se ve con ninguna ventana.
A1b = 1;
A2b = 0.01;
xb = A1b*sin(w1*n) + A2b*sin(w2*n) + ruido;   % mismo ruido que en 1.a

figure('Name', 'Experimento 1.b');
titulos_b = {'1.b Rectangular', '1.b Hanning', '1.b Hamming', '1.b Blackman'};

for i = 1:4
    subplot(2, 2, i);
    [Pxx, w] = periodogram(xb, ventanas{i}, Nfft);
    plot(w, Pxx);
    xlim([0 0.3*pi]);
    title(titulos_b{i});
    xlabel('Frecuencia (rad/muestra)');
    ylabel('Potencia');
    grid on;
end

% ---------------------------------------------------------------------
% EXPERIMENTO 2: PERIODOGRAMA PROMEDIO DE RUIDO BLANCO
% ---------------------------------------------------------------------
% Pregunta: cuanto se calma la estimacion al promediar varios
% periodogramas? La DEP verdadera del ruido blanco es PLANA, asi que una
% buena estimacion debe verse como una linea casi horizontal.
muestrasRuido = N * 16;                % ruido 16 veces mas largo que una ventana
ruidoLargo = randn(1, muestrasRuido);  % -> 16 segmentos de N muestras

% pwelch(x, ventana, solape, Nfft, Fs): corta x en segmentos del largo de la
% ventana, calcula el periodograma de cada uno y los promedia.
%   - solape = 0: segmentos sin solape (promedio de Bartlett).
%   - [CAMBIO] Fs = 2*pi como 5to argumento. Sin el, Octave usa Fs = 1: el eje
%     sale en ciclos/muestra (0 a 0.5) y no en rad/muestra, y el nivel de
%     potencia queda ~2 en vez de ~0.32, distinto al de periodogram().
%     Con Fs = 2*pi coincide con periodogram() y con MATLAB.
[PxxRect, wProm] = pwelch(ruidoLargo, wRect, 0, Nfft, 2*pi);
[PxxHann, wProm] = pwelch(ruidoLargo, wHann, 0, Nfft, 2*pi);
[PxxHamm, wProm] = pwelch(ruidoLargo, wHamm, 0, Nfft, 2*pi);
[PxxBlck, wProm] = pwelch(ruidoLargo, wBlck, 0, Nfft, 2*pi);

figure('Name', 'Experimento 2: Periodograma Promedio');
Pxx_promedios = {PxxRect, PxxHann, PxxHamm, PxxBlck};
titulos_prom = {'2. Rectangular (Promedio)', '2. Hanning (Promedio)', '2. Hamming (Promedio)', '2. Blackman (Promedio)'};

for i = 1:4
    subplot(2, 2, i);
    plot(wProm, Pxx_promedios{i});
    title(titulos_prom{i});
    xlabel('Frecuencia (rad/muestra)');
    ylabel('Densidad Espectral');
    grid on;
end

% ---------------------------------------------------------------------
% VERSIONES EN dB
% ---------------------------------------------------------------------
% En dB (10*log10 de la potencia) se pueden ver a la vez cosas muy grandes
% y muy pequenas, por ejemplo los picos y el piso de ruido.

% Experimento 1.b en dB
figure('Name', 'Experimento 1.b en dB');
titulos_b_dB = {'1.b Rectangular (dB)', '1.b Hanning (dB)', '1.b Hamming (dB)', '1.b Blackman (dB)'};

for i = 1:4
    subplot(2, 2, i);
    [Pxx, w] = periodogram(xb, ventanas{i}, Nfft);
    plot(w, 10*log10(Pxx));      % potencia a escala logaritmica (dB)
    xlim([0 0.3*pi]);
    title(titulos_b_dB{i});
    xlabel('Frecuencia (rad/muestra)');
    ylabel('Potencia (dB)');
    grid on;
end

% Experimento 2 en dB (reutiliza los datos de pwelch ya calculados)
figure('Name', 'Experimento 2: Periodograma Promedio en dB');
titulos_prom_dB = {'2. Rectangular Promedio (dB)', '2. Hanning Promedio (dB)', '2. Hamming Promedio (dB)', '2. Blackman Promedio (dB)'};

for i = 1:4
    subplot(2, 2, i);
    plot(wProm, 10*log10(Pxx_promedios{i}));
    title(titulos_prom_dB{i});
    xlabel('Frecuencia (rad/muestra)');
    ylabel('Densidad Espectral (dB)');
    grid on;
end
