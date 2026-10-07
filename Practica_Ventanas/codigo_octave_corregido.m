% Estimacion espectral clasica: efecto de la ventana (Octave / MATLAB)
% Version corregida del codigo del grupo. Cambios marcados con [CAMBIO].

% Se limpian graficas, variables y Consola de Comando
clear; clc; close all;

% Paquete de senales para OCTAVE (Si usas MatLab comentarlo)
pkg load signal;

% [CAMBIO] Semilla fija para que los resultados se puedan repetir
% (en MATLAB tambien funciona rng(1))
randn('state', 1);

N = 1024;        % Numero de muestras. Valor alto para ver los lobulos bien definidos.
Nfft = 4096;     % Puntos para la FFT. Mayor que N hace "zero-padding" (relleno de ceros).

% Generacion de las ventanas en tiempo
wRect = rectwin(N);          % Ventana Rectangular
wHann = hann(N);             % Ventana de Hanning
wHamm = hamming(N);          % Ventana de Hamming
% [CAMBIO] En Octave blackman(N) devuelve -1.4e-17 en los extremos por redondeo
% y pwelch rechaza ventanas con valores negativos. max(.,0) lo corrige.
wBlck = max(blackman(N), 0); % Ventana de Blackman

% Calculo del espectro de frecuencias
% fft(): Transformada Rapida de Fourier. fftshift(): centra la frecuencia 0.
% abs(): magnitud (la FFT arroja numeros complejos).
WfRect = abs(fftshift(fft(wRect, Nfft)));
WfHann = abs(fftshift(fft(wHann, Nfft)));
WfHamm = abs(fftshift(fft(wHamm, Nfft)));
WfBlck = abs(fftshift(fft(wBlck, Nfft)));

% Eje de frecuencia: -pi a pi (rad/muestra)
fVentana = linspace(-pi, pi, Nfft);

figure('Name', 'Paso 1: Espectro de las Ventanas');
plot(fVentana, WfRect); hold on;
plot(fVentana, WfHann);
plot(fVentana, WfHamm);
plot(fVentana, WfBlck);
xlim([-0.05 0.05]);              % Zoom al lobulo principal
title('Caracteristicas Espectrales de las Ventanas');
xlabel('Frecuencia (rad/muestra)');   % [CAMBIO] antes decia (\pi rad/muestra), pero el eje esta en rad
ylabel('Magnitud');
legend('Rectangular', 'Hanning', 'Hamming', 'Blackman');
grid on;
hold off;

% [CAMBIO, opcional] Espectro en dB normalizado: permite leer el nivel del
% lobulo secundario, que el enunciado pide y la escala lineal no muestra.
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

% 1.A (Amplitudes iguales A1 = 1, A2 = 1 + Ruido Blanco)
n = 0:N-1;               % muestras

% Frecuencias angulares dadas en el enunciado
w1 = 2*pi/12;
w2 = 2*pi/14;

% Ruido blanco de media 0 y varianza unitaria (randn: normal estandar)
ruido = randn(1, N);

% Senal del experimento 1.a
A1a = 1;
A2a = 1;
xa = A1a*sin(w1*n) + A2a*sin(w2*n) + ruido;

figure('Name', 'Experimento 1.a');
ventanas = {wRect, wHann, wHamm, wBlck};
titulos_a = {'1.a Rectangular', '1.a Hanning', '1.a Hamming', '1.a Blackman'};

for i = 1:4
    subplot(2, 2, i);
    % periodogram() devuelve la DEP (Pxx) y la frecuencia w en rad/muestra (0 a pi)
    [Pxx, w] = periodogram(xa, ventanas{i}, Nfft);
    plot(w, Pxx);
    xlim([0 0.3*pi]);            % zoom para ver w1 y w2
    title(titulos_a{i});
    xlabel('Frecuencia (rad/muestra)');
    ylabel('Potencia');
    grid on;
end

% 1.B (Amplitudes A1 = 1, A2 = 0.01 + Ruido Blanco)
A1b = 1;
A2b = 0.01;
xb = A1b*sin(w1*n) + A2b*sin(w2*n) + ruido;

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

% 2.- Periodograma promedio para ruido blanco
% Secuencia 16 veces mas larga para poder promediar (K = 16 segmentos)
muestrasRuido = N * 16;
ruidoLargo = randn(1, muestrasRuido);

% Periodogramas promediados (solapamiento = 0)
% [CAMBIO] Se agrega Fs = 2*pi como 5to argumento. Sin el, Octave usa Fs = 1:
% el eje sale en ciclos/muestra (0 a 0.5) y no en rad/muestra, y el nivel de
% potencia queda ~2 en vez de ~0.32, distinto al de periodogram().
% Con Fs = 2*pi coincide con periodogram() y con MATLAB.
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

% Experimento 2 en dB
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
