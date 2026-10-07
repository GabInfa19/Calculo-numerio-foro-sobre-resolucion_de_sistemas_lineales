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
% =====================================================================

% Se limpian graficas, variables y Consola de Comando
clear; clc; close all;

% Paquete de senales para OCTAVE (Si usas MatLab comentarlo).
% En MATLAB, rectwin/hann/hamming/blackman/periodogram/pwelch vienen del
% Signal Processing Toolbox.
pkg load signal;

% Semilla fija: asi el "ruido aleatorio" es el mismo en cada
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
% En Octave blackman(N) devuelve -1.4e-17 en los extremos por
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
xlabel('Frecuencia (rad/muestra)');   % antes decia (\pi rad/muestra),
                                      % pero el eje esta en radianes.
ylabel('Magnitud');
legend('Rectangular', 'Hanning', 'Hamming', 'Blackman');
grid on;
hold off;

% Mismo espectro pero en dB y normalizado a su maximo.
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
%   - Fs = 2*pi como 5to argumento. Sin el, Octave usa Fs = 1: el eje
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

% ---------------------------------------------------------------------
% EXPERIMENTO 1.b SIN RUIDO
% ---------------------------------------------------------------------
% Con ruido de varianza 1, la sinusoide debil (A2 = 0.01) queda ~80 dB por
% debajo de la fuerte y el ruido se la "come": no se ve con ninguna ventana.
% Para aislar el efecto de la VENTANA se repite el experimento SIN ruido. Asi
% lo unico que puede tapar a la sinusoide debil es la FUGA espectral de la
% sinusoide fuerte (los lobulos secundarios de la ventana).
xc = A1b*sin(w1*n) + A2b*sin(w2*n);   % sin ruido
xs = A1b*sin(w1*n);                   % solo la fuerte (para medir la fuga)

figure('Name', 'Experimento 1.b sin ruido en dB');
titulos_c = {'1.b sin ruido Rectangular (dB)', '1.b sin ruido Hanning (dB)', '1.b sin ruido Hamming (dB)', '1.b sin ruido Blackman (dB)'};
nombres = {'Rectangular', 'Hanning', 'Hamming', 'Blackman'};

for i = 1:4
    subplot(2, 2, i);
    [Pc, w] = periodogram(xc, ventanas{i}, Nfft);   % fuerte + debil
    [Ps, w] = periodogram(xs, ventanas{i}, Nfft);   % solo fuerte
    plot(w, 10*log10(Pc)); hold on;
    plot(w, 10*log10(Ps), 'r--');     % linea roja: fuga de la fuerte sola
    xlim([0 0.3*pi]); ylim([-120 20]);
    title(titulos_c{i});
    xlabel('Frecuencia (rad/muestra)');
    ylabel('Potencia (dB)');
    grid on; hold off;

    % Margen: cuantos dB sobresale el pico de la debil (cerca de w2) por
    % encima del MAXIMO de la fuga de la fuerte en esa misma zona (la fuga
    % oscila entre picos y valles, por eso se compara con su maximo).
    % Mientras mas grande, mas facil es ver la sinusoide debil.
    cerca = find(abs(w - w2) < 0.01);
    margen = 10*log10(max(Pc(cerca))) - 10*log10(max(Ps(cerca)));
    printf('%s: la sinusoide debil sobresale %.1f dB sobre la fuga\n', nombres{i}, margen);
end
legend('con la debil', 'solo la fuerte (fuga)');

% ---------------------------------------------------------------------
% VALORES NUMERICOS PARA EL INFORME
% ---------------------------------------------------------------------
% Imprime en la consola los numeros que se citan en el informe, asi se
% pueden reproducir con este mismo codigo.

% (a) Caracteristicas de las ventanas, medidas en su espectro. Se usa una FFT
% con mucho relleno de ceros (2^18 puntos) para leer bien los lobulos.
%   - ancho nulo a nulo: distancia entre los dos primeros ceros del lobulo
%     principal, en multiplos de 1/N ciclos/muestra.
%   - ancho a -3 dB: ancho donde la magnitud cae a 0.707 del maximo.
%   - lobulo secundario: mayor pico fuera del lobulo principal, en dB
%     respecto al pico principal.
NF = 2^18;
printf('\n--- Caracteristicas de las ventanas (N = %d) ---\n', N);
for i = 1:4
    Wm = abs(fft(ventanas{i}, NF));
    Wm = Wm(1:NF/2) / max(Wm);              % de 0 a pi, normalizado
    fc = (0:NF/2-1) / NF;                   % ciclos/muestra
    k = 2;
    while ~(Wm(k) < Wm(k-1) && Wm(k) <= Wm(k+1))   % primer nulo
        k = k + 1;
    end
    k3 = find(Wm < 10^(-3/20), 1);          % primer punto bajo -3 dB
    printf('%-12s nulo-nulo = %.2f/N | -3 dB = %.2f/N | lobulo sec. = %.1f dB\n', ...
        nombres{i}, 2*fc(k)*N, 2*fc(k3)*N, 20*log10(max(Wm(k:end))));
end

% (b) Experimento 1.b con ruido: piso de ruido lejos de las sinusoides
% (frecuencias mayores que 0.6*pi rad/muestra).
printf('\n--- Experimento 1.b con ruido: piso de ruido (w > 0.6*pi) ---\n');
for i = 1:4
    [Pxx, w] = periodogram(xb, ventanas{i}, Nfft);
    zona = 10*log10(Pxx(w > 0.6*pi));
    printf('%-12s media = %.1f dB | maximo = %.1f dB\n', nombres{i}, mean(zona), max(zona));
end

% (c) Experimento 2: cuanto se calma la estimacion al promediar.
% Se mide la VARIANZA RELATIVA = var/media^2 de la curva (para ruido blanco la
% DEP verdadera es constante, asi que toda variacion es error de estimacion).
% Es independiente de la escala. Un periodograma solo da ~1; promediar K
% segmentos sin solape deberia dar ~1/K.
K = muestrasRuido / N;
printf('\n--- Experimento 2: K = %d segmentos ---\n', K);
for i = 1:4
    P1 = periodogram(ruidoLargo(1:N), ventanas{i}, Nfft);   % un solo segmento
    Pp = Pxx_promedios{i};                                  % promedio de K segmentos
    printf('%-12s 1 segmento: media = %.3f, var.rel. = %.3f | promedio: media = %.3f, var.rel. = %.3f (1/K = %.3f)\n', ...
        nombres{i}, mean(P1), var(P1)/mean(P1)^2, mean(Pp), var(Pp)/mean(Pp)^2, 1/K);
end
