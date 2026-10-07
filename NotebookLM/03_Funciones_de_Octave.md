# Funciones de Octave usadas en el código de la práctica

Este documento explica cada función del archivo practica_ventanas.m (código de la práctica de estimación espectral clásica con ventanas, en Octave/MATLAB). Para cada función se indica qué hace, cómo se llama en el código y para qué se usa.

## Preparación del programa

- **clear; clc; close all;** Borra las variables, limpia la ventana de comandos y cierra las figuras abiertas. Se usa al inicio para empezar de cero.
- **pkg load signal;** Carga el paquete de procesamiento de señales de Octave, donde están periodogram, pwelch, hann, hamming, blackman y rectwin. En MATLAB no se usa (las mismas funciones vienen en el Signal Processing Toolbox).
- **randn('state', 1);** Fija la semilla del generador de números aleatorios. Así el ruido es el mismo en cada ejecución y los resultados (figuras y tablas) se pueden repetir. En MATLAB equivale a rng(1).

## Generar la señal

- **randn(1, N)** Genera N números aleatorios con distribución normal estándar (media 0 y varianza 1). Es el ruido blanco gaussiano. En el código: ruido = randn(1, N) para el experimento 1 y ruidoLargo = randn(1, N*16) para el experimento 2.
- **sin(x)** Calcula el seno de cada elemento de x, en radianes. En el código: sin(w1*n) y sin(w2*n) generan las dos sinusoides, con w1 = 2π/12 y w2 = 2π/14.
- **0:N-1** Operador de rango: crea el vector 0, 1, 2, …, N−1. En el código: n = 0:N-1 son los índices de las muestras.

## Ventanas

- **rectwin(N)** Ventana rectangular de N puntos: todos los valores valen 1.
- **hann(N)** Ventana de Hanning de N puntos: 0.5 − 0.5·cos(2πn/(N−1)); vale 0 en los extremos.
- **hamming(N)** Ventana de Hamming de N puntos: 0.54 − 0.46·cos(2πn/(N−1)); vale 0.08 en los extremos.
- **blackman(N)** Ventana de Blackman de N puntos: 0.42 − 0.5·cos(2πn/(N−1)) + 0.08·cos(4πn/(N−1)). En Octave devuelve −1.4·10⁻¹⁷ en los extremos por redondeo, por eso el código usa max(blackman(N), 0).

Las cuatro funciones devuelven un vector columna de N valores.

## Transformada de Fourier y espectro

- **fft(x, Nfft)** Transformada rápida de Fourier de x con Nfft puntos. Si Nfft es mayor que la longitud de x, rellena con ceros (zero-padding), lo cual interpola la curva sin mejorar la resolución. En el código: fft(wRect, Nfft) para obtener el espectro de la ventana.
- **fftshift(X)** Reordena el resultado de la FFT para que la frecuencia cero quede en el centro del gráfico, con las frecuencias negativas a la izquierda. La FFT entrega la frecuencia cero al principio del vector.
- **abs(X)** Magnitud de cada elemento. Se usa porque la FFT entrega números complejos.
- **linspace(-pi, pi, Nfft)** Crea un vector de Nfft puntos equiespaciados entre −π y π. Es el eje de frecuencia en rad/muestra para graficar el espectro de las ventanas.
- **log10(x)** Logaritmo en base 10. Se usa para pasar a dB: 10·log10(potencia) o 20·log10(magnitud). En el código se suma 1e-12 dentro del logaritmo del espectro de ventanas para evitar log10(0).

## Estimación espectral (el corazón de la práctica)

- **[Pxx, w] = periodogram(x, ventana, Nfft)** Calcula el periodograma modificado de x: multiplica la señal por la ventana, calcula la FFT de Nfft puntos y normaliza por la potencia de la ventana (Σw²). Devuelve Pxx, la densidad espectral de potencia unilateral (por rad/muestra), y w, el vector de frecuencias en rad/muestra, de 0 a π. Para ruido blanco de varianza 1, el nivel es 1/π ≈ 0.318.
- **[Pxx, w] = pwelch(x, ventana, solape, Nfft, Fs)** Calcula el periodograma promedio con el método de Welch. Corta x en segmentos del largo de la ventana, desplazados según el solape, multiplica cada segmento por la ventana, calcula su periodograma y promedia todos. Argumentos: x es la señal; ventana es el vector de la ventana (su longitud es el tamaño de cada segmento); solape es el número de muestras de solape entre segmentos vecinos (0 en el código); Nfft es el número de puntos de la FFT; Fs es la frecuencia de muestreo. En el código se usa Fs = 2π para que el eje esté en rad/muestra y el nivel sea comparable con el de periodogram; sin ese argumento, Octave usa Fs = 1. Número de segmentos: K = ⌊(N − L)/(L − solape)⌋ + 1; con 16384 muestras, L = 1024 y solape 0, K = 16.

## Gráficas

- **figure('Name', texto)** Abre una figura nueva con ese nombre.
- **plot(x, y)** Dibuja la curva y contra x. Con hold on se dibujan varias curvas en la misma figura; con hold off se termina.
- **subplot(2, 2, i)** Divide la figura en una cuadrícula de 2 filas por 2 columnas y dibuja en la posición i. Se usa para mostrar las cuatro ventanas en una sola figura.
- **xlim([a b]), ylim([a b])** Fijan el rango visible de cada eje (zoom). En el código: xlim([0 0.3*pi]) muestra la zona donde están las sinusoides.
- **title, xlabel, ylabel** Ponen el título de la gráfica y el nombre de cada eje.
- **legend('a', 'b', …)** Muestra la leyenda con el nombre de cada curva.
- **grid on** Muestra la cuadrícula.

## Cálculo y salida de resultados

- **max(x)** Máximo de un vector. Con dos argumentos, max(x, 0) compara elemento a elemento y deja el mayor; en el código recorta los valores negativos minúsculos de blackman.
- **mean(x), var(x)** Media y varianza de un vector. Se usan para medir el nivel y la varianza relativa (var/media²) de las estimaciones del experimento 2.
- **find(condición)** Devuelve los índices donde se cumple la condición. En el código: find(abs(w - w2) < 0.01) selecciona las frecuencias cercanas a w2, y find(Wm < 10^(-3/20), 1) encuentra el primer punto por debajo de −3 dB.
- **min(x)** Mínimo de un vector; con un índice de salida se usa para encontrar el punto más cercano a un valor.
- **printf('texto %.1f', valor)** Imprime en la consola texto con números formateados; se usa para mostrar las tablas de resultados.

## Sintaxis que conviene reconocer

- **ventanas = {wRect, wHann, wHamm, wBlck};** Es un arreglo de celdas (llaves). Se accede con ventanas{i}, lo que permite recorrer las cuatro ventanas con un ciclo for i = 1:4.
- **[Pxx, w] = …** Una función puede devolver más de una salida, y se guardan entre corchetes.
- **for i = 1:4 … end** Ciclo que repite el bloque con i = 1, 2, 3, 4.
- **while … end** Ciclo que repite mientras se cumpla una condición; se usa para buscar el primer nulo del lóbulo principal.
