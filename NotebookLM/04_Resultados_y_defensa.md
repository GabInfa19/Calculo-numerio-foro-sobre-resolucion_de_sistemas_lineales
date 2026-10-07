# Resultados de la práctica y preparación para la defensa

## Descripción de la práctica

Práctica 1 de Procesamiento de Señales (Escuela de Ingeniería en Telecomunicaciones, UCAB): estimación espectral clásica. Se compara el periodograma y el periodograma promedio usando cuatro ventanas: Rectangular, Hanning, Hamming y Blackman. Grupo: Gerardo Febres, Gabriel Infante y Aria Bahrami. Código en Octave/MATLAB.

La práctica pide: (1) investigar las características espectrales de las ventanas (ancho del lóbulo principal y nivel relativo del lóbulo secundario) y graficar su espectro; (2) experimento 1.a: sinusoides de amplitud A1 = A2 = 1 con frecuencias w1 = 2π/12 y w2 = 2π/14 más ruido blanco de media cero y varianza unitaria, obtener el periodograma con cada ventana; (3) experimento 1.b: igual con A1 = 1 y A2 = 0.01; (4) experimento 2: periodograma promedio de ruido blanco variando la ventana.

Parámetros del código final: N = 1024 muestras, Nfft = 4096, semilla fija randn('state', 1), ruido del experimento 2 de 16·N muestras (K = 16 segmentos sin solape), pwelch con Fs = 2π.

## Tabla 1. Características de las ventanas (N = 1024, medidas con el código)

| Ventana | Lóbulo principal (nulo a nulo) | Ancho a −3 dB | Lóbulo secundario (dB) | Decaimiento | ENBW (bins) | Ganancia coherente |
|---|---|---|---|---|---|---|
| Rectangular | 2.00/N | 0.89/N | −13.3 | 6 dB/octava | 1.00 | 1.00 |
| Hanning | 4.01/N | 1.45/N | −31.5 | 18 dB/octava | 1.50 | 0.50 |
| Hamming | 4.01/N | 1.30/N | −42.7 | 6 dB/octava | 1.36 | 0.54 |
| Blackman | 6.01/N | 1.65/N | −58.1 | 18 dB/octava | 1.73 | 0.42 |

Los anchos están en múltiplos de 1/N ciclos/muestra. La fórmula de las ventanas es w[n] = a0 − a1·cos(2πn/(N−1)) + a2·cos(4πn/(N−1)) con (a0, a1, a2) = (1, 0, 0), (0.5, 0.5, 0), (0.54, 0.46, 0) y (0.42, 0.5, 0.08).

## Experimento 1.a (A1 = A2 = 1, con ruido)

Con N = 1024 las cuatro ventanas muestran los dos picos separados, en w2 ≈ 0.449 y w1 ≈ 0.524 rad/muestra. La separación es Δf = 1/12 − 1/14 = 0.0119 ciclos/muestra = 12.2 bins. Las alturas teóricas del pico son A²·N/(4π·ENBW) = 81.5 (rectangular), 54.3 (Hanning), 59.8 (Hamming) y 47.1 (Blackman). Con N = 128 (prueba en Python) solo la rectangular separaba los picos.

## Experimento 1.b (A1 = 1, A2 = 0.01, con ruido)

Solo aparece el pico de la sinusoide fuerte con las cuatro ventanas. La débil está 40 dB por debajo de la fuerte ((0.01)² = 10⁻⁴). Su pico esperado es ≈ −21 dB y el piso de ruido está en ≈ −7 dB: queda unos 14 dB por debajo del piso.

Tabla 2. Piso de ruido (frecuencias w > 0.6π rad/muestra):

| Ventana | Piso medio (dB) | Máximo (dB) |
|---|---|---|
| Rectangular | −7.3 | 3.5 |
| Hanning | −7.1 | 2.5 |
| Hamming | −7.1 | 2.6 |
| Blackman | −7.1 | 2.0 |

## Experimento 1.b sin ruido (complemento)

Tabla 3. Cuánto sobresale el pico de la sinusoide débil sobre el máximo de la fuga de la fuerte en esa zona:

| Ventana | Margen (dB) | Resultado |
|---|---|---|
| Rectangular | 0.6 | tapada por la fuga |
| Hamming | 9.7 | visible, con poco margen |
| Hanning | 31.5 | claramente visible |
| Blackman | 39.3 | claramente visible |

## Experimento 2 (ruido blanco, K = 16 segmentos de 1024 muestras, sin solape)

Tabla 4. Media y varianza relativa (var/media²) de la estimación de la DEP (valor verdadero 1/π = 0.318):

| Ventana | Media, 1 segmento | Var. relativa, 1 segmento | Media, promedio | Var. relativa, promedio | 1/K |
|---|---|---|---|---|---|
| Rectangular | 0.320 | 0.930 | 0.323 | 0.057 | 0.062 |
| Hanning | 0.295 | 0.844 | 0.320 | 0.056 | 0.062 |
| Hamming | 0.296 | 0.848 | 0.321 | 0.056 | 0.062 |
| Blackman | 0.290 | 0.855 | 0.320 | 0.057 | 0.062 |

Prueba adicional en Python (4096 muestras, L = 128): el solape del 50 % redujo la varianza de ≈ 0.027 a ≈ 0.013 con Hanning, Hamming y Blackman, y de 0.028 a 0.022 con la rectangular.

## Problemas encontrados y cómo se resolvieron

1. pwelch fallaba con Blackman ("window vector must be real and >=0"): blackman(N) devuelve −1.4·10⁻¹⁷ en los extremos en Octave. Solución: max(blackman(N), 0).
2. Eje y nivel distintos entre periodogram y pwelch: sin Fs, pwelch usa Fs = 1 (eje 0 a 0.5, nivel ≈ 2). Solución: pasar 2*pi como quinto argumento (eje 0 a π, nivel ≈ 0.32).
3. Eje rotulado "(π rad/muestra)" cuando estaba en radianes. Solución: "(rad/muestra)".
4. No se veía el nivel de los lóbulos secundarios en escala lineal. Solución: gráfica extra en dB normalizada.
5. Resultados distintos en cada corrida. Solución: semilla fija.
6. Error de cálculo: se dijo que la sinusoide débil estaba 80 dB por debajo; son 40 dB.
7. Medición del margen de la débil con la rectangular daba 12.8 dB con un solo punto de la fuga, que no coincidía con la figura. Solución: comparar con el máximo de la fuga en la zona (0.6 dB).
8. Con N = 128 solo la rectangular separaba las sinusoides. Se usó N = 1024 (igual a una práctica de 2022 de otro grupo).

## Guion de defensa (≈ 2 minutos)

Nuestra práctica compara el periodograma y el periodograma promedio usando cuatro ventanas: Rectangular, Hanning, Hamming y Blackman. Trabajamos en Octave con N = 1024 muestras y una FFT de 4096 puntos.

Primero estudiamos las ventanas. La rectangular tiene el lóbulo principal más estrecho, 2/N, pero el lóbulo secundario más alto, −13 dB. Hanning y Hamming tienen un lóbulo principal de 4/N y laterales de −31 y −43 dB. Blackman tiene el más ancho, 6/N, y los laterales más bajos, −58 dB. Es un compromiso: lo que se gana en fuga espectral se paga en resolución.

En el experimento 1.a generamos dos sinusoides de amplitud 1, en 2π/12 y 2π/14, más ruido blanco. Con N = 1024 las cuatro ventanas separan los dos picos, porque la separación equivale a 12 bins y es mayor que el lóbulo más ancho. Con N = 128 solo la rectangular los separaba.

En el 1.b la segunda sinusoide tiene amplitud 0.01, es decir 40 dB por debajo de la primera. Con ruido de varianza 1 no se ve con ninguna ventana, porque su pico queda unos 14 dB bajo el piso de ruido. Para ver el efecto de la ventana repetimos el experimento sin ruido: ahí Blackman y Hanning la dejan ver, con 39 y 31 dB de margen, Hamming apenas con 10 dB y la rectangular la tapa. Esto se debe a que Hanning y Blackman bajan sus lóbulos laterales 18 dB por octava, y la rectangular y Hamming solo 6.

En el experimento 2 promediamos 16 segmentos de ruido blanco con pwelch. El nivel medio fue 0.32, que es 1/π, el valor teórico. La varianza relativa bajó de casi 1 en un solo periodograma a 0.057, cerca de 1/K = 0.0625. El precio del promedio es la resolución, porque cada segmento es más corto.

En conclusión: la ventana ataca la fuga espectral, el promedio ataca la varianza, y ninguna de las dos reemplaza tener una mejor relación señal a ruido.

## Preguntas probables en la defensa, con respuesta

**¿Por qué usaron Nfft = 4096 si N = 1024?** Para interpolar la curva (zero-padding): se ve más suave y los picos se localizan mejor. No mejora la resolución real, que depende de N y de la ventana.

**¿Por qué el periodograma simple no es consistente?** Porque en cada frecuencia es una variable exponencial con desviación estándar igual a su media (varianza relativa ≈ 1), sin importar N. Más muestras dan más puntos de frecuencia, no menos varianza por punto.

**¿Qué ventana elegirían para separar dos frecuencias muy cercanas y de amplitud parecida?** La rectangular, por su lóbulo principal más estrecho.

**¿Y para detectar una componente débil junto a una fuerte?** Hanning o Blackman, por sus lóbulos laterales bajos que caen 18 dB por octava. Siempre que la componente débil quede por encima del ruido.

**¿Por qué en 1.b con ruido no se ve la sinusoide débil con ninguna ventana?** Porque está 40 dB por debajo de la fuerte y su pico queda unos 14 dB bajo el nivel medio del ruido. La ventana solo cambia la fuga, no puede subir un pico por encima del ruido.

**¿Por qué el piso de ruido sale en −7 dB y no en −5 dB?** El periodograma de ruido es exponencial en cada frecuencia, y el promedio del logaritmo de una exponencial queda 2.5 dB bajo el logaritmo de su promedio (−5.0 − 2.5 ≈ −7.5 dB).

**¿Por qué pwelch lleva 2*pi?** Sin la frecuencia de muestreo, Octave usa Fs = 1 y el eje y el nivel no coinciden con los de periodogram.

**¿Por qué max(blackman(N), 0)?** Octave devuelve un valor negativo del orden de −1.4·10⁻¹⁷ en los extremos y pwelch rechaza ventanas con valores negativos.

**¿Cómo cambiarían K?** Cambiando el factor 16 en muestrasRuido = N * 16. Con solape distinto de cero, K = ⌊(N − L)/(L − solape)⌋ + 1.

**¿Qué pasa con la resolución si se promedia?** Empeora, porque cada segmento tiene L muestras y el lóbulo principal es proporcional a 1/L. Con L = 64 el ancho a −3 dB de la rectangular (0.0139) supera la separación entre las sinusoides (0.0119) y ya no se separan.

**¿Qué diferencia hay entre Bartlett y Welch?** Bartlett usa segmentos sin solape y ventana rectangular; Welch permite solape (típicamente 50 %) y usa una ventana, así que sus periodogramas son modificados. pwelch implementa Welch.

**¿Por qué se normaliza por Σw²?** Porque la DEP es una potencia, cuadrática en la señal; Σw² mide la potencia que la ventana deja pasar de un ruido blanco. Así el nivel estimado no depende de la ventana.

**¿Qué hicieron distinto a las prácticas de cursos anteriores?** Compararon ventanas (las de 2022 usaban solo la rectangular), usaron escala en dB, midieron valores en tablas y explicaron por qué sucede cada resultado.
