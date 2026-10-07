# Por qué sucede cada cosa en la práctica (y gracias a qué)

Este documento explica, para cada resultado de la práctica de estimación espectral clásica con ventanas, tres cosas: qué se observó, por qué sucede y gracias a qué propiedad o teoría sucede. Los valores numéricos fueron medidos ejecutando el código de la práctica en Octave (N = 1024 muestras, Nfft = 4096, semilla fija) o calculados con las fórmulas indicadas.

Notación: N es el número de muestras; Nfft es el tamaño de la FFT (con relleno de ceros); ω es la frecuencia en rad/muestra y f = ω/2π en ciclos/muestra; "bin" es la separación entre puntos de la DFT, 1/N ciclos/muestra; DEP es la densidad espectral de potencia; dB significa 10·log10 para potencias y 20·log10 para amplitudes.

---

## 1. El ancho del lóbulo principal de cada ventana: 2/N, 4/N, 4/N y 6/N

**Qué se observó.** Medidos en el espectro de las ventanas, el ancho entre los primeros ceros del lóbulo principal fue 2.00/N (rectangular), 4.01/N (Hanning), 4.01/N (Hamming) y 6.01/N (Blackman). A −3 dB los anchos fueron 0.89/N, 1.45/N, 1.30/N y 1.65/N.

**Por qué sucede.** Multiplicar la señal por una ventana equivale, en frecuencia, a convolucionar su espectro con el espectro de la ventana. La rectangular da un núcleo de Dirichlet (parecido a un sinc), que vale cero en ω = 2πk/N para k distinto de cero; su primer cero está en 2π/N, así que su lóbulo principal mide, de nulo a nulo, 4π/N rad/muestra, es decir 2/N ciclos/muestra. Hanning y Hamming se forman sumando a la rectangular dos cosenos, y cada coseno en el tiempo desplaza una copia del espectro de la rectangular un bin hacia cada lado; la suma de tres núcleos desplazados forma un lóbulo principal que cubre más bins (4/N). Blackman agrega un segundo coseno, con copias a uno y dos bins, y llega a 6/N.

**Gracias a que.** Gracias a la propiedad de la transformada de Fourier de que un producto en el tiempo es una convolución en frecuencia, y a que multiplicar por un coseno desplaza el espectro. Es también la razón por la que el ancho es proporcional a 1/N: más muestras dan lóbulos más estrechos y mejor resolución.

**Error típico.** Pensar que el relleno de ceros (Nfft mayor que N) estrecha el lóbulo. No lo hace: solo interpola la curva.

---

## 2. El nivel de los lóbulos secundarios: −13.3, −31.5, −42.7 y −58.1 dB

**Qué se observó.** El mayor lóbulo lateral, relativo al pico principal, fue −13.3 dB en la rectangular, −31.5 dB en Hanning, −42.7 dB en Hamming y −58.1 dB en Blackman.

**Por qué sucede.** La causa física es el corte brusco de la señal en los bordes de la ventana. La rectangular corta de golpe y produce los laterales más altos. Las otras ventanas bajan suavemente hacia cero (Hanning llega a 0, Blackman llega a 0, Hamming queda en 0.08), y en frecuencia los núcleos desplazados que se suman tienen signos alternados, así que se cancelan parcialmente en los laterales.

**Gracias a que.** El −13.3 dB de la rectangular es el primer lóbulo lateral del núcleo de Dirichlet (cerca de 0.21 del pico). Los coeficientes de Hamming, 0.54 y 0.46 (aproximación de 25/46 y 21/46), están elegidos precisamente para cancelar casi del todo ese primer lóbulo lateral, y por eso baja a −43 dB. Hanning usa 0.5 y 0.5, que llevan la ventana a cero en los bordes.

**El compromiso inevitable.** Lo que se gana cancelando laterales se paga ensanchando el lóbulo principal. No existe una ventana mejor en todo.

---

## 3. El decaimiento de los lóbulos laterales: 6 dB/octava (rectangular y Hamming) frente a 18 dB/octava (Hanning y Blackman)

**Qué se observó.** En el espectro en dB, los laterales de la rectangular y de Hamming casi no bajan al alejarse del centro; los de Hanning y Blackman caen rápidamente, hasta menos de −100 dB.

**Por qué sucede.** La velocidad de caída depende de qué tan suave llega la ventana a sus bordes. La rectangular tiene un salto; Hamming conserva un pequeño salto (0.08) en los extremos; ambas decaen 6 dB por octava. Hanning y Blackman llegan a cero con derivada cero y decaen 18 dB por octava.

**Gracias a que.** Gracias a que la transformada de Fourier de una función con una discontinuidad decae como 1/ω (6 dB/octava), mientras que una función continua con derivada continua decae como 1/ω³ (18 dB/octava). Esto es el criterio de "velocidad de caída" que Harris (1978) usa junto con el ancho de banda equivalente de ruido y el nivel máximo de lóbulo lateral para comparar ventanas.

---

## 4. Experimento 1.a: con N = 1024 las cuatro ventanas separan las dos sinusoides

**Qué se observó.** Con A1 = A2 = 1, frecuencias 2π/12 y 2π/14 y ruido, se ven dos picos separados, en ω ≈ 0.449 y 0.524 rad/muestra, con las cuatro ventanas.

**Por qué sucede.** El periodograma ve cada línea espectral como una copia del lóbulo principal de la ventana. Dos picos se distinguen si su separación es mayor que el ancho de ese lóbulo. La separación es Δf = 1/12 − 1/14 = 0.0119 ciclos/muestra, que equivale a 0.0119 · 1024 = 12.2 bins. Esa distancia supera incluso al lóbulo más ancho (Blackman, 6 bins de nulo a nulo), así que todas las ventanas resuelven los dos picos.

**Con menos muestras sería distinto.** Con N = 128 la separación sería 1.5 bins. La rectangular (ancho a −3 dB de 0.89 bins) separaría los picos; Hanning (1.45 bins) quedaría al límite y los funde; Hamming y Blackman los fundirían. Esto se comprobó en una prueba con N = 128 hecha en Python.

**Gracias a que.** Gracias a que la resolución es proporcional a 1/N y a que la ventana convoluciona cada línea con su lóbulo principal.

---

## 5. La altura de los picos en el experimento 1.a

**Qué se observó.** El pico de cada sinusoide es más alto con la rectangular (≈ 80–90) y más bajo con Blackman (≈ 50–55).

**Por qué sucede.** Una ventana que atenúa más los bordes "ve" menos amplitud de la sinusoide, porque pondera con valores menores las muestras de los extremos.

**Gracias a que.** Para una sinusoide de amplitud A, con la normalización de la función periodogram (densidad unilateral en rad/muestra), el pico vale aproximadamente A²·N/(4π·ENBW). El ENBW (ancho de banda equivalente de ruido) se calcula como N·Σw²/(Σw)² y vale 1.00 (rectangular), 1.50 (Hanning), 1.36 (Hamming) y 1.73 (Blackman). Con A = 1 y N = 1024 resulta 81.5, 54.3, 59.8 y 47.1. Los valores medidos fluctúan alrededor de esos porque el ruido se suma al pico. La derivación: el módulo de la FFT cerca de ω₀ es aproximadamente (A/2)·Σw, de modo que |X|² ≈ A²(Σw)²/4, y la normalización 2/(2π·Σw²) da el resultado.

---

## 6. Experimento 1.b: la sinusoide débil no se ve con ninguna ventana

**Qué se observó.** Con A1 = 1 y A2 = 0.01, solo aparece el pico de la sinusoide fuerte, en escala lineal y en dB, con las cuatro ventanas.

**Por qué sucede.** La potencia de una sinusoide es A²/2. La razón de potencias entre la débil y la fuerte es (0.01/1)² = 10⁻⁴, que son −40 dB (10·log10(10⁻⁴) = −40). El pico de la fuerte ronda +19 dB y el de la débil debería quedar cerca de −21 dB (0.25·A2²·N/π ≈ 0.008 con la rectangular). El piso de ruido, en cambio, está en −7 dB. Es decir, el pico de la débil queda unos 14 dB por debajo del nivel medio del ruido, y el ruido (con máximos de +2 a +3.5 dB) lo cubre.

**Gracias a que.** Gracias a que la ventana solo controla cómo se reparte la fuga espectral; no puede elevar un pico por encima del piso de ruido. Lo único que cambiaría el resultado es la relación señal a ruido: el pico de una sinusoide crece proporcionalmente a N, mientras que el nivel de la DEP del ruido no cambia, así que con muchas más muestras la débil terminaría sobresaliendo.

**Error típico corregido.** Se había dicho que la débil estaba "80 dB por debajo". Es un error de cálculo: son 40 dB. Con amplitudes se usa 20·log10(0.01) = −40 dB, o con potencias 10·log10(10⁻⁴) = −40 dB.

---

## 7. Experimento 1.b sin ruido: Blackman y Hanning la ven, Hamming a medias y la rectangular no

**Qué se observó.** La sinusoide débil sobresale de la fuga de la fuerte 39.3 dB con Blackman, 31.5 dB con Hanning, 9.7 dB con Hamming y solo 0.6 dB con la rectangular.

**Por qué sucede.** Sin ruido, lo único que puede tapar a la débil es la fuga de la fuerte, que es el nivel de los lóbulos laterales de la ventana a la distancia de la débil: 12.2 bins. La rectangular decae solo 6 dB por octava partiendo de −13 dB; a 12.2 bins la envolvente del núcleo vale 1/(π·12.2) ≈ 0.026, es decir unos −32 dB respecto al pico de la fuerte, por encima de los −40 dB de la débil: la tapa. Hamming parte de un lateral muy bajo (−43 dB) pero también decae solo 6 dB/octava, así que a esa distancia su fuga queda unos 10 dB por debajo de la débil y se ve con poco margen. Hanning y Blackman decaen 18 dB/octava y a esa distancia la fuga cae por debajo de −60 dB, así que la débil sobresale 30 a 40 dB.

**Gracias a que.** Gracias a que la velocidad de decaimiento de los laterales depende de la suavidad de la ventana en los bordes (punto 3 de este documento).

**Cómo se midió.** Se compara el pico de la débil (cerca de w2) con el máximo de la fuga de la fuerte en esa misma zona, no con un solo punto, porque la fuga oscila entre picos y valles. Una primera medición comparaba con un solo punto y dio 12.8 dB para la rectangular, lo cual no coincidía con la figura; al comparar con el máximo da 0.6 dB.

---

## 8. El piso de ruido sale en −7 dB y no en −5 dB

**Qué se observó.** El nivel medio del ruido de varianza 1, lejos de las sinusoides, fue −7.1 a −7.3 dB.

**Por qué sucede.** La DEP verdadera del ruido blanco de varianza 1 vale σ² = 1 (bilateral). La función periodogram devuelve la densidad unilateral en rad/muestra: 2σ²/(2π) = 1/π ≈ 0.318, es decir −5.0 dB. Pero el periodograma del ruido, en cada frecuencia, no es una constante sino una variable aleatoria exponencial. El promedio de los valores en dB es el promedio de los logaritmos, y el promedio del logaritmo de una exponencial queda por debajo del logaritmo del promedio.

**Gracias a que.** Gracias a una propiedad de la distribución exponencial: E[ln X] = ln(media) − γ, con γ la constante de Euler (0.5772). En dB, 10·log10(e^(−γ)) = −2.51 dB. Entonces −5.0 − 2.5 ≈ −7.5 dB, y se midió −7.1 a −7.3 dB.

**Diferencias entre ventanas en el piso.** Los máximos del piso en la ventana rectangular (3.5 dB) frente a las otras (2.0 a 2.6 dB) difieren 1 a 1.5 dB, que es del orden de la fluctuación aleatoria; no son concluyentes. La ventana, al estar normalizada por Σw², no cambia el nivel del ruido.

---

## 9. Experimento 2: la curva es casi plana y su nivel medio es 0.32

**Qué se observó.** El periodograma promedio del ruido blanco queda aproximadamente horizontal alrededor de 0.32 (−5 dB), con las cuatro ventanas.

**Por qué es plana.** La autocorrelación del ruido blanco es σ²·δ[m]: un solo pico en m = 0. La transformada de un pico en el origen es una constante, así que la DEP es plana (todas las frecuencias tienen la misma potencia, de ahí el nombre "blanco").

**Por qué vale 0.32.** Es 2σ²/(2π) = 1/π = 0.318, la DEP unilateral en rad/muestra para σ² = 1. Las medias medidas fueron 0.320 a 0.323.

**Por qué las cuatro ventanas dan lo mismo.** Porque la DEP del ruido blanco no tiene estructura (no hay líneas que ensanchar ni componentes débiles que tapar) y porque la normalización por Σw² compensa la potencia que la ventana quita. La ventana solo importa para señales con estructura espectral.

**Gracias a que.** Gracias al teorema de Wiener-Khinchin (la DEP es la transformada de la autocorrelación) y a la normalización del periodograma modificado.

---

## 10. La varianza baja de ≈ 1 a ≈ 0.057 al promediar K = 16 segmentos

**Qué se observó.** La varianza relativa (varianza dividida por la media al cuadrado) de un solo periodograma fue 0.84 a 0.93; la del promedio de 16 segmentos fue 0.056 a 0.057, cerca de 1/16 = 0.0625.

**Por qué sucede.** La FFT de ruido gaussiano tiene parte real e imaginaria gaussianas e independientes en cada frecuencia, así que |X|² es una variable exponencial (chi-cuadrado con 2 grados de libertad), cuya desviación estándar es igual a su media: varianza relativa 1. Al promediar K estimaciones independientes, la varianza se divide entre K.

**Gracias a que.** Gracias a que la varianza de un promedio de K variables independientes es la varianza individual dividida entre K (el mismo argumento por el que el estimador de la media tiene varianza σ²/L). Una simulación de Monte Carlo de comprobación dio 0.9997 para un periodograma y 0.0616 para K = 16.

**Por qué el periodograma simple no es consistente.** Porque su varianza en cada frecuencia no baja al aumentar N; con más muestras se obtienen más puntos de frecuencia, no un promedio mayor por punto.

---

## 11. El costo de promediar: se pierde resolución

**Qué se observó.** Con N fijo y sin solape, K·L ≈ N, y más segmentos significa segmentos más cortos.

**Por qué sucede.** Cada segmento tiene solo L muestras, así que el ancho del lóbulo principal es proporcional a 1/L en vez de 1/N. Con N = 1024: L = 1024 da un ancho a −3 dB de la rectangular de 0.00087; L = 256 da 0.0035; L = 64 da 0.0139. Esta última supera la separación entre las sinusoides (0.0119), así que con L = 64 ya no se pueden separar.

**Gracias a que.** Gracias a que la resolución depende de la longitud de la ventana. La resolución limitada es un sesgo; la varianza es el otro componente del error cuadrático medio. Promediar cambia varianza por sesgo.

**Aclaración sobre una explicación incorrecta de una práctica anterior.** Una práctica de 2022 afirma que con L = 64 y K = 31 "la varianza aumenta". Es lo contrario: la varianza baja (hay más segmentos que promediar). Lo que empeora es la resolución.

---

## 12. Por qué en Octave hay que escribir 2*pi en pwelch

**Qué se observó.** Sin ese argumento, el eje de frecuencia de pwelch iba de 0 a 0.5 y el nivel del ruido era ≈ 2.0, en vez de eje 0 a π y nivel ≈ 0.32.

**Por qué sucede.** Cuando no se indica la frecuencia de muestreo Fs, pwelch en Octave usa Fs = 1: el eje queda en ciclos/muestra y la densidad en potencia por (ciclo/muestra), y para ruido de varianza 1 la densidad unilateral es 2σ²/Fs = 2. En cambio periodogram devuelve el eje en rad/muestra (de 0 a π) y la densidad por rad/muestra (nivel 1/π).

**Gracias a que.** Gracias a que la densidad espectral depende de la unidad de frecuencia: con Fs = 2π (rad/muestra) el nivel es 2σ²/(2π) = 0.318, igual al de periodogram. En MATLAB, con Fs omitida, pwelch ya usa rad/muestra.

---

## 13. Por qué max(blackman(N), 0) en Octave

**Qué se observó.** pwelch daba el error "window vector must be real and >=0".

**Por qué sucede.** blackman(N) en Octave devuelve −1.4·10⁻¹⁷ en los extremos por redondeo del punto flotante (la fórmula 0.42 − 0.5 + 0.08 da cero en teoría, pero no exactamente en la computadora), y pwelch rechaza ventanas con algún valor negativo.

**Solución.** max(blackman(N), 0) reemplaza esos valores minúsculos negativos por cero, sin cambiar la ventana en la práctica.

---

## 14. Por qué se normaliza por Σw² y no por Σw

**Por qué sucede.** La DEP es una potencia, cuadrática en la señal. Σw² mide la potencia que la ventana deja pasar de un ruido blanco. Sin esa normalización, una ventana que atenúa mucho los bordes (Blackman, Σw² ≈ 0.30·N) daría un nivel unas 3.3 veces menor que la rectangular (Σw² = N).

**Gracias a que.** Gracias a que para ruido blanco E[|X|²] = σ²·Σw². Dividir entre Σw² hace que el nivel de la DEP estimada no dependa de la ventana. Normalizar por Σw (la ganancia coherente) serviría para corregir la amplitud de los picos, no el nivel del ruido.

---

## 15. Por qué el relleno de ceros no mejora la resolución

**Por qué sucede.** La DFT de N puntos de una señal de N muestras solo muestrea la transformada de la señal enventanada en puntos separados 1/N ciclos/muestra. Si se agregan ceros hasta Nfft puntos, se muestrea la misma transformada en puntos más juntos: la curva se ve más suave y los picos se localizan mejor, pero la información es la misma.

**Gracias a que.** Gracias a que la transformada de la señal enventanada, V(e^jω), no cambia al agregar ceros; solo cambia cuántas muestras de ella se calculan. Su ancho de lóbulo lo fija la ventana de N puntos.

---

## 16. Por qué ninguna de las frecuencias cae justo en un bin

**Qué se observó.** Las sinusoides tienen fuga incluso con ventana rectangular.

**Por qué sucede.** Una sinusoide de frecuencia f = k₀/N con k₀ entero concentraría toda su energía en un solo punto de la DFT. Con N = 1024, f₁·N = 1024/12 = 85.33 y f₂·N = 1024/14 = 73.14: no son enteros, así que la energía se reparte entre varios bins.

**Gracias a que.** Gracias a que el núcleo de Dirichlet vale cero justo en los demás bins solo cuando la frecuencia cae exactamente en uno; si no, las muestras de la DFT caen en lóbulos del núcleo y no en sus ceros.

---

## Resumen en una frase por tema

- La ventana reduce la fuga espectral a cambio de resolución, porque el espectro de una ventana suave son núcleos desplazados que se cancelan en los laterales y se acumulan en el centro.
- La resolución es proporcional a 1/N, porque el lóbulo principal es proporcional a 1/N.
- El periodograma simple no es consistente, porque cada frecuencia es una variable exponencial con desviación igual a su media.
- El promedio reduce la varianza en 1/K, porque la varianza de un promedio de K estimaciones independientes es la individual dividida entre K, y cuesta resolución porque los segmentos son más cortos.
- Una sinusoide débil se tapa con el ruido o con la fuga: con ruido nada la salva salvo más muestras; sin ruido la salvan las ventanas cuyos laterales decaen 18 dB/octava.
