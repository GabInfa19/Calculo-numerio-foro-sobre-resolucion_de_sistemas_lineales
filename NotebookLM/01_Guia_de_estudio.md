# Guía de estudio: estimación espectral clásica (Procesamiento de Señales, Unidad I)


## 0. El mapa: de dónde viene todo

Toda la Unidad I responde a una sola pregunta: **si solo tengo un pedazo corto de una señal aleatoria, ¿cómo estimo en qué frecuencias está su potencia?** Cada tema de la clase es un peldaño para poder contestarla, y la práctica de las ventanas es el último peldaño.

Esta guía sigue el material de clase de tu Drive (*Estimadores de momentos estadísticos*, de María Cristi Stefanelli, UCAB) y lo completa con la parte investigativa de la práctica: las ventanas y los estimadores de Bartlett y Welch. Lo que viene de la clase está marcado como tal; lo que agrego yo (demostraciones, ejemplos, investigación) también.

La cadena de ideas, en orden, es esta:

1. **Una señal aleatoria no se puede escribir con una fórmula.** Se describe con un proceso aleatorio y con sus promedios (sección 1).
2. **Para que los promedios tengan sentido hace falta que no cambien con el tiempo** (estacionariedad) **y que un solo pedazo de señal represente a todo el proceso** (ergodicidad). Sin esto no se puede estimar nada (sección 2).
3. **La autocorrelación dice qué tan parecida es la señal consigo misma desplazada.** Su transformada de Fourier es la densidad espectral de potencia, la DEP: la respuesta a nuestra pregunta (sección 3).
4. **Pero nunca tenemos el proceso completo, solo L muestras.** Entonces la DEP se *estima*, y un estimador se juzga por su sesgo y su varianza (sección 4).
5. **Analizar L muestras con la DFT equivale a multiplicar por una ventana** y eso deforma el espectro (sección 5).
6. **El periodograma es el estimador más directo**, y tiene dos defectos: fuga y varianza alta (sección 6).
7. **La ventana ataca la fuga** (sección 7) **y promediar periodogramas ataca la varianza** (sección 8).
8. **La práctica prueba todo esto con simulaciones** (sección 9).

Si algo no se entiende, casi siempre es porque falta el peldaño anterior: vuelve uno atrás en la lista.

**Cómo usar la guía.** Cada sección tiene definiciones, fórmulas con su origen, un ejemplo y preguntas con respuesta. Intenta contestar cada pregunta *antes* de leer la respuesta: eso es lo que de verdad fija las ideas. Al final (secciones 10 y 11) hay preguntas de repaso, un glosario y los errores más comunes.

## 1. Procesos aleatorios y sus promedios

**Idea central:** una señal aleatoria no se puede predecir muestra a muestra, pero sus *promedios* sí son estables y describen su comportamiento.

### 1.1 Variable aleatoria (clase, tema 01)

Una **variable aleatoria (v.a.)** es una función que asigna un número a cada resultado posible de un experimento aleatorio. Ejemplo de la clase: lanzar una moneda; hay dos resultados (cara y sello) y se les asigna 0 y 1, así que X toma valores en {0, 1}. Es **discreta** si toma un número finito o numerable de valores (moneda) y **continua** si toma valores en un intervalo (la temperatura de una persona elegida al azar).

Una v.a. se describe con dos funciones:

```latex
F(x) = P[X \le x] \qquad\qquad p(x) = \frac{dF(x)}{dx}
```

- **F(x)**, la función de distribución acumulada (FDA): la probabilidad de que X sea menor o igual que x. Siempre va de 0 a 1 y nunca decrece.
- **p(x)**, la función de densidad de probabilidad (fdp): la derivada de F. El área bajo p(x) entre x1 y x2 es la probabilidad de que X caiga en ese intervalo.

### 1.2 Proceso aleatorio (clase, tema 02)

Una señal discreta aleatoria se modela como un **proceso estocástico**: una familia de v.a., una por cada instante n.

```latex
X[n] = \{X_n\}
```

Piensa en 1000 micrófonos grabando a la vez el mismo ruido de fondo. Cada grabación es una **realización** (una señal concreta, x₃\[n\], x₇\[n\], …). El conjunto de todas las realizaciones posibles es el **ensamble**. Si fijas un instante n = k y miras los 1000 valores de ese instante, tienes una v.a. X\[k\].

Como manejar las funciones F y p de cada instante es incómodo, se usan **promedios estadísticos**.

### 1.3 Promedios estadísticos (de ensamble)

| Promedio | Fórmula | Qué mide |
| --- | --- | --- |
| Media (valor esperado) | E\[X\_n\] = ∫ x p(x, n) dx | El valor central; en una señal eléctrica, el nivel DC |
| Valor cuadrático medio | E\[X\_n²\] = ∫ x² p(x, n) dx | Potencia promedio total |
| Varianza | Var\[X\_n\] = E\[\|X\_n − E\[X\_n\]\|²\] = σ² | Potencia de las variaciones alrededor de la media (potencia AC) |
| Autocorrelación | φxx\[n, m\] = E\[X\_n X\_m\*\] | Qué tan relacionado está el valor en n con el valor en m |
| Autocovarianza | γxx\[n, m\] = E\[(X\_n − m\_n)(X\_m − m\_m)\*\] | Igual, pero quitando antes la media |

**De dónde sale la media:** es el promedio de todos los valores posibles, cada uno pesado por su probabilidad (como una nota ponderada). La **varianza** nace de preguntar “¿cuánto se aleja de la media, en promedio?”; se eleva al cuadrado para que las desviaciones positivas y negativas no se cancelen. Su raíz, σ, es la **desviación estándar** (el valor RMS de la parte AC).

Propiedades que se usan siempre (clase):

- E\[X + Y\] = E\[X\] + E\[Y\], E\[k\] = k, E\[kX\] = k·E\[X\].
- Si X e Y son independientes: E\[XY\] = E\[X\]·E\[Y\] y Var\[X + Y\] = Var\[X\] + Var\[Y\].
- Var\[kX\] = k²·Var\[X\] y Var\[constante\] = 0.
- Se relacionan así: autocorrelación = autocovarianza + (media)². Si la media es cero, **coinciden**.

### 1.4 Promedios temporales

En la práctica no se tienen 1000 micrófonos, sino una sola realización. Sobre ella se calculan promedios *en el tiempo*:

```latex
\langle x \rangle = \lim_{L\to\infty} \frac{1}{2L+1}\sum_{n=-L}^{L} x[n] \qquad\qquad \langle x[n+m]\,x^*[n] \rangle = \lim_{L\to\infty} \frac{1}{2L+1}\sum_{n=-L}^{L} x[n+m]\,x^*[n]
```

La pregunta clave de la siguiente sección es: ¿cuándo el promedio en el tiempo de *una* señal coincide con el promedio del ensamble?

### 1.5 Ejemplo resuelto: ruido blanco

El ruido blanco de la práctica (randn) cumple que cada muestra es independiente de las demás, con media 0 y varianza σ² = 1. Su autocorrelación sale directamente de las propiedades de arriba:

- Si n ≠ m: E\[X\_n X\_m\] = E\[X\_n\]·E\[X\_m\] = 0·0 = 0 (independientes).
- Si n = m: E\[X\_n²\] = Var\[X\_n\] + (media)² = σ² + 0 = σ².

```latex
\varphi_{xx}[m] = \sigma^2\,\delta[m]
```

Es un “pico” en m = 0 y cero en todo lo demás: el ruido blanco no se parece a sí mismo desplazado, ni siquiera una muestra.

### Preguntas

**P1. ¿Qué diferencia hay entre una realización y un proceso?** R: el proceso es la familia completa de señales posibles (el ensamble); una realización es una señal concreta de esa familia. El código de la práctica genera una realización de ruido cada vez que llama a randn.

**P2. ¿Por qué la autocovarianza y la autocorrelación coinciden cuando la media es cero?** R: porque la autocovarianza es la autocorrelación menos el producto de las medias; si la media es 0, ese producto vale 0.

**P3. ¿Qué es σ² físicamente en una señal eléctrica?** R: la potencia promedio de la parte AC (lo que oscila alrededor del nivel DC). Por eso el ruido de varianza 1 “tiene potencia 1”.

## 2. Estacionariedad, ergodicidad y autocorrelación

**Idea central:** para estimar algo a partir de *una* señal, ese algo tiene que ser el mismo en todos los instantes (estacionariedad) y la señal tiene que ser “típica” de todo el proceso (ergodicidad).

### 2.1 Estacionariedad (clase, tema 03)

- **En sentido estricto:** todas las propiedades estadísticas son independientes del tiempo. La densidad conjunta no cambia si se desplazan todos los instantes en τ: p(x₁,…,x\_k; n₁,…,n\_k) = p(x₁,…,x\_k; n₁+τ,…,n\_k+τ).
- **En sentido amplio (WSS):** solo se pide que los promedios de primer orden no dependan del tiempo y que los de segundo orden dependan únicamente de la *diferencia* de tiempos. Es la que realmente se usa en la materia.

Un proceso WSS cumple tres cosas:

```latex
E[X_n] = m_x \ (\text{constante}) \qquad \text{Var}[X_n] = \sigma_x^2 \ (\text{constante}) \qquad \varphi_{xx}[n+m,\,n] = E[X_{n+m}X_n^*] = \varphi_{xx}[m]
```

La última línea es importante: la autocorrelación pasa de depender de dos variables (n y m) a depender solo del **desfase m**. Es una secuencia unidimensional, mucho más manejable.

**Ejemplos.** El ruido blanco de la práctica es estacionario. Una señal de voz *no* lo es si se la mira durante segundos (cambia de vocal a vocal), pero sí es aproximadamente estacionaria en tramos de 5 a 100 ms; por eso en las prácticas de voz se cortaba en tramos y se calculaba un espectro por tramo.

### 2.2 Propiedades de la autocorrelación de un proceso WSS (clase)

| Propiedad | Fórmula | Qué significa |
| --- | --- | --- |
| Es par | φxx\[m\] = φxx\*\[−m\] (real: φxx\[m\] = φxx\[−m\]) | Comparar x\[n\] con x\[n+m\] es lo mismo que comparar con x\[n−m\] |
| Máximo en el origen | \|φxx\[m\]\| ≤ φxx\[0\] | Nada se parece más a la señal que ella misma sin desplazar |
| Valor en cero | φxx\[0\] = E\[\|X\|²\] | Valor cuadrático medio = **potencia promedio total** |
| Autocovarianza en cero | γxx\[0\] = Var\[X\] | Potencia de la parte AC |
| Decorrelación | γxx\[m\] → 0 cuando m → ∞ (en muchos procesos) | Las muestras muy separadas dejan de estar relacionadas |

La autocorrelación y la autocovarianza de un proceso WSS son secuencias **determinísticas** (no aleatorias) y de energía finita, así que tienen transformada de Fourier y Z. Este es el puente hacia la DEP en la sección 3.

### 2.3 Ergodicidad (clase)

Un proceso es **ergódico** si cumple dos condiciones:

1. Sus promedios en el tiempo, ⟨X\_n⟩ y ⟨X\_{n+m}X\_n\*⟩, **no dependen de qué realización** se escoja.
2. Esos promedios en el tiempo **coinciden** con los del ensamble, m\_x y φxx\[m\].

De aquí se concluye (clase): **todo proceso ergódico es WSS, pero no todo WSS es ergódico.**

**Por qué importa tanto.** Gracias a la ergodicidad se pueden medir los promedios del ensamble (que no tenemos) con promedios en el tiempo de una sola señal (que sí tenemos). Con una señal ergódica, las cantidades eléctricas de siempre coinciden con las estadísticas (clase):

| Cantidad eléctrica | Expresión |
| --- | --- |
| Valor DC | E\[x\] = m\_x |
| Potencia promedio total | φxx\[0\] = E\[x²\] = ⟨x²⟩ |
| Potencia promedio DC | (E\[x\])² |
| Potencia promedio AC | Var(x) = E\[x²\] − (E\[x\])² |
| Valor RMS de la parte AC | σx = √Var(x) |

### 2.4 Dos contraejemplos que aclaran todo

- **Estacionario pero no ergódico.** Sea X\[n\] = A para todo n, donde A es una constante *aleatoria* con media 0 y varianza σ² que se sortea una sola vez por realización. Entonces E\[X\_n\] = 0 y φxx\[m\] = σ² (no depende de n, así que es WSS). Pero el promedio en el tiempo de una realización es A, que cambia de realización en realización y no vale 0: no es ergódico. Con una sola señal jamás podrías descubrir que la media del ensamble es 0.
- **Estacionario y ergódico.** El ruido blanco: cualquier realización larga tiene media temporal cercana a 0 y potencia temporal cercana a σ². Por eso en la práctica basta una sola realización de ruido para estimar su espectro.

### Preguntas

**P1. ¿Por qué se necesita estacionariedad para usar un promedio en el tiempo?** R: si las propiedades cambian con n, el promedio temporal mezcla comportamientos distintos y no estima ningún parámetro fijo. Estacionariedad significa “el parámetro que quiero estimar es el mismo en todo instante”.

**P2. ¿Por qué φxx\[0\] es la potencia promedio?** R: porque φxx\[0\] = E\[X\_n·X\_n\*\] = E\[|X\_n|²\], y el valor cuadrático medio es, por definición, la potencia promedio de la señal.

**P3. ¿El ruido blanco gaussiano es ergódico?** R: sí. Es WSS y sus muestras se decorrelacionan por completo (φxx\[m\] = 0 para m ≠ 0), así que los promedios temporales convergen a los del ensamble.

**P4. ¿Qué pasa si mi señal de voz dura 3 segundos?** R: no es estacionaria en todo ese intervalo. Hay que cortarla en tramos de 5 a 100 ms, donde se comporta como estacionaria, y estimar un espectro por tramo.

## 3. La densidad espectral de potencia (DEP)

**Idea central:** la DEP reparte la potencia de la señal entre las frecuencias. Es lo que queremos estimar en toda la unidad.

### 3.1 Definición y de dónde sale (clase, tema 04)

Como la autocorrelación de un proceso WSS es una secuencia determinística de energía finita, tiene transformada de Fourier. A esa transformada se le llama **densidad espectral de potencia** (teorema de Wiener-Khinchin):

```latex
\Phi_{xx}(e^{j\omega}) = \sum_{m=-\infty}^{\infty} \varphi_{xx}[m]\,e^{-j\omega m} \qquad\qquad \varphi_{xx}[m] = \frac{1}{2\pi}\int_{-\pi}^{\pi}\Phi_{xx}(e^{j\omega})\,e^{j\omega m}\,d\omega
```

**Por qué se llama así.** Si evaluamos la fórmula inversa en m = 0, el exponencial vale 1 y queda:

```latex
\varphi_{xx}[0] = E[|x|^2] = \frac{1}{2\pi}\int_{-\pi}^{\pi}\Phi_{xx}(e^{j\omega})\,d\omega
```

O sea: el **área** bajo la DEP es la potencia promedio total. Por eso es una *densidad de potencia* (al integrarla da potencia) y es *espectral* porque es función de la frecuencia. Para un proceso ergódico, esa potencia coincide con la calculada sobre una señal cualquiera del proceso (clase).

Si el proceso es real, la DEP es **real, par y no negativa**. Si la media es cero, la DEP de la autocorrelación y la de la autocovarianza son la misma.

### 3.2 Unidades y dos formas de medir la frecuencia

- En tiempo continuo, la DEP se mide en W/Hz. En tiempo discreto es potencia por rad/muestra (o por ciclo/muestra).
- La frecuencia angular ω (rad/muestra) y la frecuencia normalizada f (ciclos/muestra) se relacionan por ω = 2πf. Las sinusoides de la práctica, ω₁ = 2π/12 y ω₂ = 2π/14, son f₁ = 1/12 ≈ 0.0833 y f₂ = 1/14 ≈ 0.0714 ciclos/muestra, o ω₁ ≈ 0.524 y ω₂ ≈ 0.449 rad/muestra.
- Para señales reales basta mirar de 0 a π (la parte negativa es un espejo). Octave/MATLAB devuelven la versión **unilateral** y por eso duplican el nivel.

### 3.3 Ejemplos resueltos (las tres DEP de la práctica)

**Ejemplo 1. Ruido blanco.** Ya vimos que φxx\[m\] = σ²·δ\[m\]. Su transformada es una suma con un solo término distinto de cero:

```latex
\Phi_{xx}(e^{j\omega}) = \sum_m \sigma^2\,\delta[m]\,e^{-j\omega m} = \sigma^2
```

La DEP es **constante**: todas las frecuencias tienen la misma potencia. De ahí el nombre “blanco”, por analogía con la luz blanca (todos los colores por igual). Su área es (1/2π)·σ²·2π = σ², la potencia, como debe ser.

**Ejemplo 2. Una sinusoide.** Sea x\[n\] = A·cos(ω₀n + θ), con la fase θ uniforme en \[0, 2π) (así el proceso es estacionario). Su autocorrelación es:

```latex
\varphi_{xx}[m] = E\left[A\cos(\omega_0 (n+m)+\theta)\,A\cos(\omega_0 n+\theta)\right] = \frac{A^2}{2}\cos(\omega_0 m)
```

Sale de la identidad cos(a)·cos(b) = ½\[cos(a−b) + cos(a+b)\]: el primer término da cos(ω₀m), que no depende de θ, y el segundo promedia 0 sobre θ. El coseno tiene como transformada dos deltas:

```latex
\Phi_{xx}(e^{j\omega}) = \frac{\pi A^2}{2}\left[\delta(\omega-\omega_0)+\delta(\omega+\omega_0)\right]
```

Comprobación: la potencia es (1/2π)·(πA²/2)·2 = A²/2, que es la potencia conocida de una sinusoide. Toda su potencia está concentrada en una sola frecuencia: **una línea espectral**.

**Ejemplo 3. Sinusoides más ruido (la señal de la práctica).** Si las sinusoides y el ruido son independientes, las DEP se suman: líneas espectrales en ω₁ y ω₂ encima de un piso constante σ² (por eso la DEP de la suma es la suma de las DEP “solo si no están correlacionados”, como dice la clase). Un buen estimador debe mostrar dos picos sobre un piso plano. La potencia de cada línea es A²/2:

| Componente | Potencia | Forma en la DEP |
| --- | --- | --- |
| Sinusoide de amplitud 1 | 1/2 = 0.5 | Línea (delta) en ω₀ |
| Sinusoide de amplitud 0.01 | 0.0001/2 = 0.00005 | Línea (delta) mucho más pequeña |
| Ruido de varianza 1 | 1 | Piso plano de valor σ² = 1 |

### Preguntas

**P1. ¿Por qué la DEP del ruido blanco es plana?** R: porque su autocorrelación es un solo pico en m = 0, y la transformada de un pico en el origen es una constante.

**P2. ¿Qué representa el área bajo la DEP?** R: la potencia promedio total de la señal, φxx\[0\] (dividida por 2π en la convención de rad/muestra).

**P3. ¿Por qué la DEP no puede ser negativa?** R: porque es una potencia por unidad de frecuencia, y una potencia no puede ser negativa. Matemáticamente, Φxx es el límite de |X(ω)|²/N, que es un cuadrado.

**P4. ¿Cuántas veces es más pequeña la potencia de la sinusoide de amplitud 0.01 que la de amplitud 1?** R: la potencia depende de A², así que la razón es (0.01/1)² = 10⁻⁴, es decir 10 000 veces menos, o 40 dB menos (10·log₁₀(10⁻⁴) = −40 dB). Este número es la clave del experimento 1.b.

## 4. Estimar con datos finitos: sesgo, varianza y consistencia

**Idea central:** nunca tenemos el proceso completo ni infinitas muestras, así que todo lo que calculamos es un *estimador*, y un estimador es él mismo una variable aleatoria. Hay que saber juzgar si es bueno.

### 4.1 Qué es un estimador (clase, “Estimación de parámetros”)

Estimar un parámetro es preguntarse cuánta información del parámetro tienen las muestras. Como no se conoce el valor exacto, se habla de la *probabilidad* de que el valor estimado coincida con el verdadero. El estimado es una v.a. que, ojalá, tiene dos virtudes:

- que **en promedio** acierte (sin sesgo);
- que su **variación** sea pequeña y tienda a cero al aumentar las muestras (consistencia).

Un criterio para escoger estimador es el **error cuadrático medio (ECM)** mínimo:

```latex
\text{ECM} = E\left[(\hat{\theta}-\theta)^2\right] = \underbrace{\left(E[\hat\theta]-\theta\right)^2}_{\text{sesgo}^2} + \underbrace{\text{Var}[\hat\theta]}_{\text{varianza}}
```

Esta descomposición (que sale de sumar y restar E\[θ̂\] dentro del cuadrado) dice que el error tiene dos fuentes independientes: apuntar mal (sesgo) y apuntar bien pero con mala puntería (varianza). Imagina lanzar dardos: sesgo es que el grupo de dardos esté corrido del centro; varianza es que estén dispersos.

### 4.2 Las tres clasificaciones (clase)

| Tipo | Condición | En palabras |
| --- | --- | --- |
| **No sesgado** | E\[θ̂\] = θ | En promedio acierta |
| **Asintóticamente no sesgado** | E\[θ̂\] → θ cuando L → ∞ | Puede errar un poco con pocas muestras, pero el error desaparece al crecer L |
| **Consistente** | Var\[θ̂\] → 0 cuando L → ∞ | Con muchas muestras se concentra en un valor |

### 4.3 Dos estimadores concretos (clase)

```latex
\hat{m}_x = \frac{1}{L}\sum_{n=0}^{L-1}x[n] \qquad\qquad \hat{\sigma}_x^2 = \frac{1}{L}\sum_{n=0}^{L-1}\left(x[n]-\hat{m}_x\right)^2
```

La clase afirma que el estimador de la media es **no sesgado y consistente**, y el de la varianza es **asintóticamente no sesgado y consistente**.

**Demostración de que la media es no sesgada (clase, con todos los pasos):**

```latex
E[\hat m_x] = E\left[\frac{1}{L}\sum_{n=0}^{L-1}x[n]\right] = \frac{1}{L}\sum_{n=0}^{L-1}E[x[n]] = \frac{1}{L}\cdot L\cdot E[x] = E[x]
```

Se usan dos propiedades: el valor esperado de una suma es la suma de los valores esperados, y como el proceso es estacionario E\[x\[n\]\] no depende de n y sale de la sumatoria.

**Y es consistente (complemento mío).** Si las muestras son independientes con varianza σ², la varianza de la suma es la suma de varianzas, y al dividir por L el factor se eleva al cuadrado:

```latex
\text{Var}[\hat m_x] = \frac{1}{L^2}\sum_{n=0}^{L-1}\sigma^2 = \frac{\sigma^2}{L}\ \xrightarrow[L\to\infty]{}\ 0
```

**Ejemplo numérico.** Con ruido blanco de σ² = 1 y L = 1024 muestras, la media estimada tiene desviación estándar √(1/1024) = 1/32 ≈ 0.031. O sea, casi siempre sale entre −0.06 y 0.06 aunque la media real sea 0. Con L = 128 la desviación sería 0.088. Cuantas más muestras, más fiable.

**La varianza con 1/L en vez de 1/(L−1).** Con 1/L el estimador tiene un sesgo pequeño, E\[σ̂²\] = (L−1)/L·σ², que desaparece cuando L crece: por eso es “asintóticamente no sesgado”, como dice la clase. (Con 1/(L−1) es exactamente no sesgado.)

### 4.4 Por qué esto es la base de la práctica

La DEP también se estima, y el periodograma se juzga con esta misma vara: ¿tiene sesgo? ¿tiene varianza que se anule al crecer L? La respuesta, que vemos en la sección 6, es que el periodograma es (asintóticamente) no sesgado pero **no consistente**: su varianza no baja aunque aumentes las muestras.

> **En nuestra tarea.** Medimos esas dos propiedades con Octave sobre ruido blanco de varianza 1 (experimento 2). La **media** de la estimación fue 0.320 a 0.323 en las cuatro ventanas, igual al valor verdadero 1/π ≈ 0.318: sin sesgo apreciable. La **varianza relativa** (varianza dividida por media al cuadrado, que no depende de la escala) de un solo periodograma fue 0.84 a 0.93, es decir ≈ 1; al promediar K = 16 segmentos bajó a 0.056 – 0.057, cerca de 1/16 = 0.0625.

### Preguntas

**P1. ¿Qué es mejor: un estimador sin sesgo y con mucha varianza, o uno con un poco de sesgo y poca varianza?** R: depende del ECM, que suma sesgo² y varianza. A veces conviene aceptar un sesgo pequeño si baja mucho la varianza. Eso es justo lo que hace el periodograma promedio: pierde un poco de resolución (sesgo) para ganar una varianza mucho menor.

**P2. ¿Si un estimador es consistente, también es no sesgado?** R: no necesariamente. Puede tener sesgo que se vaya a cero (asintóticamente no sesgado) o incluso sesgo permanente con varianza que se anula. Son propiedades distintas; la clase las lista por separado.

**P3. ¿Por qué la varianza de la media baja como 1/L?** R: porque sumar L muestras independientes multiplica la varianza por L, y dividir entre L la divide por L²; el cociente neto es 1/L.

## 5. Del análisis con DFT a las ventanas

**Idea central:** para calcular un espectro en el computador hay que usar una cantidad finita de muestras, y eso equivale a multiplicar la señal por una ventana. Esa multiplicación tiene un efecto preciso y calculable en frecuencia.

### 5.1 Antes: sistemas lineales con entrada aleatoria (clase, tema 05)

Un sistema lineal e invariante (LTI) tiene respuesta al impulso h\[n\] y respuesta en frecuencia H(e^jω). Si la entrada es aleatoria, la salida también, y no se puede describir con fórmulas; se describe con la autocorrelación y la DEP, que sí son determinísticas:

```latex
\varphi_{yy}[m] = \varphi_{xx}[m] * h[m] * h^*[-m] \qquad\Longleftrightarrow\qquad \Phi_{yy}(e^{j\omega}) = |H(e^{j\omega})|^2\,\Phi_{xx}(e^{j\omega})
```

```latex
\varphi_{yx}[m] = h[m]*\varphi_{xx}[m] \qquad\Longleftrightarrow\qquad \Phi_{yx}(e^{j\omega}) = H(e^{j\omega})\,\Phi_{xx}(e^{j\omega})
```

**Ejemplo.** Si se filtra ruido blanco (Φxx = σ²) con un filtro, la salida tiene DEP σ²|H(e^jω)|²: el filtro “colorea” el ruido dando forma a su espectro. Esta idea es la base de los modelos AR de la siguiente unidad. En nuestra práctica no se filtra nada, pero explica por qué la DEP es la herramienta natural para señales aleatorias: pasa limpiamente por los sistemas.

### 5.2 El proceso de análisis de Fourier de una señal (clase, figuras 4 y 5, Oppenheim 2011)

La clase presenta los pasos para obtener el espectro de una señal continua con la DFT:

1. **Filtro pasa bajo antisolapamiento.** Limita el ancho de banda para que, al muestrear, las frecuencias altas no se “doblen” sobre las bajas (solapamiento). Como el filtro real no es ideal, queda algo de solapamiento, pero se hace despreciable con un buen filtro.
2. **Muestreo (conversión A/D).** Se pasa a la secuencia x\[n\]. Debe cumplirse el teorema de Nyquist-Shannon: muestrear al menos al doble de la frecuencia máxima. La cuantificación introduce un ruido que se modela como una secuencia sumada a x\[n\] y se hace despreciable con cuantificación fina.
3. **Multiplicación por una ventana w\[n\].** La DFT exige longitud finita, así que se toma un tramo: v\[n\] = w\[n\]·x\[n\].
4. **DFT de v\[n\].** Da V\[k\], muestras equiespaciadas de la transformada de Fourier de v\[n\].

**En nuestra tarea.** Los pasos 1 y 2 no existen: la señal la creamos directamente en el computador, ya muestreada (sinusoides con sin y ruido con randn). Los pasos 3 y 4 son exactamente lo que hace periodogram(x, ventana, Nfft).

### 5.3 El efecto de la ventana: multiplicar en el tiempo es convolucionar en frecuencia

La clase lo dice así: el efecto de multiplicar por una ventana de longitud finita es, en frecuencia, una **convolución periódica** entre el espectro de x\[n\] y el espectro de w\[n\].

```latex
v[n] = w[n]\,x[n] \qquad\Longleftrightarrow\qquad V(e^{j\omega}) = \frac{1}{2\pi}\int_{-\pi}^{\pi}X(e^{j\theta})\,W\!\left(e^{j(\omega-\theta)}\right)d\theta
```

La clase describe el espectro de una ventana típica: un **lóbulo principal** concentrado alrededor de ω = 0 y **lóbulos laterales** pequeños. La convolución con W reparte cada componente de X como una copia de W: una línea fina como una delta se convierte en una copia de W desplazada a esa frecuencia. Dos consecuencias:

- Cada línea espectral queda **ensanchada** al ancho del lóbulo principal (la resolución se limita).
- Cada línea “derrama” energía hacia los lados por los lóbulos laterales (la **fuga espectral**).

**Lo que sale para la ventana rectangular** (la que usa un periodograma sin ventana, solo por tomar N muestras):

```latex
W_R(e^{j\omega}) = \sum_{n=0}^{N-1}e^{-j\omega n} = e^{-j\omega(N-1)/2}\,\frac{\sin(\omega N/2)}{\sin(\omega/2)}
```

Es el núcleo de Dirichlet. Se anula en ω = 2πk/N para k ≠ 0. El primer cero está en ω = 2π/N, así que el lóbulo principal mide de nulo a nulo 4π/N rad/muestra = **2/N ciclos/muestra**. Esto es de donde sale el “2/N” de la tabla de ventanas.

### 5.4 La DFT: muestras de la transformada

```latex
X[k] = \sum_{n=0}^{N-1}x[n]\,e^{-j2\pi kn/N},\qquad k=0,\ldots,N-1
```

Según la clase, V\[k\] = V(e^jω) con ω = 2πk/N: la DFT solo *muestrea* la transformada de v\[n\] en puntos separados 2π/N rad/muestra (1/N ciclos/muestra). Si la ventana tiene longitud L ≤ N, la DFT de N puntos agrega N − L ceros (**zero-padding**): la curva se ve más suave porque hay más puntos, pero la información no aumenta ni mejora la resolución; solo la interpola.

**Ejemplo clave: ¿la frecuencia cae justo en un bin?** Una sinusoide de frecuencia f = k₀/N con k₀ entero cae exactamente en un punto de la DFT y ahí concentra toda su energía (los demás puntos caen en ceros del núcleo de Dirichlet). Si no cae en un bin, la energía se reparte entre varios.

> **En nuestra tarea.** Con N = 1024: f₁ · N = 1024/12 = 85.33 y f₂ · N = 1024/14 = 73.14. Ninguno es entero, así que ambas sinusoides tienen fuga *incluso con ventana rectangular*. Con N = 128 pasaría lo mismo (10.67 y 9.14). Por eso la ventana importa en el experimento 1.

### Preguntas

**P1. ¿Por qué hace falta una ventana para calcular la DFT?** R: porque la DFT trabaja con un número finito de muestras. Tomar N muestras ya es multiplicar por una ventana rectangular; la cuestión es solo cuál ventana usar.

**P2. ¿Qué hace el zero-padding y qué no hace?** R: agrega puntos entre los bins originales, así que la curva se ve más suave y los picos se localizan mejor. No agrega información ni separa picos que el lóbulo principal ya fundió.

**P3. ¿Por qué una delta se ve ancha en el periodograma?** R: porque el periodograma ve la delta convolucionada con el espectro de la ventana; en lugar de una línea infinitamente fina aparece una copia del lóbulo principal.

**P4. ¿Cuántos puntos de la DFT separa un par de frecuencias f₁ = 1/12 y f₂ = 1/14 con N = 1024?** R: Δf · N = (1/12 − 1/14) · 1024 = 12.2 bins. Con N = 128 serían solo 1.5 bins.

## 6. El periodograma

**Idea central:** el periodograma es el estimador más directo de la DEP: calcular la FFT de los datos y elevar al cuadrado. Es rápido y natural, pero tiene dos defectos que motivan el resto de la unidad.

### 6.1 Definición

```latex
\hat{P}_{per}(\omega) = \frac{1}{N}\left|\sum_{n=0}^{N-1}x[n]\,e^{-j\omega n}\right|^2 = \frac{1}{N}|X(e^{j\omega})|^2
```

En el computador se evalua en los puntos de la FFT: P\[k\] = |X\[k\]|²/N, con X\[k\] la FFT de las N muestras (sección 5.4). Ese es el estimador que calcula `periodogram` cuando no se le da ventana.

### 6.2 De dónde sale: es la transformada de una autocorrelación estimada

La DEP es la transformada de la autocorrelación (sección 3). Con N muestras no podemos promediar infinito, así que estimamos la autocorrelación con un promedio temporal finito:

```latex
\hat{r}[m] = \frac{1}{N}\sum_{n=0}^{N-1-|m|}x[n+|m|]\,x^*[n], \qquad |m| \le N-1
```

Se divide entre N (no entre N − |m|), lo que da un estimador con sesgo pero que garantiza que su transformada nunca sea negativa. Ahora el paso clave: expandir |X|² como doble suma y agrupar los términos con el mismo desfase m = n − k:

```latex
|X(e^{j\omega})|^2 = \sum_{n}\sum_{k}x[n]x^*[k]\,e^{-j\omega(n-k)} = \sum_{m=-(N-1)}^{N-1}\Big(\sum_k x[k+m]x^*[k]\Big)e^{-j\omega m} = N\sum_{m}\hat r[m]\,e^{-j\omega m}
```

Es decir: **periodograma = transformada de Fourier de la autocorrelación estimada**. Es la definición de DEP con r̂ en lugar de φ. (Por eso en las prácticas de cursos anteriores se comparaba el periodograma con “un estimador usando la FFT de la autocorrelación”: son lo mismo si se usa la versión con 1/N.)

### 6.3 Media: casi sin sesgo, pero con fuga

Como E\[r̂\[m\]\] = (1 − |m|/N)·φ\[m\], la autocorrelación estimada es la verdadera multiplicada por una ventana triangular (Bartlett). En frecuencia, el valor esperado del periodograma es la DEP verdadera convolucionada con el espectro de esa ventana:

```latex
E\left[\hat P_{per}(\omega)\right] = \frac{1}{2\pi}\int_{-\pi}^{\pi}\Phi_{xx}(\theta)\,W_B(\omega-\theta)\,d\theta
```

Cuando N crece, W\_B se vuelve cada vez más estrecho y se parece a una delta, así que E\[P\] → Φ: es **asintóticamente no sesgado**. Con N finito, el sesgo es justamente el ensanchamiento y la fuga de la sección 5.3.

### 6.4 Varianza: no baja al aumentar N

Para un proceso gaussiano, la FFT en cada frecuencia es un número complejo con partes real e imaginaria gaussianas e independientes. Su módulo cuadrado, |X|², es entonces una variable **exponencial** (chi-cuadrado con 2 grados de libertad), y una exponencial tiene desviación estándar igual a su media:

```latex
\text{Var}\left[\hat P_{per}(\omega)\right] \approx \Phi_{xx}^2(\omega) \qquad\Longrightarrow\qquad \frac{\text{Var}}{(\text{media})^2} \approx 1
```

Lo grave: esto es cierto **sin importar N**. Si aumentas las muestras no baja la varianza en cada frecuencia; lo que aumenta es el número de frecuencias que calculas (la curva se vuelve más fina y más “peluda”). Se necesita información *nueva* por cada punto de frecuencia y no la hay. Por definición de la sección 4, un estimador cuya varianza no tiende a 0 **no es consistente**.

### 6.5 Resolución

La resolución es el menor Δf que permite ver dos picos separados. Para la ventana rectangular es del orden de 1/N ciclos/muestra (ancho a −3 dB ≈ 0.89/N). Con N = 128 es ≈ 0.0070; con N = 1024, ≈ 0.00087. Como Δf entre nuestras sinusoides es 0.0119, N = 1024 las separa con mucho margen (sección 5.4, P4).

### 6.6 El periodograma modificado (con ventana)

Si antes de la FFT se multiplica por una ventana w\[n\], se obtiene el **periodograma modificado**:

```latex
\hat P_M(\omega) = \frac{1}{N\,U}\left|\sum_{n=0}^{N-1}w[n]\,x[n]\,e^{-j\omega n}\right|^2, \qquad U = \frac{1}{N}\sum_{n=0}^{N-1}w^2[n]
```

La división por U compensa la potencia que la ventana le quita a la señal, para que el *nivel* de la DEP estimada no dependa de la ventana. Sin ella, una ventana como Blackman (que atenua mucho los bordes) daría niveles más bajos que la rectangular. **El problema de la fuga se ataca aquí; el de la varianza no** (la ventana no cambia que cada punto sea una exponencial).

**Detalle de escala en Octave/MATLAB.** La función `periodogram` devuelve densidad por rad/muestra y *unilateral*: Pxx = 2·|X|² / (2π·Σw²) para 0 ≤ ω ≤ π. Para ruido blanco de varianza 1 eso da σ²/π ≈ **0.318** (−5 dB), no 1. Ese es el nivel que verás en las gráficas de la práctica (no es un error: el 2 viene de ser unilateral y el 1/2π de la convención de rad/muestra).

> **En nuestra tarea.** `[Pxx, w] = periodogram(xa, ventanas{i}, Nfft);` es literalmente el periodograma modificado de la ecuación anterior, con Nfft = 4096 (zero-padding sobre N = 1024). El nivel de ruido en nuestras gráficas en dB salió a ≈ −7 dB y no −5 dB por el detalle explicado en la pregunta P3 de esta sección.

### Preguntas

**P1. ¿Qué es el periodograma en una frase?** R: el módulo cuadrado de la FFT de los datos, dividido entre N; equivale a la transformada de la autocorrelación estimada.

**P2. ¿Por qué duplicar las muestras no suaviza el periodograma?** R: porque cada punto de frecuencia sigue siendo una variable exponencial con desviación igual a su media; con más muestras solo tienes más puntos, no más promedio por punto.

**P3. ¿Por qué el piso de ruido salió en −7 dB y no en −5 dB si la DEP es 0.318?** R: el periodograma de ruido en cada frecuencia es exponencial, y el promedio del *logaritmo* de una exponencial queda 2.5 dB por debajo del logaritmo de su promedio: 10·log₁₀(e^{−γ}) = −2.51 dB (γ es la constante de Euler). −5.0 − 2.5 ≈ −7.5 dB, y medimos −7.1 a −7.3 dB.

**P4. ¿Qué corrige el periodograma modificado respecto al simple y qué no?** R: reduce la fuga (eligiendo una ventana de lóbulos laterales bajos), pero no reduce la varianza. Para eso hay que promediar (sección 8).

## 7. Las ventanas (parte investigativa)

**Idea central:** una ventana es una curva que pondera las muestras antes de la FFT. No existe la ventana perfecta: mejorar la fuga siempre cuesta resolución.

### 7.1 El problema que resuelven

Tomar N muestras sin más es multiplicar por una ventana rectangular, cuyo espectro (núcleo de Dirichlet, sección 5.3) tiene lóbulos laterales altos: el primero a solo −13 dB del principal. Cada línea espectral fuerte “derrama” energía por esos lóbulos y puede tapar a una línea débil vecina. La causa física es el **corte brusco** de la señal en los extremos. Si la ventana baja suavemente hacia cero en los bordes, el corte deja de ser brusco y los lóbulos laterales bajan.

### 7.2 Las cuatro ventanas y su fórmula

Todas son del tipo “suma de cosenos”, w\[n\] = a₀ − a₁·cos(2πn/(N−1)) + a₂·cos(4πn/(N−1)), con 0 ≤ n ≤ N−1:

| Ventana | a₀ | a₁ | a₂ | Fórmula |
| --- | --- | --- | --- | --- |
| Rectangular | 1 | 0 | 0 | w\[n\] = 1 |
| Hanning (Hann) | 0.5 | 0.5 | 0 | 0.5 − 0.5·cos(2πn/(N−1)) |
| Hamming | 0.54 | 0.46 | 0 | 0.54 − 0.46·cos(2πn/(N−1)) |
| Blackman | 0.42 | 0.5 | 0.08 | 0.42 − 0.5·cos(2πn/(N−1)) + 0.08·cos(4πn/(N−1)) |

Detalles que suelen confundir:

- **Hann y Hanning** son la misma ventana (lleva el nombre de Julius von Hann; “Hanning” es el nombre que se popularizó por analogía con “Hamming”). Según la documentación de MATLAB, `hann(N)` vale 0 en los extremos mientras que `hanning(N)` usa N+1 en el denominador y no los incluye: la diferencia es mínima.
- **Hamming** es una Hann “elevada”: no llega a 0 en los bordes (vale 0.08). Los coeficientes 0.54 y 0.46 son una aproximación de 25/46 y 21/46, elegidos para anular el primer lóbulo lateral.
- **Blackman**: 0.42, 0.5 y 0.08 son una aproximación redondeada de los coeficientes “exactos” (≈ 0.4266, 0.4966 y 0.0769). Octave y MATLAB usan los redondeados.
- Con la fórmula de Blackman, Octave devuelve −1.4·10⁻¹⁷ en los extremos por redondeo; si luego se usa `pwelch`, hay que recortar con `max(blackman(N), 0)` porque rechaza valores negativos (error que corregimos en la tarea).

### 7.3 Por qué esa suma de cosenos mejora los lóbulos laterales

Un coseno en el tiempo desplaza el espectro: cos(2πn/(N−1)) multiplica a w en el tiempo, y en frecuencia suma copias del espectro de la rectangular desplazadas ±1 bin (y ±2 bins para el segundo coseno). El espectro de la ventana es la suma de núcleos de Dirichlet desplazados y con pesos a₀, −a₁/2, … de signos alternados:

- En el **centro**, las copias se suman y forman un lóbulo principal “montado” sobre varios bins, por eso es **más ancho**: 2/N (rectangular), 4/N (Hann y Hamming: copias a ±1 bin), 6/N (Blackman: copias a ±1 y ±2 bins).
- En los **laterales**, las copias tienen signos opuestos y se cancelan parcialmente, por eso son **más bajos**. Hamming escoge sus pesos justo para que el primer lóbulo lateral se cancele casi por completo.

Por eso el compromiso es inevitable: lo que se gana cancelando laterales se paga ensanchando el principal.

### 7.4 Qué se mide en una ventana (los criterios de Harris, 1978)

El artículo clásico de Harris (“On the Use of Windows for Harmonic Analysis with the Discrete Fourier Transform”, Proc. IEEE, 1978) compara ventanas con tres criterios: el ancho de banda equivalente de ruido (ENBW), el nivel máximo de lóbulo lateral y la velocidad de caída de los laterales. La práctica pide dos: ancho del lóbulo principal y nivel del lóbulo secundario. Todos se miden así:

- **Ancho del lóbulo principal (nulo a nulo):** distancia entre los dos primeros ceros alrededor de ω = 0.
- **Ancho a −3 dB:** ancho donde la magnitud cae a 0.707 del pico (ancho “útil” para decidir si dos picos se separan).
- **Lóbulo secundario (sidelobe):** el mayor pico fuera del lóbulo principal, en dB respecto al pico principal.
- **Decaimiento:** cuántos dB por octava (cada vez que se duplica la distancia al centro) bajan los laterales. Depende de qué tan suave llega la ventana a los bordes: una ventana que salta (rectangular, Hamming) cae 6 dB/octava; una que llega a 0 con derivada cero (Hann, Blackman) cae 18 dB/octava.
- **Ganancia coherente (CG):** CG = (1/N)·Σw\[n\]. Es cuánto se atenua una sinusoide. Hann 0.5, Hamming 0.54, Blackman 0.42.
- **ENBW (ancho de banda equivalente de ruido):** el ancho, en bins, de un filtro rectangular ideal que dejaría pasar la misma potencia de ruido blanco que la ventana:

```latex
\text{ENBW} = \frac{N\sum_n w^2[n]}{\left(\sum_n w[n]\right)^2}
```

### 7.5 Tabla de características (valores medidos por nosotros, N = 1024)

Los anchos y niveles los medimos en el código de Octave con una FFT de 2^18 puntos; ENBW y CG se calcularon de la fórmula con N = 1024:

| Ventana | Lóbulo principal (nulo a nulo) | Ancho a −3 dB | Lóbulo secundario (dB) | Decaimiento | ENBW (bins) | CG |
| --- | --- | --- | --- | --- | --- | --- |
| Rectangular | 2.00/N | 0.89/N | −13.3 | 6 dB/oct | 1.00 | 1.00 |
| Hanning | 4.01/N | 1.45/N | −31.5 | 18 dB/oct | 1.50 | 0.50 |
| Hamming | 4.01/N | 1.30/N | −42.7 | 6 dB/oct | 1.36 | 0.54 |
| Blackman | 6.01/N | 1.65/N | −58.1 | 18 dB/oct | 1.73 | 0.42 |

Los anchos de “nulo a nulo” están en ciclos/muestra; en rad/muestra se multiplican por 2π (el rectangular: 4π/N). Los valores teóricos 2/N, 4/N y 6/N coinciden con los medidos.

**Nota de honestidad sobre las fuentes.** Confirmé en línea que los tres criterios de arriba son los de Harris (1978), pero no pude abrir su tabla original ni otras tablas de ventanas (los sitios estaban bloqueados desde este entorno). Por eso la tabla de arriba son **nuestras mediciones y cálculos**, no valores copiados de un libro. En textos aparecen valores un poco distintos (por ejemplo −43 dB o −41 dB para Hamming) según si se mide el pico del primer lóbulo lateral o el nivel de la ventana “exacta” o la de coeficientes redondeados.

### 7.6 Ejemplos resueltos

**Ejemplo 1. Altura del pico de una sinusoide.** Una sinusoide A·cos(ω₀n) es (A/2)(e^{jω₀n} + e^{−jω₀n}). Cerca de ω₀, |X| ≈ (A/2)·Σw, así que |X|² ≈ A²(Σw)²/4. Con la escala de `periodogram` (sección 6.6):

```latex
P_{pico} \approx \frac{2}{2\pi\sum w^2}\cdot\frac{A^2(\sum w)^2}{4} = \frac{A^2}{4\pi}\cdot\frac{(\sum w)^2}{\sum w^2} = \frac{A^2\,N}{4\pi\,\text{ENBW}}
```

Para A = 1 y N = 1024: Rectangular ≈ 81.5, Hanning ≈ 54.3, Hamming ≈ 59.8, Blackman ≈ 47.1. Eso explica por qué en la gráfica 1.a el pico es más alto con la rectangular (ENBW menor) y más bajo con Blackman (ENBW mayor): una ventana que atenua más los bordes “ve” menos amplitud.

**Ejemplo 2. ¿Se separan los dos picos?** Regla práctica: dos picos se ven separados si Δf es mayor que el ancho a −3 dB de la ventana (y claramente si supera el ancho de nulo a nulo). Con Δf = 0.0119:

| N | Δf·N | Rectangular (0.89/N) | Hanning (1.45/N) | Blackman (6/N nulo a nulo) |
| --- | --- | --- | --- | --- |
| 128 | 1.5 | Sí: 0.89 < 1.5 | En el límite: 1.45 ≈ 1.5, se funden | No: 6 > 1.5 |
| 1024 | 12.2 | Sí | Sí | Sí: 6 < 12.2 |

Por eso con N = 128 solo la rectangular separaba los dos picos, y con N = 1024 (el código final) todas lo hacen.

**Ejemplo 3. ¿Qué ventana elijo?**

| Si lo que quiero es… | Conviene | Por qué |
| --- | --- | --- |
| Separar dos frecuencias muy cercanas de amplitud parecida | Rectangular | Lóbulo principal más estrecho |
| Ver una componente débil junto a una fuerte | Hanning o Blackman | Lóbulos laterales bajos que caen 18 dB/oct |
| Un compromiso general | Hanning o Hamming | Resolución moderada y fuga razonable |
| Medir la amplitud exacta de una sinusoide | Ventana de cima plana (no vista en la unidad) | La ganancia del pico varía poco si la frecuencia cae entre bins |

### Preguntas

**P1. ¿Por qué una ventana más suave tiene lóbulos laterales más bajos?** R: porque deja de cortar la señal bruscamente; la discontinuidad es lo que produce los laterales altos y de lento decaimiento. En frecuencia, la suma de núcleos desplazados cancela los laterales.

**P2. ¿Por qué las ventanas Hann y Blackman “decaen” más rápido que Hamming si Hamming tiene el lateral más bajo?** R: Hamming no llega a 0 en los bordes (queda un salto de 0.08), así que lejos del centro sus laterales caen solo 6 dB/octava. Hann y Blackman llegan suavemente a 0 y caen 18 dB/octava. Cerca del lóbulo principal gana Hamming; lejos, ganan Hann y Blackman. Eso lo vimos en el experimento 1.b sin ruido.

**P3. ¿Qué pasa con el ancho del lóbulo principal si duplico N?** R: se reduce a la mitad (es proporcional a 1/N), y con eso mejora la resolución. El nivel de los lóbulos laterales no cambia: depende solo de la forma de la ventana.

**P4. ¿Por qué se normaliza por Σw² y no por Σw?** R: la DEP es una *potencia* (cuadrática en la señal), así que la corrección también debe ser cuadrática: Σw² mide la potencia que la ventana deja pasar de un ruido blanco. Normalizar por Σw serviría para corregir amplitudes de picos, no niveles de ruido.

## 8. El periodograma promedio (Bartlett y Welch)

**Idea central:** como un solo periodograma tiene varianza que no baja, se calculan varios con trozos distintos de la señal y se promedian, igual que promediar varias mediciones reduce el error. El precio es perder resolución.

### 8.1 El método

Se corta la señal de N muestras en K segmentos de longitud L, separados D muestras entre inicios. Se calcula el periodograma modificado de cada uno y se promedian:

```latex
x_i[n] = x[n+iD],\ n=0,\ldots,L-1,\ i=0,\ldots,K-1 \qquad \hat P_i(\omega)=\frac{1}{LU}\left|\sum_{n=0}^{L-1}w[n]\,x_i[n]\,e^{-j\omega n}\right|^2 \qquad \hat P_{avg}(\omega)=\frac{1}{K}\sum_{i=0}^{K-1}\hat P_i(\omega)
```

**Número de segmentos.** Con solape de S muestras entre segmentos vecinos, D = L − S y:

```latex
K = \left\lfloor \frac{N-L}{L-S} \right\rfloor + 1
```

| N | L | Solape | K | Origen |
| --- | --- | --- | --- | --- |
| 1024 | 256 | L/2 = 128 | 7 | Práctica de 2022 de cursos anteriores |
| 1024 | 64 | L/2 = 32 | 31 | Práctica de 2022 de cursos anteriores |
| 16384 | 1024 | 0 | 16 | Nuestra tarea (experimento 2, código de Octave) |
| 16384 | 1024 | 512 | 31 | Misma señal con solape del 50 % |

Las dos variantes tienen nombre:

- **Bartlett** (1948/1950): segmentos sin solape y ventana rectangular.
- **Welch** (1967): segmentos que pueden solaparse (típicamente 50 %) y con una ventana (periodogramas modificados). Es el que implementa `pwelch`. Por defecto, según la documentación citada en las prácticas de cursos anteriores, `pwelch` usa ocho segmentos con 50 % de solape y ventana Hamming.

### 8.2 Por qué la varianza baja como 1/K

Es el mismo argumento de la sección 4.3 (varianza de la media). Si los K periodogramas fueran independientes y cada uno tuviera varianza relativa 1, su promedio tendría:

```latex
\text{Var}\left[\hat P_{avg}\right] = \frac{\text{Var}\left[\hat P_i\right]}{K} \qquad\Longrightarrow\qquad \frac{\text{Var}}{(\text{media})^2} \approx \frac{1}{K}
```

Además, ahora sí hay **consistencia** si dejamos crecer K (con N → ∞ y L fijo, K → ∞ y la varianza → 0). Con solape los segmentos están correlacionados, así que el beneficio por segmento es menor que el de segmentos independientes, pero se compensa porque caben más segmentos; en la práctica el solape reduce un poco más la varianza, especialmente con ventanas suaves (que atenuan las muestras de los bordes, y el solape las recupera).

### 8.3 El precio: resolución

Cada segmento tiene solo L muestras, así que el ancho del lóbulo principal pasa de ser proporcional a 1/N a ser proporcional a 1/L. Con N fijo y sin solape, K·L ≈ N: **más segmentos significa segmentos más cortos, es decir menos resolución.** Es un compromiso entre varianza y sesgo (la resolución limitada es un sesgo, sección 6.3):

| L | K (N = 1024, sin solape) | Varianza relativa ≈ 1/K | Ancho a −3 dB de la rectangular ≈ 0.89/L |
| --- | --- | --- | --- |
| 1024 | 1 | 1.00 | 0.00087 |
| 256 | 4 | 0.25 | 0.0035 |
| 64 | 16 | 0.0625 | 0.0139 |

Con Δf = 0.0119 ciclos/muestra entre las sinusoides de la práctica, el caso L = 64 tiene un ancho (0.0139) mayor que Δf: **ya no puede separar las dos sinusoides**, aunque la curva sea mucho más suave. Con L = 256 (0.0035) las separa sin problema.

> **Ojo con un error de las prácticas de cursos anteriores.** En la práctica de 2022 que revisamos se dice que con L = 64 y K = 31 la estimación es menos precisa “porque al reducir las muestras de cada segmento la varianza aumenta”. Es lo contrario: la varianza *baja* (hay más segmentos que promediar). Lo que empeora es la **resolución**. No copies esa explicación.

### 8.4 Relación con las señales no estacionarias

El mismo corte en segmentos sirve cuando la señal no es estacionaria (voz). Si en vez de promediar se grafican los K espectros uno al lado del otro en función del tiempo, se obtiene la **evolución temporal del espectro** (`waterfall` o espectrograma), como en la práctica de voz del Grupo 5.

> **En nuestra tarea (experimento 2).** `[PxxRect, wProm] = pwelch(ruidoLargo, wRect, 0, Nfft, 2*pi);` corta 16384 muestras de ruido en K = 16 segmentos de 1024, sin solape, con la ventana indicada. La varianza relativa de la curva bajó de ≈ 0.85–0.93 (un segmento) a 0.056–0.057 (promedio), contra el 1/16 = 0.0625 esperado. En una prueba adicional hecha en Python (4096 muestras, L = 128) el solape de 50 % redujo la varianza de ≈ 0.027 a ≈ 0.013 con Hanning, Hamming y Blackman, y de 0.028 a 0.022 con la rectangular.

### Preguntas

**P1. ¿Cuántos segmentos salen con N = 1024, L = 256 y solape L/2?** R: K = ⌊(1024 − 256)/(256 − 128)⌋ + 1 = ⌊768/128⌋ + 1 = 6 + 1 = 7.

**P2. ¿Por qué promediar reduce la varianza pero ensancha los picos?** R: reduce la varianza porque promedia K estimaciones casi independientes; ensancha los picos porque cada estimación usa solo L < N muestras y el lóbulo principal es proporcional a 1/L.

**P3. ¿Diferencia entre Bartlett y Welch?** R: Bartlett usa segmentos sin solape y ventana rectangular; Welch permite solape y usa una ventana, así que sus periodogramas son modificados.

**P4. Si quiero detectar dos sinusoides muy cercanas y también una curva suave, ¿puedo tener ambas cosas?** R: no del todo. Resolución (L grande) y suavidad (K grande) compiten porque K·L ≈ N. Con una señal más larga (N mayor) se pueden mejorar las dos a la vez.

## 9. La práctica explicada: lo que hicimos y por qué

**Idea central:** la práctica es el laboratorio de todo lo anterior. Cada experimento comprueba una frase de la teoría.

### 9.1 Qué pide el enunciado y con qué parte de la teoría se conecta

| Lo que pide la práctica | Qué teoría comprueba | Sección de esta guía |
| --- | --- | --- |
| Investigar y graficar el espectro de las ventanas Rectangular, Hanning, Hamming y Blackman | Ancho del lóbulo principal y nivel del lóbulo secundario | 7 |
| 1.a: dos sinusoides de amplitud 1 + ruido, periodograma con cada ventana | La ventana limita la resolución (separar dos frecuencias) | 5.3, 6.5, 7 |
| 1.b: amplitudes 1 y 0.01 + ruido | La fuga y el ruido pueden tapar una componente débil | 3.3, 6.4, 7 |
| 2: periodograma promedio de ruido blanco con cada ventana | Promediar reduce la varianza; la DEP del ruido blanco es plana | 3.3, 4, 8 |
| Comparar los estimadores | Periodograma vs periodograma promedio | 6 y 8 |

### 9.2 Cómo funciona el código (Octave/MATLAB), paso a paso

El código es el mismo patrón en todos los experimentos: **crear la señal, enventanar, FFT, normalizar, graficar en dB.**

**Paso 1. Parámetros y ventanas.**

```matlab
N = 1024;  Nfft = 4096;
wRect = rectwin(N);   wHann = hann(N);   wHamm = hamming(N);   wBlck = max(blackman(N), 0);
```

N es el número de muestras y el largo de las ventanas. Nfft > N hace zero-padding (sección 5.4): solo interpola la curva. `max(blackman(N), 0)` corrige el −1.4·10⁻¹⁷ de redondeo que hace fallar a `pwelch`.

**Paso 2. Espectro de las ventanas.**

```matlab
WfRect = abs(fftshift(fft(wRect, Nfft)));     % FFT, centrar el cero, magnitud
fVentana = linspace(-pi, pi, Nfft);           % eje en rad/muestra
```

Es la transformada de la ventana (sección 5.3). `fftshift` solo reordena para que ω = 0 quede al centro; `abs` porque la FFT es compleja.

**Paso 3. La señal del experimento 1.**

```matlab
n = 0:N-1;   w1 = 2*pi/12;   w2 = 2*pi/14;
ruido = randn(1, N);                          % ruido blanco, media 0, varianza 1
xa = 1*sin(w1*n) + 1*sin(w2*n) + ruido;       % 1.a  (en 1.b: 0.01 en la segunda)
```

Son los tres bloques de la DEP del ejemplo 3 (sección 3.3): dos líneas y un piso plano. `randn` es el ruido blanco gaussiano. La *misma* realización de ruido se usa con las cuatro ventanas para que las diferencias se deban solo a la ventana.

**Paso 4. Periodograma con cada ventana.**

```matlab
[Pxx, w] = periodogram(xa, ventanas{i}, Nfft);
```

Es el periodograma modificado de la sección 6.6: multiplica por la ventana, FFT, normaliza por Σw² y devuelve la DEP unilateral en rad/muestra.

**Paso 5. Periodograma promedio.**

```matlab
ruidoLargo = randn(1, N*16);                                  % K = 16 segmentos
[PxxRect, wProm] = pwelch(ruidoLargo, wRect, 0, Nfft, 2*pi);  % solape 0, Fs = 2*pi
```

Corta en segmentos del largo de la ventana (1024), calcula el periodograma de cada uno y promedia (sección 8). El `2*pi` final pone el eje en rad/muestra y el nivel igual al de `periodogram`; sin él, Octave usa Fs = 1 y los dos experimentos no se pueden comparar (problema 2 de la tabla 9.4).

### 9.3 Resultados y qué explican

**Las ventanas.** La figura 1 muestra el espectro en escala lineal: el pico de cada ventana es la suma de sus coeficientes (1024, 512, ≈ 553 y ≈ 430) y la rectangular es la más estrecha. La figura 2 (dB) deja ver los lóbulos laterales; por eso se graficó en dB.

Figura 1 (descripción): espectro de las cuatro ventanas en escala lineal con zoom al lóbulo principal. El pico de la rectangular llega a 1024, Hanning ≈ 512, Hamming ≈ 553 y Blackman ≈ 430; la rectangular es la más estrecha.

Figura 2 (descripción): espectro de las ventanas en dB normalizado a 0 dB. La rectangular tiene laterales cerca de −13 dB que casi no decaen; Hamming ≈ −43 dB que decae poco; Hanning y Blackman decaen rápido, llegando bajo −100 dB lejos del lóbulo principal.

**Experimento 1.a (A₁ = A₂ = 1).** Con N = 1024 las cuatro ventanas muestran los dos picos separados, porque Δf·N = 12.2 es mayor que el lóbulo más ancho (6). Las alturas siguen la fórmula A²N/(4π·ENBW): 81, 54, 60 y 47.

Figura 3 (descripción): periodograma de la señal con A1 = A2 = 1 y ruido, N = 1024, con cada ventana. Se ven dos picos separados, en ω ≈ 0.449 y 0.524 rad/muestra, con las cuatro ventanas; el pico es más alto con la rectangular (≈ 80–90) y más bajo con Blackman (≈ 50–55).

**Experimento 1.b (A₁ = 1, A₂ = 0.01).** Solo aparece el pico de la sinusoide fuerte. La débil está **40 dB** por debajo de la fuerte, porque (A₂/A₁)² = 10⁻⁴. Su pico queda cerca de −21 dB, mientras el piso de ruido está en −7 dB: la débil queda unos 14 dB por debajo del piso y el ruido la tapa, con *cualquier* ventana.

Figura 4 (descripción): periodograma con A1 = 1 y A2 = 0.01 más ruido, escala lineal. Solo se ve el pico de la sinusoide fuerte en ω ≈ 0.524; la débil es invisible con las cuatro ventanas.

Figura 5 (descripción): el mismo experimento en dB. Se ve el pico fuerte (≈ +19 dB) sobre un piso de ruido en torno a −7 dB con fluctuaciones; la sinusoide débil (pico esperado ≈ −21 dB) queda dentro del ruido.

**Experimento 1.b sin ruido (complemento nuestro).** Para ver el efecto de la *ventana* aislado, se quita el ruido. Entonces solo la fuga de la sinusoide fuerte puede tapar a la débil. La línea roja de la figura es esa fuga; la débil se ve cuando el pico azul sobresale de ella:

| Ventana | La débil sobresale de la fuga | Causa |
| --- | --- | --- |
| Rectangular | 0.6 dB (tapada) | Lóbulos laterales altos, caen solo 6 dB/oct |
| Hamming | 9.7 dB | Lateral bajo, pero cae solo 6 dB/oct |
| Hanning | 31.5 dB | Cae 18 dB/oct |
| Blackman | 39.3 dB | Cae 18 dB/oct, laterales más bajos |

Figura 6 (descripción): experimento 1.b sin ruido, en dB, con N = 1024. Línea azul: señal con las dos sinusoides; línea roja discontinua: solo la fuerte (su fuga). Con Hanning y Blackman aparece un pico claro de la débil a la izquierda del fuerte; con Hamming un pico pequeño; con la rectangular queda tapado por la fuga.

**Experimento 2 (ruido blanco, promedio de K = 16).** La curva queda casi horizontal alrededor de 0.32 (= 1/π, −5 dB), como pide la DEP plana del ruido blanco. Con un solo segmento, la varianza relativa es ≈ 0.9; promediando baja a ≈ 0.057 (esperado 1/16 = 0.0625).

Figura 7 (descripción): periodograma promedio (pwelch, K = 16) de ruido blanco con cada ventana, escala lineal. Curvas casi horizontales alrededor de 0.32 (= 1/π) con pequeñas variaciones.

### 9.4 Problemas que encontramos (y por qué pasaron)

| # | Problema | Causa | Solución |
| --- | --- | --- | --- |
| 1 | `pwelch` fallaba con Blackman: “window vector must be real and >=0” | `blackman(N)` devuelve −1.4·10⁻¹⁷ en los extremos por redondeo y `pwelch` rechaza valores negativos | `max(blackman(N), 0)` |
| 2 | Eje y nivel distintos entre `periodogram` y `pwelch` | Sin frecuencia de muestreo, `pwelch` usa Fs = 1: eje en ciclos/muestra (0 a 0.5) y nivel ≈ 2, no ≈ 0.32 | Pasar `2*pi` como quinto argumento |
| 3 | Eje del espectro de ventanas rotulado “(π rad/muestra)” | El eje estaba en radianes (−π a π) | Rotular “(rad/muestra)” |
| 4 | No se ve el nivel de los lóbulos laterales | La escala lineal los aplasta | Gráfica extra en dB, normalizada |
| 5 | Los resultados cambiaban en cada corrida | El ruido es aleatorio | Semilla fija `randn('state', 1)` |
| 6 | Se dijo que la débil estaba “80 dB por debajo” | Error de cálculo: (0.01)² = 10⁻⁴ son **40 dB** | Corregido en todos los documentos |
| 7 | El margen de la débil en la rectangular salió 12.8 dB, pero en la figura no se ve | Se comparó con un solo punto de la fuga, que oscila entre picos y valles | Comparar con el *máximo* de la fuga en esa zona: 0.6 dB |
| 8 | Con N = 128 solo la rectangular separaba las sinusoides; con N = 1024, todas | La resolución depende de N (lóbulo ∝ 1/N) | Se documentó y se usó N = 1024 (igual a la práctica de 2022) |

### 9.5 Qué entregamos y cómo se compara con las prácticas de cursos anteriores

- **Informe en Word** con el formato de las prácticas de la UCAB: encabezado, objetivos, actividad (descripción de funciones), procedimiento paso a paso con código y figuras, análisis de resultados, conclusiones y código completo.
- **Código** (`practica_ventanas.m`, limpio para entregar; `codigo_octave_corregido.m`, con las marcas de cambios).
- **Figuras de Octave** (`figuras_octave/`).

Las prácticas de 2022 (sinusoides + ruido, N = 1024, periodograma y `pwelch` con L = 256 y L = 64) son las más cercanas a la nuestra, pero **no comparan ventanas**; la del Grupo 5 de 2024 estima el espectro de una señal de voz. Lo nuevo en nuestra práctica es la comparación de las cuatro ventanas, la escala en dB, las tablas con valores medidos y la explicación de por qué ocurre cada resultado.

### Preguntas

**P1. ¿Por qué en el experimento 1.b no se ve la segunda sinusoide con ninguna ventana?** R: porque está 40 dB por debajo de la primera y su pico queda unos 14 dB por debajo del nivel medio del ruido. Una ventana solo cambia la fuga; no puede subir un pico por encima del ruido.

**P2. ¿Qué se logra quitando el ruido en 1.b?** R: se ve el efecto de la ventana por sí solo. Sin ruido, el único obstáculo es la fuga de la sinusoide fuerte, y Hanning y Blackman, con laterales que caen 18 dB/octava, dejan ver a la débil; la rectangular no.

**P3. ¿Por qué las cuatro ventanas dan casi la misma curva en el experimento 2?** R: la DEP del ruido blanco es plana, así que no hay líneas que puedan ensancharse ni fugas que importen, y la normalización por Σw² hace que todas den el mismo nivel medio. La ventana solo importa cuando la señal tiene estructura espectral.

**P4. ¿Qué cambiaría si en vez de N = 1024 usamos N = 128?** R: los lóbulos principales serían 8 veces más anchos, así que solo la rectangular separaría las dos sinusoides de 1.a; la sinusoide débil de 1.b quedaría aún más por debajo del ruido relativo (el pico crece ∝ N).

## 10. Preguntas para autoevaluarte

Tapa la respuesta, intenta contestar y luego compara. Están en orden: conceptos, cálculos y razonamiento sobre la práctica.

### Conceptos

**1. ¿Qué significa que un proceso sea estacionario en sentido amplio?** R: que su media y su varianza no dependen del tiempo, y que su autocorrelación depende solo del desfase m = n₂ − n₁, no de los instantes concretos.

**2. ¿Qué pide la ergodicidad que no pida la estacionariedad?** R: que los promedios temporales de una realización coincidan con los del ensamble. Es la propiedad que permite estimar con una sola señal. Ergódico implica WSS, pero no al revés.

**3. ¿Qué es la DEP y cómo se obtiene?** R: la transformada de Fourier de la autocorrelación (Wiener-Khinchin). Su área (dividida entre 2π) es la potencia promedio total.

**4. ¿Cuál es la diferencia entre sesgo y varianza de un estimador?** R: el sesgo es cuánto se aleja el *promedio* del estimador del valor verdadero; la varianza es cuánto se dispersa el estimador alrededor de su propio promedio. Un estimador puede tener uno, otro, ambos o ninguno.

**5. ¿Por qué se dice que el periodograma no es consistente?** R: porque su varianza en cada frecuencia es ≈ (DEP)² aunque aumente N, así que no tiende a 0.

**6. ¿Qué problema resuelve la ventana y cuál no?** R: reduce la fuga espectral (lóbulos laterales). No reduce la varianza: cada punto del periodograma sigue siendo una variable exponencial.

**7. ¿Qué problema resuelve el promedio y qué cuesta?** R: reduce la varianza en un factor ≈ 1/K. Cuesta resolución, porque los segmentos son más cortos (el lóbulo principal ∝ 1/L).

**8. ¿Qué ventana tiene el lóbulo principal más estrecho y cuál los laterales más bajos?** R: la rectangular (2/N) tiene el principal más estrecho; Blackman (−58 dB) tiene los laterales más bajos, a costa de un principal de 6/N.

### Cálculos

**9. Una sinusoide de amplitud 3, ¿qué potencia tiene?** R: A²/2 = 9/2 = 4.5.

**10. Convierte a dB: razón de potencias 100; razón de amplitudes 0.1; razón de amplitudes 0.01.** R: potencia: 10·log₁₀(100) = **20 dB**. Amplitud 0.1: 20·log₁₀(0.1) = **−20 dB**. Amplitud 0.01: 20·log₁₀(0.01) = **−40 dB** (la razón de potencias es 10⁻⁴, y 10·log₁₀(10⁻⁴) = −40; mismo resultado). Error típico: usar 10 en vez de 20 con amplitudes y obtener −20 en vez de −40.

**11. La media estimada de L = 400 muestras de ruido de varianza 4, ¿qué desviación estándar tiene?** R: Var = σ²/L = 4/400 = 0.01, así que la desviación es √0.01 = 0.1.

**12. ¿Cuántos segmentos salen con N = 2000, L = 200 y solape de 100 muestras?** R: K = ⌊(2000 − 200)/(200 − 100)⌋ + 1 = ⌊1800/100⌋ + 1 = 19.

**13. ¿Cuántos segmentos necesito para que la varianza relativa del promedio sea 1 %?** R: 1/K = 0.01, así que K = 100 (si son independientes).

**14. Con fs = 8000 Hz y N = 512, ¿cuál es la resolución aproximada de la ventana rectangular?** R: ancho a −3 dB ≈ 0.89/N = 0.89/512 = 0.00174 ciclos/muestra. En Hz: 0.00174 · 8000 ≈ **13.9 Hz**. La relación es f = k·fs/N: cada bin de la DFT mide fs/N = 15.6 Hz.

**15. Dos tonos a 1000 y 1030 Hz con fs = 8000 Hz. ¿Se separan con N = 256? ¿Y con N = 1024? (rectangular y Hanning)** R: Δf = 30/8000 = 0.00375 ciclos/muestra. Con N = 256, Δf·N = 0.96 bins: la rectangular (0.89 bins) los separa justo; Hanning (1.45 bins) los funde. Con N = 1024, Δf·N = 3.84 bins, mayor que ambos anchos: los dos los separan.

**16. ¿Cuánto vale el nivel de ruido blanco de varianza 1 que verás en `periodogram`, en lineal y en dB?** R: σ²/π ≈ 0.318, es decir −5.0 dB como valor verdadero; el promedio de la curva en dB sale ≈ −7.5 dB (−2.5 dB por promediar un logaritmo de una exponencial).

### Razonamiento sobre la práctica

**17. En 1.a con N = 1024 se ven dos picos con todas las ventanas. Si bajamos a N = 128, ¿qué pasa con Hanning y por qué?** R: Δf·N = 1.5 bins y el ancho a −3 dB de Hanning es 1.45 bins: están al límite, y los picos se funden en una meseta. Es consecuencia de que la ventana convoluciona cada línea con su lóbulo principal.

**18. En 1.b, ¿por qué una ventana de laterales muy bajos como Blackman no ayuda?** R: porque lo que tapa a la débil es el *ruido*, no la fuga; la fuga solo importa cuando es mayor que el ruido o que la componente débil. Sin ruido, Blackman sí ayuda (margen de 39 dB).

**19. ¿Por qué en el experimento 2 no hace falta una ventana especial?** R: porque la DEP del ruido blanco es plana y no hay líneas que ensanchar ni componentes débiles que tapar. La ventana importa para señales con estructura espectral.

**20. ¿Qué ocurriría con la media del periodograma promedio si no se normalizara por Σw²?** R: dependería de la ventana: Blackman (que atenua mucho los bordes, Σw² ≈ 0.3·N) daría un nivel ≈ 3 veces menor que la rectangular (Σw² = N). La normalización por Σw² hace que todas den ≈ 0.32.

**21. Un compañero dice: “con L = 64 y K = 31 la varianza aumenta”. ¿Está bien?** R: no. Con más segmentos la varianza *baja* (≈ 1/K). Lo que empeora con L = 64 es la resolución (0.89/64 = 0.0139 > Δf = 0.0119), por eso ya no se separan las dos sinusoides.

**22. Una señal de voz de 3 s, ¿por qué se divide en tramos de 5 a 100 ms?** R: porque la voz no es estacionaria en intervalos largos, pero sí aproximadamente en intervalos cortos; en cada tramo la DEP tiene sentido y se puede estimar.

## 11. Glosario, errores comunes y fuentes

### 11.1 Glosario

| Término | Definición |
| --- | --- |
| Variable aleatoria | Función que asigna un número a cada resultado de un experimento aleatorio |
| Proceso estocástico | Familia de variables aleatorias, una por cada instante n |
| Realización | Una señal concreta de un proceso aleatorio |
| Ensamble | El conjunto de todas las realizaciones posibles |
| Valor esperado (media) | Promedio ponderado por probabilidad de los valores posibles: E\[X\] |
| Varianza | E\[\|X − E\[X\]\|²\]: potencia de la parte AC |
| Autocorrelación | E\[X\_{n+m} X\_n\*\]: qué tan parecida es la señal a sí misma desplazada m muestras |
| Estacionario (sentido amplio) | Media y varianza constantes; autocorrelación depende solo del desfase |
| Ergódico | Los promedios temporales de una realización coinciden con los del ensamble |
| DEP | Densidad espectral de potencia: transformada de Fourier de la autocorrelación |
| Ruido blanco | Muestras no correlacionadas; DEP constante σ² |
| Estimador | Función de las muestras que aproxima un parámetro; es una variable aleatoria |
| Sesgo | E\[estimador\] − valor verdadero |
| Consistente | Su varianza tiende a 0 al aumentar las muestras |
| ECM | Error cuadrático medio = sesgo² + varianza |
| FFT / DFT | Algoritmo rápido / transformada que calcula muestras equiespaciadas del espectro |
| Bin | Cada punto de la DFT; separados fs/N Hz (o 1/N ciclos/muestra) |
| Zero-padding | Agregar ceros antes de la FFT; interpola la curva, no mejora la resolución |
| Ventana | Curva que pondera las muestras antes de la FFT |
| Lóbulo principal | Pico central del espectro de la ventana; fija la resolución |
| Lóbulo secundario | Picos laterales; fijan la fuga |
| Fuga espectral | Energía de una frecuencia que se derrama a las vecinas por los lóbulos laterales |
| Resolución | Capacidad de separar dos frecuencias cercanas (∝ 1/N) |
| ENBW | Ancho de banda equivalente de ruido: NΣw²/(Σw)² bins |
| Ganancia coherente | (1/N)Σw: cuánto atenua la ventana a una sinusoide |
| Periodograma | Módulo cuadrado de la FFT dividido entre N: estimador directo de la DEP |
| Periodograma modificado | Periodograma con ventana, normalizado por la potencia de la ventana |
| Bartlett / Welch | Promedio de periodogramas de segmentos sin solape y rectangular / con solape y ventana |
| dB | Escala logarítmica: 10·log₁₀ en potencia, 20·log₁₀ en amplitud |
| rad/muestra vs ciclos/muestra | ω = 2πf; ω va de 0 a π y f de 0 a 0.5 en señales reales |

### 11.2 Errores comunes (algunos los cometimos nosotros)

1. **Confundir 10·log₁₀ con 20·log₁₀.** Con potencias es 10, con amplitudes es 20. La sinusoide de amplitud 0.01 está a −40 dB de la de amplitud 1, no a −80 dB (error que cometimos y corregimos).
2. **Creer que más relleno de ceros (Nfft) mejora la resolución.** Solo interpola. La resolución depende de N y de la ventana.
3. **Creer que más muestras suavizan el periodograma simple.** Dan más puntos de frecuencia, no menos varianza por punto.
4. **Decir que el periodograma promedio “mejora todo”.** Reduce la varianza pero baja la resolución.
5. **Pensar que una ventana es “mejor” en general.** Depende de qué se busque: resolución (rectangular) o rango dinámico (Hann, Blackman).
6. **Mezclar los ejes.** `periodogram` devuelve rad/muestra (0 a π); `pwelch` sin Fs en Octave devuelve ciclos/muestra (0 a 0.5) y otro nivel. Usar `2*pi` como Fs.
7. **Esperar que el piso de ruido esté en la DEP verdadera.** En dB el promedio queda 2.5 dB por debajo (promedio de logaritmos).
8. **Comparar la débil con un solo punto de la fuga.** La fuga oscila; hay que compararla con su máximo en la zona.
9. **Olvidar la normalización por Σw².** Sin ella, el nivel depende de la ventana.
10. **Copiar explicaciones de prácticas anteriores sin revisarlas.** Una de 2022 dice que la varianza aumenta con L = 64; es la resolución la que empeora.

### 11.3 Fuentes

**Material de clase y prácticas (de tu Drive, leídos para esta guía):**

- [UIT1\_EB\_301.pdf, *Estimadores de momentos estadísticos*, María Cristi Stefanelli, UCAB](https://drive.google.com/file/d/1hn8w1jb7WTBPEwXQy1-_00_5doFfbYep/view?usp=drivesdk): secciones 1 a 5 de esta guía.
- [Práctica 1 de 2022, Estimador de la densidad espectral clásica usando MatLab](https://drive.google.com/file/d/1iwe-NVL3uVF8lzuupyCvn-fBmZ2z3ElU/view?usp=drivesdk): comparación y K = 7 y 31.
- [P1 Procesamiento de Señales, Grupo 5 (2024), estimación de la DEP de una señal de voz](https://drive.google.com/file/d/10Pgp7DVl3kVGYpT1kDzgWOWFqP8pKEXu/view?usp=drivesdk): formato y descripción de `pwelch`/`periodogram`.
- Los PDF `UIT1_Tutorial_EB_301.pdf` y `UIT2_EB_301.pdf` de la misma carpeta no tienen texto legible (son imágenes o están vacíos); no pude usarlos, revísalos a mano.

**Literatura (parte investigativa):**

- Harris, F. J. (1978). *On the Use of Windows for Harmonic Analysis with the Discrete Fourier Transform*. Proc. IEEE, 66(1). Los tres criterios (ENBW, nivel máximo de lóbulo lateral y velocidad de caída) los verifiqué en [este resumen](https://arxiv.org/pdf/1404.6979); no pude abrir la tabla original.
- Welch, P. D. (1967). *The use of fast Fourier transform for the estimation of power spectra: a method based on time averaging over short, modified periodograms*. IEEE Trans. Audio Electroacoust., 15, 70–73. Cita verificada en [IBM Research](https://research.ibm.com/publications/the-use-of-fast-fourier-transform-for-the-estimation-of-power-spectra-a-method-based-on-time-averaging-over-short-modified-periodograms).
- Bartlett, M. S. (1948, 1950), promedio de periodogramas sin solape. Blackman, R. B. y Tukey, J. W. (1958), *The Measurement of Power Spectra*. Schuster, A. (1898), origen del periodograma. Oppenheim, A. V. y Schafer, R. W., *Discrete-Time Signal Processing* (la clase cita la edición de 2011). Stoica, P. y Moses, R., *Spectral Analysis of Signals*. Estas cinco las cito de memoria, como lecturas recomendadas, sin haberlas abierto aquí.

**Lo que es mío y no viene de ninguna fuente:** las demostraciones de las secciones 4.3, 6.2 y 7.6, los ejemplos, y las mediciones de la práctica (hechas por nosotros en Octave y Python). Los valores numéricos de las tablas los verifiqué ejecutando el código.
