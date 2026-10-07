# Cómo usar este paquete en NotebookLM

NotebookLM es una herramienta de Google que hay que usar desde tu cuenta (notebooklm.google.com). Yo no puedo crear el notebook por ti, pero este paquete tiene las fuentes listas para subirlas.

## Paso 1. Crear el notebook y subir las fuentes

1. Entra a notebooklm.google.com y crea un notebook nuevo. Ponle un nombre, por ejemplo "Estimación espectral con ventanas (Procesamiento de Señales)".
2. En "Añadir fuentes", sube estos archivos de la carpeta NotebookLM del repositorio (rama practica-ventanas-periodograma):
   - **NotebookLM_Fuente_Completa.pdf**: todo el contenido en un solo archivo (guía de estudio, explicaciones causales, funciones de Octave, resultados y defensa). Con solo esta fuente ya puedes trabajar.
   - Opcional, para separar el contenido en varias fuentes y citar por tema: 01_Guia_de_estudio.md, 02_Por_que_sucede_cada_cosa.md, 03_Funciones_de_Octave.md y 04_Resultados_y_defensa.md (si NotebookLM no acepta .md en tu cuenta, usa el PDF).
   - **practica_ventanas.txt**: el código de Octave de la práctica.
   - Opcional: el informe en Word (Informe_Estimacion_Espectral.docx), que está en la carpeta Practica_Ventanas.
3. Material de clase, si quieres que las respuestas lo citen: agrega desde tu Google Drive el PDF UIT1_EB_301.pdf (Estimadores de momentos estadísticos, de la profesora Stefanelli), que está en la carpeta "Procesamiento De Señales".

Nota: no subas datos personales que no quieras compartir; las fuentes de este paquete no los contienen.

## Paso 2. Preguntas para hacerle al notebook (en el chat)

- "Explícame paso a paso por qué la ventana rectangular tiene lóbulos laterales de −13 dB y Blackman de −58 dB, usando solo las fuentes."
- "¿Por qué el periodograma simple no es consistente y cómo lo resuelve el periodograma promedio?"
- "Dame un ejemplo numérico con N = 1024 para decidir si dos sinusoides cercanas se separan con cada ventana."
- "¿De dónde sale el piso de ruido de −7 dB? Explícalo con la propiedad de la distribución exponencial."
- "Hazme 15 preguntas de práctica sobre la Unidad I con respuesta, de fácil a difícil."
- "Explícame la derivación de que el periodograma es la transformada de la autocorrelación estimada."
- "Simula una defensa oral: hazme una pregunta a la vez sobre la práctica y corrígeme."

## Paso 3. Herramientas de NotebookLM sugeridas

- **Resumen de audio (Audio Overview).** Pulsa "Personalizar" y pega esta instrucción: "Explica a fondo, como a un estudiante a quien le cuesta la materia, de dónde viene cada concepto: procesos aleatorios, estacionariedad, ergodicidad, DEP, estimadores, periodograma, ventanas y periodograma promedio. Después explica por qué ocurre cada resultado de los experimentos 1.a, 1.b y 2 y gracias a qué propiedad. Usa ejemplos numéricos y analogías."
- **Guía de estudio y preguntas frecuentes.** Genera ambas desde el panel de estudio y compáralas con las preguntas del documento 04.
- **Mapa mental.** Útil para ver la cadena: señal aleatoria → promedios → estacionariedad y ergodicidad → autocorrelación → DEP → estimadores → periodograma → ventana → promedio.
- **Tarjetas de estudio (flashcards) y cuestionario.** Genéralas con el glosario de la guía de estudio.

## Qué contiene cada fuente

| Archivo | Contenido |
|---|---|
| 01_Guia_de_estudio | Teoría desde cero con fórmulas, orígenes, ejemplos resueltos y preguntas con respuesta, basada en la clase UIT1_EB_301 y en la investigación de ventanas y estimadores |
| 02_Por_que_sucede_cada_cosa | Para cada resultado: qué se observó, por qué sucede y gracias a qué teoría |
| 03_Funciones_de_Octave | Explicación de cada función del código |
| 04_Resultados_y_defensa | Tablas de resultados, problemas resueltos, guion de defensa y preguntas probables |
| practica_ventanas.txt | Código de la práctica |

## Aviso sobre la fiabilidad

Los valores numéricos de las tablas se obtuvieron ejecutando el código en Octave y Python (y verificaciones por cálculo). Las citas de Harris (1978) y Welch (1967) se verificaron en línea; otras referencias bibliográficas (Bartlett, Blackman y Tukey, Oppenheim y Schafer, Stoica y Moses) se citan de memoria. Si NotebookLM contradice algún dato, revisa la fuente original.
