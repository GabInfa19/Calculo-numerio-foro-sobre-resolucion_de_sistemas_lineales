"""
Estimacion espectral clasica: efecto de la ventana en el periodograma.
Ventanas: Rectangular, Hanning, Hamming, Blackman.
Genera las figuras en ./figuras y las tablas de resultados.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

np.random.seed(1)                      # semilla fija para reproducibilidad
OUT = "figuras/"

# ---------------------------------------------------------------- ventanas
def ventana(nombre, N):
    n = np.arange(N)
    if nombre == "Rectangular":
        return np.ones(N)
    if nombre == "Hanning":            # Hann simetrica
        return 0.5 - 0.5*np.cos(2*np.pi*n/(N-1))
    if nombre == "Hamming":
        return 0.54 - 0.46*np.cos(2*np.pi*n/(N-1))
    if nombre == "Blackman":
        return 0.42 - 0.5*np.cos(2*np.pi*n/(N-1)) + 0.08*np.cos(4*np.pi*n/(N-1))
VENT = ["Rectangular", "Hanning", "Hamming", "Blackman"]
COL = dict(zip(VENT, ["C0", "C1", "C2", "C3"]))

def espectro_ventana(w, NFFT=2**16):
    """|W(f)| normalizado en dB, f en ciclos/muestra (-0.5..0.5)."""
    W = np.fft.fftshift(np.fft.fft(w, NFFT))
    f = np.fft.fftshift(np.fft.fftfreq(NFFT))
    return f, 20*np.log10(np.abs(W)/np.abs(W).max() + 1e-12)

def caracteristicas(w, NFFT=2**18):
    """Ancho de lobulo principal (nulo a nulo, 3 dB) y nivel del lobulo secundario."""
    f, dB = espectro_ventana(w, NFFT)
    m = f >= 0
    f, dB = f[m], dB[m]
    mag = 10**(dB/20)
    k = 1
    while not (mag[k] < mag[k-1] and mag[k] <= mag[k+1]):   # primer nulo
        k += 1
    sl = dB[k:].max()                                        # mayor lobulo lateral
    i3 = np.argmax(dB < -3.0)
    return 2*f[k], 2*f[i3], sl

# ------------------------------------------------- 1) espectro de ventanas
N = 64
fig, ax = plt.subplots(1, 2, figsize=(12, 4))
filas = []
for v in VENT:
    w = ventana(v, N)
    ax[0].plot(w, label=v, color=COL[v])
    f, dB = espectro_ventana(w)
    ax[1].plot(f, dB, label=v, color=COL[v])
    filas.append((v,) + caracteristicas(w))
ax[0].set(title="Fig. 1a. Ventanas en el tiempo (N=64)", xlabel="n (muestras)", ylabel="w[n]")
ax[1].set(title="Fig. 1b. Espectro de las ventanas (N=64)", xlabel="Frecuencia (ciclos/muestra)",
          ylabel="|W(f)| normalizado (dB)", ylim=(-120, 5), xlim=(-0.5, 0.5))
ax[0].legend(); ax[1].legend(); ax[0].grid(); ax[1].grid()
plt.tight_layout(); plt.savefig(OUT+"fig1_ventanas.png", dpi=150); plt.close()

with open("tabla_ventanas.txt", "w") as fh:
    fh.write("Ventana | Ancho nulo-nulo (x1/N) | Ancho 3dB (x1/N) | Lobulo sec. (dB)\n")
    for v, bw, b3, sl in filas:
        fh.write(f"{v} | {bw*N:.2f} | {b3*N:.2f} | {sl:.1f}\n")
print(open("tabla_ventanas.txt").read())

# ---------------------------------------------------------- periodograma
def periodograma(x, w, NFFT=4096):
    """Periodograma modificado: |FFT(x.w)|^2 / (N * U),  U = mean(w^2)."""
    U = np.mean(w**2)
    X = np.fft.fft(x*w, NFFT)
    return np.abs(X)**2 / (len(x)*U)

# ---------------------------------------------- 2) experimentos 1.a y 1.b
Nx = 128
n = np.arange(Nx)
w1, w2 = 2*np.pi/12, 2*np.pi/14
ruido = np.random.randn(Nx)            # ruido blanco N(0,1); misma realizacion en todos
f = np.fft.fftfreq(4096)[:2048]

def experimento(A1, A2, nombre, num, tag):
    x = A1*np.cos(w1*n) + A2*np.cos(w2*n) + ruido
    fig, ax = plt.subplots(2, 2, figsize=(12, 7), sharex=True, sharey=True)
    for a, v in zip(ax.ravel(), VENT):
        P = periodograma(x, ventana(v, Nx))[:2048]
        a.plot(f, 10*np.log10(P+1e-12), color=COL[v])
        for wk in (w1, w2):
            a.axvline(wk/(2*np.pi), color="k", ls=":", lw=0.8)
        a.set(title=f"Ventana {v}"); a.grid()
    for a in ax[1]: a.set_xlabel("Frecuencia (ciclos/muestra)")
    for a in ax[:, 0]: a.set_ylabel("DEP estimada (dB)")
    fig.suptitle(f"Fig. {num}. Periodograma, {nombre}, N={Nx} (lineas punteadas: f1 y f2 reales)")
    plt.tight_layout(); plt.savefig(OUT+f"fig{num}_periodograma_{tag}.png", dpi=150); plt.close()

experimento(1, 1,    "A1=A2=1",       2, "1a")
experimento(1, 0.01, "A1=1, A2=0.01", 3, "1b")

x = np.cos(w1*n) + 0.01*np.cos(w2*n) + ruido
plt.figure(figsize=(9, 4.5))
for v in VENT:
    plt.plot(f, 10*np.log10(periodograma(x, ventana(v, Nx))[:2048]+1e-12), color=COL[v], label=v, lw=1)
plt.axvline(w2/(2*np.pi), color="k", ls=":"); plt.legend(); plt.grid()
plt.xlabel("Frecuencia (ciclos/muestra)"); plt.ylabel("DEP estimada (dB)")
plt.title("Fig. 4. Periodograma, A1=1, A2=0.01 (superpuesto; punteada: f2)")
plt.tight_layout(); plt.savefig(OUT+"fig4_superpuesto_1b.png", dpi=150); plt.close()

# nivel de piso de ruido lejos de las sinusoides (f>0.3) para 1.b
with open("tabla_piso.txt", "w") as fh:
    fh.write("Ventana | piso de ruido medio f>0.3 (dB) | maximo en f>0.3 (dB)\n")
    for v in VENT:
        P = 10*np.log10(periodograma(x, ventana(v, Nx))[:2048])
        fh.write(f"{v} | {P[f>0.3].mean():.1f} | {P[f>0.3].max():.1f}\n")
print(open("tabla_piso.txt").read())

# ----------------------------- 3) periodograma promedio, ruido blanco
Ntot, L = 4096, 128
x = np.random.randn(Ntot)
fL = np.fft.fftfreq(L)[:L//2]

def promedio(x, w, L, solape):
    paso = L - solape
    P = [periodograma(x[i:i+L], w, L) for i in range(0, len(x)-L+1, paso)]
    return np.mean(P, axis=0)[:L//2], len(P)

fig, ax = plt.subplots(2, 2, figsize=(12, 7), sharex=True, sharey=True)
res = []
for a, v in zip(ax.ravel(), VENT):
    w = ventana(v, L)
    Pi, _ = promedio(x[:L], w, L, 0)
    Pb, K = promedio(x, w, L, 0)
    Pw, Kw = promedio(x, w, L, L//2)
    a.plot(fL, 10*np.log10(Pi), color="0.7", label="1 segmento")
    a.plot(fL, 10*np.log10(Pb), color=COL[v], label=f"promedio K={K}")
    a.plot(fL, 10*np.log10(Pw), "k", lw=0.8, label=f"Welch 50% K={Kw}")
    a.axhline(0, color="r", ls="--", lw=0.8)
    a.set(title=f"Ventana {v}"); a.grid(); a.legend(fontsize=8)
    res.append((v, Pi.mean(), Pi.var(), Pb.mean(), Pb.var(), Pw.mean(), Pw.var(), K, Kw))
for a in ax[1]: a.set_xlabel("Frecuencia (ciclos/muestra)")
for a in ax[:, 0]: a.set_ylabel("DEP estimada (dB)")
fig.suptitle("Fig. 5. Periodograma promedio de ruido blanco (varianza 1; DEP real = 0 dB, linea roja)")
plt.tight_layout(); plt.savefig(OUT+"fig5_promedio_ruido.png", dpi=150); plt.close()

with open("tabla_promedio.txt", "w") as fh:
    fh.write("Ventana | media 1seg | var 1seg | media prom | var prom | media Welch50 | var Welch50 | K | Kw\n")
    for r in res:
        fh.write(f"{r[0]} | " + " | ".join(f"{t:.4f}" if isinstance(t, float) else str(t) for t in r[1:]) + "\n")
print(open("tabla_promedio.txt").read())

# -------- complemento: 1.b SIN ruido, para ver la fuga espectral de cada ventana
Nc = 512; nc = np.arange(Nc)
xc = np.cos(w1*nc) + 0.01*np.cos(w2*nc)
fc = np.fft.fftfreq(4096)[:2048]
plt.figure(figsize=(9, 4.5))
for v in VENT:
    plt.plot(fc, 10*np.log10(periodograma(xc, ventana(v, Nc))[:2048]+1e-15), color=COL[v], label=v, lw=1)
for wk in (w1, w2): plt.axvline(wk/(2*np.pi), color="k", ls=":", lw=0.8)
plt.ylim(-120, 40); plt.legend(); plt.grid()
plt.xlabel("Frecuencia (ciclos/muestra)"); plt.ylabel("DEP estimada (dB)")
plt.title("Fig. 6. Periodograma sin ruido, A1=1, A2=0.01, N=512")
plt.tight_layout(); plt.savefig(OUT+"fig6_sin_ruido_1b.png", dpi=150); plt.close()
