import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

st.set_page_config(page_title="Simulation Laser", layout="centered")

st.title("🔦 Simulation Laser : Profil, Pertes & Émission")

# Paramètres utilisateur
st.sidebar.header("Paramètres du laser")
λ = st.sidebar.number_input("Longueur d’onde λ (nm)", value=1064.0) * 1e-9
w0 = st.sidebar.slider("Rayon du faisceau w0 (mm)", 0.1, 2.0, 0.5) * 1e-3
L = st.sidebar.slider("Longueur cavité L (m)", 0.1, 1.0, 0.3)
R1 = st.sidebar.slider("Réflectivité miroir 1", 0.90, 1.00, 0.99)
R2 = st.sidebar.slider("Réflectivité miroir 2", 0.5, 0.99, 0.90)
alpha_int = st.sidebar.slider("Pertes internes α", 0.0, 0.05, 0.01)
G0 = st.sidebar.slider("Gain G0 (1/m)", 0.0, 0.2, 0.05)
Psat = st.sidebar.slider("Puissance de saturation (W)", 0.1, 5.0, 1.0)

# Calculs
def profil_gaussien(x, w0, P0=1.0):
    return P0 * np.exp(-2 * x*2 / w0*2)

def gain_seuil(R1, R2, alpha_int, L):
    return (1/(2*L)) * (np.log(1/(R1*R2)) + 2 * alpha_int * L)

def puissance_sortie(G0, g_seuil, R2, Psat):
    if G0 <= g_seuil:
        return 0
    return Psat * R2 * ((G0 - g_seuil) / g_seuil)

g_seuil = gain_seuil(R1, R2, alpha_int, L)
P_out = puissance_sortie(G0, g_seuil, R2, Psat)

st.subheader("📈 Profil gaussien")
x = np.linspace(-2*w0, 2*w0, 500)
intensite = profil_gaussien(x, w0)

fig1, ax1 = plt.subplots()
ax1.plot(x*1e3, intensite)
ax1.set_title("Profil gaussien (transverse)")
ax1.set_xlabel("x (mm)")
ax1.set_ylabel("Intensité relative")
st.pyplot(fig1)

st.subheader("📊 Puissance en fonction du gain")
gains = np.linspace(0, 2*g_seuil, 200)
P_list = [puissance_sortie(g, g_seuil, R2, Psat) for g in gains]

fig2, ax2 = plt.subplots()
ax2.plot(gains, P_list, label="P(G)")
ax2.axvline(g_seuil, color='r', linestyle='--', label="Gain seuil")
ax2.set_xlabel("Gain G (1/m)")
ax2.set_ylabel("Puissance (W)")
ax2.set_title("Puissance émise vs Gain")
ax2.legend()
st.pyplot(fig2)

st.success(f"Gain seuil : {g_seuil:.4f} 1/m")
st.success(f"Puissance émise : {P_out:.3f} W")