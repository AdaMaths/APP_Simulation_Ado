import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

st.set_page_config(page_title="Navier-Stokes", layout="centered")

st.title("🌊 Simulation 2D - Équation de Navier-Stokes (incompressible)")

# Paramètres utilisateur
st.sidebar.header("Paramètres de simulation")
nx = st.sidebar.slider("Points en x", 20, 100, 41)
ny = st.sidebar.slider("Points en y", 20, 100, 41)
nt = st.sidebar.slider("Pas de temps (itérations)", 100, 1000, 500, step=100)
viscosite = st.sidebar.slider("Viscosité ν", 0.01, 1.0, 0.1, step=0.01)
dt = st.sidebar.slider("Pas de temps dt", 0.0005, 0.01, 0.001, step=0.0005)
nit = st.sidebar.slider("Itérations pression", 10, 100, 50)

if st.button("▶ Lancer la simulation"):
    dx = 2 / (nx - 1)
    dy = 2 / (ny - 1)
    rho = 1
    nu = viscosite

    u = np.zeros((ny, nx))
    v = np.zeros((ny, nx))
    p = np.zeros((ny, nx))
    b = np.zeros((ny, nx))

    def build_up_b(b, rho, dt, u, v, dx, dy):
        b[1:-1,1:-1] = rho * (1/dt *
            ((u[1:-1,2:] - u[1:-1,0:-2]) / (2*dx) +
             (v[2:,1:-1] - v[0:-2,1:-1]) / (2*dy)) -
            ((u[1:-1,2:] - u[1:-1,0:-2]) / (2*dx))**2 -
            2 * ((u[2:,1:-1] - u[0:-2,1:-1]) / (2*dy) *
                 (v[1:-1,2:] - v[1:-1,0:-2]) / (2*dx)) -
            ((v[2:,1:-1] - v[0:-2,1:-1]) / (2*dy))**2)
        return b

    def pressure_poisson(p, dx, dy, b):
        for _ in range(nit):
            pn = p.copy()
            p[1:-1,1:-1] = (((pn[1:-1,2:] + pn[1:-1,0:-2]) * dy**2 +
                             (pn[2:,1:-1] + pn[0:-2,1:-1]) * dx**2) /
                            (2 * (dx*2 + dy*2)) -
                            dx*2 * dy2 / (2 * (dx2 + dy*2)) *
                            b[1:-1,1:-1])
            p[:,-1] = p[:,-2]
            p[:,0] = p[:,1]
            p[0,:] = p[1,:]
            p[-1,:] = 0
        return p

    # Simulation
    for n in range(nt):
        un = u.copy()
        vn = v.copy()

        b = build_up_b(b, rho, dt, u, v, dx, dy)
        p = pressure_poisson(p, dx, dy, b)

        u[1:-1,1:-1] = (un[1:-1,1:-1] -
                        un[1:-1,1:-1] * dt / dx *
                       (un[1:-1,1:-1] - un[1:-1,0:-2]) -
                        vn[1:-1,1:-1] * dt / dy *
                       (un[1:-1,1:-1] - un[0:-2,1:-1]) -
                        dt / (2 * rho * dx) *
                       (p[1:-1,2:] - p[1:-1,0:-2]) +
                        nu * (dt / dx**2 *
                       (un[1:-1,2:] - 2 * un[1:-1,1:-1] + un[1:-1,0:-2]) +
                             dt / dy**2 *
                       (un[2:,1:-1] - 2 * un[1:-1,1:-1] + un[0:-2,1:-1])))

        v[1:-1,1:-1] = (vn[1:-1,1:-1] -
                        un[1:-1,1:-1] * dt / dx *
                       (vn[1:-1,1:-1] - vn[1:-1,0:-2]) -
                        vn[1:-1,1:-1] * dt / dy *
                       (vn[1:-1,1:-1] - vn[0:-2,1:-1]) -
                        dt / (2 * rho * dy) *
                       (p[2:,1:-1] - p[0:-2,1:-1]) +
                        nu * (dt / dx**2 *
                       (vn[1:-1,2:] - 2 * vn[1:-1,1:-1] + vn[1:-1,0:-2]) +
                             dt / dy**2 *
                       (vn[2:,1:-1] - 2 * vn[1:-1,1:-1] + vn[0:-2,1:-1])))

        # Conditions limites
        u[0,:] = u[-1,:] = u[:,0] = 0
        u[:,-1] = 1  # Paroi supérieure mobile
        v[0,:] = v[-1,:] = v[:,0] = v[:,-1] = 0

    # Tracé
    X, Y = np.meshgrid(np.linspace(0, 2, nx), np.linspace(0, 2, ny))
    fig, ax = plt.subplots(figsize=(8,6))
    ax.quiver(X, Y, u, v)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_title("Champ de vitesse (Navier-Stokes)")
    st.pyplot(fig)