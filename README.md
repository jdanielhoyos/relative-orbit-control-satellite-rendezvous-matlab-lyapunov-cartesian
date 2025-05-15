# 🛰️ Satellite Rendezvous Simulation in MATLAB

Simulates a **nonlinear deputy-chaser rendezvous mission** using **Lyapunov-based feedback control**, with 3D Earth-centered visualization and full orbit dynamics. Built in **pure MATLAB**, no toolboxes required.

![orbit-demo](https://yourdomain.com/demo.gif) <!-- Optional: add GIF if available -->

---

## 📽️ Demo Video

👉 [Watch on YouTube]([https://www.youtube.com/watch?v=YOUR_VIDEO_LINK](https://www.youtube.com/watch?v=uv7CveoUqE0))

---

## 📌 Features

- 🚀 Nonlinear relative orbital dynamics (with gravity gradient)
- ⚙️ Lyapunov control law: `u = -Kr·δr - Δa - P·δv`
- 📉 Control saturation (e.g. 1 mm/s²)
- 🌍 Animated 3D orbits with `plot_Earth`
- 📊 Time-history plots: position, velocity, control input
- 🧮 Computes control effort: ∫‖u‖ dt
- 🌑 Dark-themed visualizations
