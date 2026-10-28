#  Satellite Rendezvous Simulation in MATLAB

Simulates a **nonlinear deputy-chaser rendezvous mission** using **Lyapunov-based feedback control**, with 3D Earth-centered visualization and full orbit dynamics. Built in **pure MATLAB**, no toolboxes required.


---

##  Demo Video

 [Watch on YouTube]https://www.youtube.com/watch?v=uv7CveoUqE0


https://github.com/user-attachments/assets/cfa3160f-e7e4-4069-8d7c-42047dfd8262


---

##  Features

-  Nonlinear relative orbital dynamics (with gravity gradient)
-  Lyapunov control law: `u = -Kr·δr - Δa - P·δv`
-  Control saturation (e.g. 1 mm/s²)
-  Animated 3D orbits with `plot_Earth`
-  Time-history plots: position, velocity, control input
-  Computes control effort: ∫‖u‖ dt
-  Dark-themed visualizations

##  Files

| File                          | Description                                    |
|-------------------------------|------------------------------------------------|
| `run_relative_orbit_sim.m`    | Main script with editable parameters           |
| `simulateRelativeOrbit.m`     | All-in-one engine with full simulation logic   |
| `plot_Earth.m`                | Earth rendering function                       |
