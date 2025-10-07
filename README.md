# 🌌 FCC Metal Simulation using Lennard–Jones Potentials

This project explores how simple pair potentials can capture the complex behaviour of real metals.  
Inspired by **Kanhaiya et al. (2021)** - the goal was to reproduce equilibrium lattice constants for ten FCC metals using **Lennard–Jones interatomic potentials**.

---

## 🔬 Overview

Interatomic potentials describe how atoms interact, bind, and move within a crystal.  
Here, both the **12–6** and **9–6** Lennard–Jones forms were implemented in Python to simulate the structural stability of metallic lattices.

The **12–6 potential** was first tested on Argon as a benchmark, validating the Runge–Kutta–style integration and optimization procedure.  
The **9–6 potential** was then applied to real FCC metals — including Ir, Fe, Rh, Ca, Sr, and others — using parameters derived from the reference paper.

Each metal’s **equilibrium lattice constant** was obtained by minimizing the total potential energy per atom with respect to lattice spacing via the **Nelder–Mead** algorithm.

---

## ⚙️ Results & Insights

The simulation successfully reproduces the expected lattice parameters with strong agreement to experimental data.  
While the Lennard–Jones model is a simplified approximation, this work shows how it can still capture essential trends in metallic bonding and provide a foundation for exploring more advanced many-body potentials.

---

Developed as part of **Green Project 2** under the supervision of **Alessandro Lunghi**  
🎓 *MSc Quantum Science and Technology, Trinity College Dublin*  
💡 *Computational physics • Materials simulation • Atomistic modelling*
