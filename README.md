# PyFemHeat
Python library for solving heat equation using finite element method.

This repository is an home-made project to implement finite element method applied to heat equation. The pdf file summarize the mathematical aspect behind technical implementation. For now, only time-independent heat equation can be solve on simple 2-D geometry. To see a full demo of what has already been implemented, check out the Demo_FEM Jupyter notebook.

Repo structure :
<ul>
  <li> Demo_FEM.ipynb : demo jupyter notebook </li>
  <li> src : source files for matrix calculation and shapes generation. </li>
  <li> utils : utility files to extract property from current simulated geometry (shape functions, quadrature points, ...) </li>
</ul>

This work has been done mainly on my free-time so future improvement may take a while.

To be added:
<ul>
  <li> Time-dependent FEM solver. </li>
  <li> Extend shapes generation capabilities (2-D and 3-D). </li>
  <li> Add Neumann boundary conditions. </li>
</ul>

Contact me at: robin.feuillassier@gmail.com
