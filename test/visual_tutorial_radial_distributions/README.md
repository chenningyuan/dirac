this tests the "Plotting radial distributions" tutorial on the DIRAC website (https://www.diracprogram.org/doc/release-25/tutorials/visual/radial_distributions/tutorial.html)

Note:
* in the "old" VISUAL code, the final result was a sum of integrated values corresponding to each response vector;
* in the "new" VISUAL code (this test), you need to sum these values, for example, in the `velocity_Mg.out` output:

```
Final grid_function values (grid_function%id=  5) integrated on grid (grid#id=  2):
   1.00000133E+000  -4.10394029E-003  -4.10394103E-003
Final grid_function values (grid_function%id=  6) integrated on grid (grid#id=  2):
  -1.33060120E-006   5.62397734E-007   5.62397794E-007
```

where the `grid_function%id=  5` and the `grid_function%id=  6` correspond to the 'e-e' and 'e-p' response vectors, respectively, and the grid functions were integrated over the DFT grid;

(2) the final result, which includes contributions from both 'e-e' and 'e-p' is therefore:

result = (1.00000133E+000 - 1.33060120E-006, -4.10394029E-003 + 5.62397734E-007, -4.10394103E-003 + 5.62397794E-007)

