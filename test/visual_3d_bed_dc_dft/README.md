

Test description:

* testing `bedsin`:

  * first, we run response calculations (`bed_anisotropic.inp`), from which we identify:
    * the frequency of the first s->d (B3g) excitation (full interaction operator) is 48.96760356
    * the corresponding transition moments are:

```
Exc. Frequency (au)                                e-e                  e-p                 Total
   1      48.9676035639 real part (sin[kr]) :     -2.2092173501E-08     2.0398571932E-10    -2.1888187781E-08
                        imag part (cos[kr]) :      0.0000000000E+00     0.0000000000E+00     0.0000000000E+00
```

    * we look for the index of response vector in PAMXVC (or RESPONSE.h5) file corresponding to each contribution, and find that
      * e-e contribution comes from response vector #1 (E+) and #7 (E-)
      * e-p contribution comes from response vector #13 (P+) and #19 (P-)

  * in `vis_bedsin_anisotropic.inp`, we then set:
    * the same excitation energy in .BEDFIX (note that it needs to be specified for each 'bedsin' gridfunction)
    * we numerically calculate the transition moments (for the solution vectors 1,7, 13, 19)

    * in the output, we find:

```
  (#1)
  Final grid function values (id=  1) integrated on grid (id=  2):
  -1.32584399E-022

  (#7)
  Final grid function values (id=  2) integrated on grid (id=  2):
  -2.20921733E-008

  (#13)
  Final grid function values (id=  3) integrated on grid (id=  2):
   2.97689285E-027

  (#19)
  Final grid function values (id=  4) integrated on grid (id=  2):
   2.03985719E-010
```

  * so (#7) well reproduces real e-e part and (#19) well reproduces real e-p part

