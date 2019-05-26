General documentation for exacorr module

Remarks:
- To use the standalone version of exacorr use the exacorr.x executable i.c.w. exacorr.inp (only need to specify EXACC section)

To Do:
- [ ] Finish ReSpect interface (@stanpa)
- [ ] Implement EOM (@stanpa)
  - [ ] Implement Davidson for EOM (@gomes.asp)
- [ ] Create module for parameters

<table>
  <tr>
    <th>Input parameters</th><th>Comment</th>
  </tr>
  <tr>
  	<td>.OCCUPIED</td><td>As list of MOs or energy range</td>
  </tr>
  <tr>
  	<td>.VIRTUAL</td><td>As list of MOs or energy range</td>
  </tr>
  <tr>
  	<td>.EXATENSOR</td><td>Calculation with ExaTensor</td>
  </tr>
  <tr>
  	<td>.CCDOUBLES</td><td>CCD calculation</td>
  </tr>
  <tr>
  	<td>.TCONVERG</td><td>Set convergence criterria</td>
  </tr>
  <tr>
  	<td>.LAMBDA</td><td>Solve Λ-equations</td>
  </tr>
  <tr>
  	<td>.NOINTS</td><td>Debugging. Set intergrals to 0</td>
  </tr>
  <tr>
  	<td>.EXA_BLOCKSIZE</td><td>Number to branch the tensors</td>
  </tr>
  <tr>
  	<td>.MOINT_SCHEME</td><td>Only 1 working implementation now</td>
  </tr>
  <tr>
  	<td> .NCYCLES</td><td>Set the number of CC cycles</td>
  </tr>
</table>