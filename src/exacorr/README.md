#General documentation for exacorr module

Remarks:
- The exacorr MOLTRA routine is bugged. Please do **not** use!

To Do:
1. Fix MOLTRA bug
2. Finish ReSpect interface
3. Implement MO to AO transformation for DMs
4. Implement EOM

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