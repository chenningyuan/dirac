

12.8
----

- Bugfix for screening in MP2. Now you can use
  screening in conjunction with batched processing (for large systems) by
  specifying .IJTSK in input.


13.0
----

- pam can handle full paths of files which are put/get to/from scratch directory,
for example, pam ... --put  "/tmp/milias/testrun/DFCOEF_MI=DFCOEF  /tmp/milias/testrun/MRCONEE  /tmp/milias/testrun/MDCINT_MI=MDCINT"  --get "/tmp/milias/testrun/DFCOEF_MI=DFCOEF_ab"
