#Info needed from ReSpect restart file to function as DFCOEF

- Description of calculation (optional)
- number of symmetries (nsym)
- number of basis functions (nao)
- Number of "electronic" (=positive energy) solutions nmo
- Number of "positronic" (=negative energy) solutions npo
(nto = npo + nmo)
- Coefficients, real(8),                            format coeff_comp(nao,nto,spin) ? Or seperate alpha and beta array -> might be useful for open shell
                                                    => need to transform complex to quaternion, format coeff(nao,nto,nz)
- eigenvalues, real(8),                             format eig(nto)
- boson irrep to which spinor belongs, integer,     format ibeig(nto)
- number of shells
- number of exponents
- number of contractions (if none, then ncnt=1)
- orb momentum,                                     format nhkt(nshells)
- atom number,                                      format ncent(nshells)
- number of primitives,                             format nuco(nshells)
- start of exponent values,                         format jstrt(nshells)
- number of contracted functions,                   format numcf(nshells)
- coordinates, real(8),                             format cent(nshells,3)
- values of the exponents, real(8),                 format priexp(npriexp)
- values of the contraction coefficients, real(8),  format priccf(npriexp,ncnt)
