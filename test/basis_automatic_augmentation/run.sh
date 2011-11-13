#!/bin/sh

export MENU=hf.inp

export MOL=aug-cc-pVDZ.mol
../run.sh
export MOL=d-aug-cc-pVDZ.mol
../run.sh
export MOL=t-aug-cc-pVDZ.mol
../run.sh

export MOL=Turbomole-DZP.mol
../run.sh
export MOL=s-aug-Turbomole-DZP.mol
../run.sh
export MOL=d-aug-Turbomole-DZP.mol
../run.sh

export MOL=Turbomole-TZVPP.mol
../run.sh
export MOL=s-aug-Turbomole-TZVPP.mol
../run.sh

export MOL=dyall.v2z.mol
../run.sh
export MOL=d-aug-dyall.v2z.mol
../run.sh

exit 0
