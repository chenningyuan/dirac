#!/bin/sh

export MENU=hf.inp

export MOL=aug-cc-pVDZ.mol
../SCRIPTS/run.sh
export MOL=d-aug-cc-pVDZ.mol
../SCRIPTS/run.sh
export MOL=t-aug-cc-pVDZ.mol
../SCRIPTS/run.sh

export MOL=Turbomole-DZP.mol
../SCRIPTS/run.sh
export MOL=s-aug-Turbomole-DZP.mol
../SCRIPTS/run.sh
export MOL=d-aug-Turbomole-DZP.mol
../SCRIPTS/run.sh

export MOL=Turbomole-TZVPP.mol
../SCRIPTS/run.sh
export MOL=s-aug-Turbomole-TZVPP.mol
../SCRIPTS/run.sh

export MOL=dyall.v2z.mol
../SCRIPTS/run.sh
export MOL=d-aug-dyall.v2z.mol
../SCRIPTS/run.sh

exit 0
