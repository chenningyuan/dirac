#!/bin/sh

./create_files.py

MOL=co; export MOL

for file in response*inp; do
    MENU=${file%.inp}; export MENU
    ../run.sh -incmo -outcmo -get PAMXVC
    mv PAMXVC xvc${MENU#response}
done

MENU=visual_ee; export MENU
../run.sh -incmo -copy xvc_ex -copy xvc_ey -copy xvc_ez -get PLOT_3D; mv PLOT_3D $MENU.3d
exit 0
MENU=visual_1b; export MENU
../run.sh -incmo -copy xvc_bx -copy xvc_by -copy xvc_bz -get PLOT_3D; mv PLOT_3D $MENU.3d
MENU=visual_b1; export MENU
../run.sh -incmo -copy xvc_1x -copy xvc_1y -copy xvc_1z -get PLOT_3D; mv PLOT_3D $MENU.3d
MENU=visual_21; export MENU
../run.sh -incmo -copy xvc_1x -copy xvc_1y -copy xvc_1z -get PLOT_3D; mv PLOT_3D $MENU.3d
MENU=visual_12; export MENU
../run.sh -incmo -copy xvc_2x -copy xvc_2y -copy xvc_2z -get PLOT_3D; mv PLOT_3D $MENU.3d

exit 0
