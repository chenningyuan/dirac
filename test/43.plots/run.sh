#!/bin/sh

MOL=co; export MOL

for file in response*inp; do
    MENU=${file%.inp}; export MENU
    ../run.sh -incmo -outcmo -get PAMXVC
    mv PAMXVC xvc${MENU#response}
done

MENU=visual_ee; export MENU
../run.sh -incmo -copy "xvc_ex xvc_ey xvc_ez" -get PLOT_3D; mv PLOT_3D $MENU.3d
MENU=visual_1b; export MENU
../run.sh -incmo -copy "xvc_bx xvc_by xvc_bz" -get PLOT_3D; mv PLOT_3D $MENU.3d
MENU=visual_b1; export MENU
../run.sh -incmo -copy "xvc_1x xvc_1y xvc_1z" -get PLOT_3D; mv PLOT_3D $MENU.3d
MENU=visual_21; export MENU
../run.sh -incmo -copy "xvc_1x xvc_1y xvc_1z" -get PLOT_3D; mv PLOT_3D $MENU.3d
MENU=visual_12; export MENU
../run.sh -incmo -copy "xvc_2x xvc_2y xvc_2z" -get PLOT_3D; mv PLOT_3D $MENU.3d

exit 0
