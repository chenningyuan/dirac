      logical          dofit, lhsfit, rhsfit, robust_e, dobox,          &
     &                 df_constr, df_gaudmp, df_attenu, diatom, lin_df, &
     &                 save_cfit, save_vab, df_fitvxc, using_diffdens,  &
     &                 dodflr

      integer          ndima, ndimb, ndimc, ndimd, nbaseaux,            &
     &                 nodclsnorm, nodclsaux, nodbtnorm, nodbtaux

      real(8)          rob_e_corr, df_attpar, df_gaupar, df_scrthr

      common /densftl/ dofit, lhsfit, rhsfit, robust_e, dobox,          &
     &                 df_constr, df_gaudmp, df_attenu, diatom, lin_df, &
     &                 save_cfit, save_vab, df_fitvxc, using_diffdens,  &
     &                 dodflr

      common /densfti/ ndima, ndimb, ndimc, ndimd, nbaseaux,            &
     &                 nodclsnorm, nodclsaux, nodbtnorm, nodbtaux

      common /densftr/ rob_e_corr, df_attpar, df_gaupar, df_scrthr
