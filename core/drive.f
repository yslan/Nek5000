      program NEKTON

      USE parallel_mod, ONLY : parallel_init => init
      USE scrns_mod, ONLY : scrns_init => init
      USE ctmp0_mod, ONLY : ctmp0_init => init
      USE ctmp1_mod, ONLY : ctmp1_init => init
      USE soln_mod, ONLY : soln_init => init
      USE hsmg_mod, ONLY : hsmg_init => init
      USE mvgeom_mod, ONLY : mvgeom_init => init
      USE vproj_mod, ONLY : vproj_init => init
      USE vrthov_mod, ONLY : vrthov_init => init
      USE orthot_mod, ONLY : orthot_init => init
      USE orthov_mod, ONLY : orthov_init => init
      USE screv_mod, ONLY : screv_init => init
      USE scrcg_mod, ONLY : scrcg_init => init
      USE scrch_mod, ONLY : scrch_init => init
      USE scrmg_mod, ONLY : scrmg_init => init
      USE scrsf_mod, ONLY : scrsf_init => init
      USE scruz_mod, ONLY : scruz_init => init
      USE scrvh_mod, ONLY : scrvh_init => init
      USE orthostrs_mod, ONLY : orthostrs_init => init
      USE orthop_mod, ONLY : orthop_init => init
      USE noncon_mod, ONLY : noncon_init => init
      USE neknek_mod, ONLY : neknek_init => init
      USE mass_mod, ONLY : mass_init => init
      USE gmres_mod, ONLY : gmres_init => init
      USE geom_mod, ONLY : geom_init => init
      USE dealias_mod, ONLY : dealias_init => init
      USE cvode_mod, ONLY : cvode_init => init
      USE avg_mod, ONLY : avg_init => init
      USE adjoint_mod, ONLY : adjoint_init => init
      USE cbplan_vol_ms_mod, ONLY : cbplan_vol_ms_init => init
      USE cvflow_a_mod, ONLY : cvflow_a_init => init
      USE scrdg_mod, ONLY : scrdg_init => init
      USE outtmp_mod, ONLY : outtmp_init => init
      USE orthox_mod, ONLY : orthox_init => init
      USE cvflow_nn_mod, ONLY : cvflow_nn_init => init
      USE c_is1_mod, ONLY : c_is1_init => init
      USE fastd_mod, ONLY : fastd_init => init
      USE scrxxti_mod, ONLY : scrxxti_init => init
      USE fdmh1_mod, ONLY : fdmh1_init => init
      USE scrpre_mod, ONLY : scrpre_init => init
      USE swaplengths_mod, ONLY : swaplengths_init => init
      USE weightop_mod, ONLY : weightop_init => init
      USE scrhi_mod, ONLY : scrhi_init => init
      USE input_mod, ONLY : input_init => init
      USE topol_mod, ONLY : topol_init => init
      USE scrct_mod, ONLY : scrct_init => init
      USE ctmpf_mod, ONLY : ctmpf_init => init
      USE scrxxt_mod, ONLY : scrxxt_init => init
      USE scrpr2_mod, ONLY : scrpr2_init => init
      USE fastg_mod, ONLY : fastg_init => init
      USE ivrtx_mod, ONLY : ivrtx_init => init
      USE domain_mod, ONLY : domain_init => init
      include 'mpif.h'

      integer comm
      comm = MPI_COMM_WORLD

      call parallel_init()
      call scrns_init()
      call ctmp0_init()
      call ctmp1_init()
      call soln_init()
      call hsmg_init()
      call mvgeom_init()
      call vproj_init()
      call vrthov_init()
      call orthot_init()
      call orthov_init()
      call screv_init()
      call scrcg_init()
      call scrch_init()
      call scrmg_init()
      call scrsf_init()
      call scruz_init()
      call scrvh_init()
      call orthostrs_init()
      call orthop_init()
      call noncon_init()
      call neknek_init()
      call mass_init()
      call gmres_init()
      call geom_init()
      call dealias_init()
      call cvode_init()
      call avg_init()
      call adjoint_init()
      call cbplan_vol_ms_init()
      call cvflow_a_init()
      call scrdg_init()
      call outtmp_init()
      call orthox_init()
      call cvflow_nn_init()
      call c_is1_init()
      call fastd_init()
      call scrxxti_init()
      call fdmh1_init()
      call scrpre_init()
      call swaplengths_init()
      call weightop_init()
      call scrhi_init()
      call input_init()
      call topol_init()
      call scrct_init()
      call ctmpf_init()
      call scrxxt_init()
      call scrpr2_init()
      call fastg_init()
      call ivrtx_init()
      call domain_init()

      call nek_init(comm)
      call nek_solve()
      call nek_end()

      end
