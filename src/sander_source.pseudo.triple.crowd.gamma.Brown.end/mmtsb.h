!+ MMTSB Replica Exchange


! public members

integer mmtsb_switch       ! Replica Exchange Control:
integer mmtsb_off          ! No Replica Exchange
parameter ( mmtsb_off = 0 )
integer mmtsb_temp_rex     ! Temperature Replica Exchange
parameter ( mmtsb_temp_rex = 1 )
integer mmtsb_lambda_rex   ! Lambda Replica Exchange
parameter ( mmtsb_lambda_rex = 2 )
integer mmtsb_iterations   ! Replica Exchange Frequency in Iterations
logical mmtsb_is_exchanged ! replica exchange occurred this step

common /mmtsb_public/ mmtsb_switch, mmtsb_iterations, mmtsb_is_exchanged

! These may cause compilation errors in the subroutines themselves.
!external Mmtsb_init
!external Mmtsb_newtemp
!external Mmtsb_print_banner

#ifdef MMTSB

! private members; these may move to file mmtsb.inc

character(len=80) datadir     ! local directory for intermediate coordinates
character(len=80) jobid       ! id of this client
integer      sendfiles   ! the value 1 requests file send
character(len=80) servername
character(len=80) serverport
character(len=80) serverid

common /mmtsb_private_server/ servername, serverport, serverid, &
      jobid, datadir, sendfiles

#endif

