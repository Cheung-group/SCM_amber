!------------------------------------------------------------------
!Ross Walker (TSRI, 2005)
!This header file contains definitions that can be used for 
!offsets into arrays in place of simple integers. This makes
!it easier to ensure you are accessing the correct element
!out of the data array.
!------------------------------------------------------------------

!QMMM RIJ OFFSETS
#define QMMMONERIJ        1 
#define QMMMRIJ           2 
#define QMMMEXP1          3
#define QMMMEXP2          4
#define QMMMSQRTAEE       5
#define QMMMSQRTRRADDADE  6
#define QMMMSQRTRRMDDADE  7
#define QMMMSQRTRR2AQE    8
#define QMMMSQRTRRAQQAQE  9
#define QMMMSQRTRRMQQAQE  10
#define QMMMSQRTRRAQQ2AQE 11

#define QMMMNORIJ         11

!QMQM RIJ OFFSETS
#define QMQMRIJBOHRS2   1
#define QMQMONERIJ      2
#define QMQMRIJ         3
#define QMQMRIJBOHRS    4
#define QMQMEXP1I       5
#define QMQMEXP1J       6
#define QMQMSQRTAEE     7

#define QMQMNORIJ       7

!Note, these 4 PDDG options should go last as QMQMNORIJ will be
!reduced by 4 if PDDG hamiltonian is not in use.
#define QMQMPDDG_EXP1   8
#define QMQMPDDG_EXP2   9
#define QMQMPDDG_EXP3   10
#define QMQMPDDG_EXP4   11

#define QMQMNOPDDG      4


