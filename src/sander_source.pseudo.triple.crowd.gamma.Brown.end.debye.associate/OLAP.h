!Created by Antonios Samiotakis 9.11.2010 to include the 
!data structures needed to read the overlap info.
!All these were originally in LAN.h (amber6 version)

integer iolap,Q,count_olap,i_olap,j_olap,dummy_olap,ii_olap,jj_olap,M_olap
_REAL_ xdist_olap(350000),xolapcut,xener,xoverlap,temp_olap,xx_olap
_REAL_ xdummy_olap,dist


COMMON/xolap/xdist_olap,xolapcut,xener,xoverlap
COMMON/olap/iolap

