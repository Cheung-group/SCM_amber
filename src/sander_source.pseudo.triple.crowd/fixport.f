!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C  PUPIL program ITR-MEDIUM PROJECT                       C
!C                                                         C
!C  Author: 	J. Torras                                  C
!C  Version:	0.0.8.6                                    C
!C  Date:	09-08-2005                                     C
!C                                                         C
!C  torras@qtp.ufl.edu                                     C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE FIXPORT
        CHARACTER*30 ARGV,HOST,PORT,OPTPRINT
		INTEGER IOPTPRINT
        INTEGER I,IARGC, M
        M = IARGC()
        DO I = 1, M
             CALL GETARG ( I, ARGV )
             IF(ARGV .EQ. "-ORBInitialPort") THEN
           		CALL GETARG(I+1,PORT)
           		WRITE(*,*)'FOUND PORT ',PORT
             ELSE IF(ARGV .EQ. '-ORBInitialHost') THEN
           		CALL GETARG(I+1,HOST)
           		WRITE(*,*)'FOUND HOST ',HOST
             ELSE IF(ARGV .EQ. '-OptPrint') THEN
           		CALL GETARG(I+1,OPTPRINT)
           		READ(OPTPRINT,'(I30)') IOPTPRINT
           		WRITE(*,*)'FOUND LOGLEVEL ',IOPTPRINT
             END IF
        END DO
!C	ASSING HOST, PORT AND PRINTLOG TO PUPIL SYSTEM
		CALL SETCORBANAMESERVER(HOST,PORT,IOPTPRINT)             
!C       WRITE( *, '( I2, 1X, A )' ) I, ARGV
        RETURN
        END
        
