PROGRAM MAIN
  USE DEF_MOD, ONLY : Nrun, problem
  USE ReadInput_MOD
  USE NEUTRON_MOD
  USE SOLVER_MOD
  IMPLICIT NONE
  REAL :: t1, t2, delt
  INTEGER :: size, rank, ierr, i, outgoing_idx
  
  CALL CPU_TIME(t1)

  CALL Scaninput
  CALL Readinput
  CALL SetGeom
  CALL Alloc
  DO i = 1, Nrun
    CALL GenNeutrons
    print *, "i = ", i
    CALL Solver
  END DO
  CALL OutputProcess
  
  CALL CPU_TIME(t2)
  delt = t2-t1
  PRINT *, "Total time : ", t2, " - ", t1, " = ", delt, "s"
  WRITE(outfid,*), "Total time : ", t2, " - ", t1, " = ", delt, "s"
  
  
  END PROGRAM