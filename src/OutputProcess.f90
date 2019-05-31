SUBROUTINE OutputProcess
  USE DEF_MOD
  USE Config_MOD
  
  IMPLICIT NONE
  
  INTEGER :: outgoing_idx, phi_idx, ig
  CHARACTER*100 :: filename_out = " "
  
  IF (Problem .EQ. 1) THEN
    Response = Response / Nrun
    filename_out = TRIM(filename) // ".ss"
    OPEN(outfid_SS, file = filename_out)
    DO outgoing_idx = 1, Ng*Nazi*16
      WRITE(outfid_SS, '(<Ng*Nazi*4>ES16.5)'), Response(outgoing_idx,:)
    END DO
    CLOSE(outfid_SS)
  ELSE IF (Problem .EQ. 2) THEN
    Response = Response / Nrun
    filename_out = TRIM(filename) // ".sphi"
    OPEN(outfid_Sphi, file = filename_out)
    DO phi_idx = 1, Ng*3
      WRITE(outfid_Sphi, '(<Ng*Nazi*4>ES16.5)'), Response(phi_idx,:)
    END DO
    CLOSE(outfid_Sphi) 
  ELSE IF (Problem .EQ. 3) THEN
    Response = Response / Nrun; Response2 = Response2 / Nrun
    filename_out = TRIM(filename) // ".vs"
    OPEN(outfid_VS, file = filename_out)
    DO outgoing_idx = 1, Ng*Nazi*16
      WRITE(outfid_VS, '(<Ng*3>ES16.5)'), Response(outgoing_idx,:)
    END DO
    CLOSE(outfid_VS)
    
    filename_out = TRIM(filename) // ".vphi"
    OPEN(outfid_Vphi, file = filename_out)
    DO ig = 1, Ng*3
      phi_idx = ig
      WRITE(outfid_Vphi, '(<Ng*3>ES16.5)'), Response2(phi_idx,:)
    END DO
    CLOSE(outfid_Vphi)
  ELSE
    PRINT *, "Define the type of problem."
  END IF
  
  END SUBROUTINE
  
    