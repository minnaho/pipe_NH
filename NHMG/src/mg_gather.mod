	  7O     k820309    h          14.0        Ü­]                                                                                                           
       mg_gather.f90 MG_GATHER                                                    
                                                          
                                                          
                                                          
                                                          
                                                              u #GATHER_2D    #GATHER_3D                      @                               '             .      #NX 	   #NY 
   #NZ    #NPX    #NPY    #INCX    #INCY    #NG2D    #NG    #NGP    #GATHER    #NGX    #NGY    #KEY    #LOCALCOMM    #NEIGHB    #CA    #P    #B    #R    #PX    #PY    #PZ    #DX     #DY !   #DXU "   #DYV #   #DZ $   #DZW %   #ARX &   #ARY '   #ARZ (   #BETA )   #GAMU *   #GAMV +   #ZXDY ,   #ZYDX -   #ALPHA .   #DUMMY3DNZ /   #DUMMY3DNZP 0   #DUMMY3 1   #GATHERBUFFER2D 2   #GATHERBUFFER 3   #DU 4   #DV 5   #DW 6                                              	                                                              
                                                                                                                                                                                                                                                                                                                                                                                                                                                        	                                                      $       
                                                      (                                                             ,                                                             0                                                             4                                                             8                                                                    <                   p          p            p                                                                            `                
            &                   &                   &                   &                                                                   È              y
                                                                                                  ð                
            &                   &                   &                                                                   È              y
                                                                                                  h               
            &                   &                   &                                                                   È              y
                                                                                                  à               
            &                   &                   &                                                                   È              y
                                                                                                  X               
            &                   &                   &                                                                   È              y
                                                                                                  Ð               
            &                   &                   &                                                                   È              y
                                                                                                  H               
            &                   &                   &                                                                   È              y
                                                                                                   À               
            &                   &                                                                   È              y
                                                                                      !                            
            &                   &                                                                   È              y
                                                                                      "                           
            &                   &                                                                   È              y
                                                                                      #            à               
            &                   &                                                                   È              y
                                                                                      $            @               
            &                   &                   &                                                                   È              y
                                                                                      %            ¸               
            &                   &                   &                                                                   È              y
                                                                                      &            0               
            &                   &                   &                                                                   È              y
                                                                                      '            ¨               
            &                   &                   &                                                                   È              y
                                                                                      (                             
            &                   &                                                                   È              y
                                                                                      )                         !  
            &                   &                                                                   È              y
                                                                                      *            à             "  
            &                   &                                                                   È              y
                                                                                      +            @             #  
            &                   &                                                                   È              y
                                                                                      ,                          $  
            &                   &                   &                                                                   È              y
                                                                                      -            	             %  
            &                   &                   &                                                                   È              y
                                                                                      .            	             &  
            &                   &                   &                                                                   È              y
                                                                                      /            
             '  
            &                   &                   &                                                                   È              y
                                                                                      0            
             (  
            &                   &                   &                                                                   È              y
                                                                                      1            ø
             )  
            &                   &                   &                                                                   È              y
                                                                                      2            p             *  
            &                   &                   &                   &                                                                   È              y
                                                                                      3                          +  
            &                   &                   &                   &                   &                                                                   È              y
                                                                                      4            ¨             ,  
            &                   &                   &                                                                   È              y
                                                                                      5                          -  
            &                   &                   &                                                                   È              y
                                                                                      6                         .  
            &                   &                   &                                                                   È              y
                                                                                           7                                                      4                                            8                                                      8          @ P                               9                                   &                                           #GRID_TYPE                                                 :                                
        L            1275070495%         @                                ;                           #LHS <   #RHS >             
                                  <                   #MPI_GROUP =             
                                  >                   #MPI_GROUP =   %         @                                ?                           #LHS @   #RHS B             
                                  @                   #MPI_INFO A             
                                  B                   #MPI_INFO A   %         @                                C                           #LHS D   #RHS F             
                                  D                   #MPI_MESSAGE E             
                                  F                   #MPI_MESSAGE E   %         @                                G                           #LHS H   #RHS J             
                                  H                   #MPI_WIN I             
                                  J                   #MPI_WIN I   %         @                                K                           #LHS L   #RHS N             
                                  L                   #MPI_DATATYPE M             
                                  N                   #MPI_DATATYPE M   %         @                                O                           #LHS P   #RHS R             
                                  P                   #MPI_REQUEST Q             
                                  R                   #MPI_REQUEST Q   %         @                                S                           #LHS T   #RHS V             
                                  T                   #MPI_FILE U             
                                  V                   #MPI_FILE U   %         @                                W                           #LHS X   #RHS Z             
                                  X                   #MPI_ERRHANDLER Y             
                                  Z                   #MPI_ERRHANDLER Y   %         @                                [                           #LHS \   #RHS ^             
                                  \                   #MPI_COMM ]             
                                  ^                   #MPI_COMM ]   %         @                                _                           #LHS `   #RHS b             
                                  `                   #MPI_OP a             
                                  b                   #MPI_OP a   %         @                                c                           #LHS d   #RHS e             
                                  d                   #MPI_GROUP =             
                                  e                   #MPI_GROUP =   %         @                                f                           #LHS g   #RHS h             
                                  g                   #MPI_INFO A             
                                  h                   #MPI_INFO A   %         @                                i                           #LHS j   #RHS k             
                                  j                   #MPI_MESSAGE E             
                                  k                   #MPI_MESSAGE E   %         @                                l                           #LHS m   #RHS n             
                                  m                   #MPI_WIN I             
                                  n                   #MPI_WIN I   %         @                                o                           #LHS p   #RHS q             
                                  p                   #MPI_DATATYPE M             
                                  q                   #MPI_DATATYPE M   %         @                                r                           #LHS s   #RHS t             
                                  s                   #MPI_REQUEST Q             
                                  t                   #MPI_REQUEST Q   %         @                                u                           #LHS v   #RHS w             
                                  v                   #MPI_FILE U             
                                  w                   #MPI_FILE U   %         @                                x                           #LHS y   #RHS z             
                                  y                   #MPI_ERRHANDLER Y             
                                  z                   #MPI_ERRHANDLER Y   %         @                                {                           #LHS |   #RHS }             
                                  |                   #MPI_COMM ]             
                                  }                   #MPI_COMM ]   %         @                                ~                           #LHS    #RHS              
                                                     #MPI_OP a             
                                                     #MPI_OP a   #         @      X                                                 #LEV    #X    #Y              
                                                     
@ P                                                 
              &                   &                                                    D                                                  
               &                   &                                           #         @      X                                                 #LEV    #X    #Y              
                                                     
@ P                                                 
              &                   &                   &                                                    D                                                  
               &                   &                   &                                           #         @                                                      #SPLIT%MOD    #LEV    #X    #Y                                                    MOD           
                                                     
                                                   
              &                   &                   &                                                    D                                                  
               &                   &                   &                                                            @                           a     '                    #MPI_VAL                                                                                    @                           ]     '                    #MPI_VAL                                                                                    @                           Y     '                    #MPI_VAL                                                                                    @                           U     '                    #MPI_VAL                                                                                    @                           Q     '                    #MPI_VAL                                                                                    @                           M     '                    #MPI_VAL                                                                                    @                           I     '                    #MPI_VAL                                                                                    @                           E     '                    #MPI_VAL                                                                                    @                           A     '                    #MPI_VAL                                                                                    @                           =     '                    #MPI_VAL                                                                                 fn#fn    À   @   J   MG_CST       @   J   MG_MPI    @  @   J   MG_TICTOC      @   J   MG_NAMELIST    À  @   J   MG_GRIDS       ^       gen@GATHER #   ^        GRID_TYPE+MG_GRIDS &   r  H   a   GRID_TYPE%NX+MG_GRIDS &   º  H   a   GRID_TYPE%NY+MG_GRIDS &     H   a   GRID_TYPE%NZ+MG_GRIDS '   J  H   a   GRID_TYPE%NPX+MG_GRIDS '     H   a   GRID_TYPE%NPY+MG_GRIDS (   Ú  H   a   GRID_TYPE%INCX+MG_GRIDS (   "  H   a   GRID_TYPE%INCY+MG_GRIDS (   j  H   a   GRID_TYPE%NG2D+MG_GRIDS &   ²  H   a   GRID_TYPE%NG+MG_GRIDS '   ú  H   a   GRID_TYPE%NGP+MG_GRIDS *   B  H   a   GRID_TYPE%GATHER+MG_GRIDS '     H   a   GRID_TYPE%NGX+MG_GRIDS '   Ò  H   a   GRID_TYPE%NGY+MG_GRIDS '     H   a   GRID_TYPE%KEY+MG_GRIDS -   b  H   a   GRID_TYPE%LOCALCOMM+MG_GRIDS *   ª     a   GRID_TYPE%NEIGHB+MG_GRIDS &   F	  <  a   GRID_TYPE%CA+MG_GRIDS %   
  $  a   GRID_TYPE%P+MG_GRIDS %   ¦  $  a   GRID_TYPE%B+MG_GRIDS %   Ê  $  a   GRID_TYPE%R+MG_GRIDS &   î  $  a   GRID_TYPE%PX+MG_GRIDS &     $  a   GRID_TYPE%PY+MG_GRIDS &   6  $  a   GRID_TYPE%PZ+MG_GRIDS &   Z    a   GRID_TYPE%DX+MG_GRIDS &   f    a   GRID_TYPE%DY+MG_GRIDS '   r    a   GRID_TYPE%DXU+MG_GRIDS '   ~    a   GRID_TYPE%DYV+MG_GRIDS &     $  a   GRID_TYPE%DZ+MG_GRIDS '   ®  $  a   GRID_TYPE%DZW+MG_GRIDS '   Ò  $  a   GRID_TYPE%ARX+MG_GRIDS '   ö  $  a   GRID_TYPE%ARY+MG_GRIDS '       a   GRID_TYPE%ARZ+MG_GRIDS (   &    a   GRID_TYPE%BETA+MG_GRIDS (   2    a   GRID_TYPE%GAMU+MG_GRIDS (   >    a   GRID_TYPE%GAMV+MG_GRIDS (   J  $  a   GRID_TYPE%ZXDY+MG_GRIDS (   n  $  a   GRID_TYPE%ZYDX+MG_GRIDS )      $  a   GRID_TYPE%ALPHA+MG_GRIDS -   ¶!  $  a   GRID_TYPE%DUMMY3DNZ+MG_GRIDS .   Ú"  $  a   GRID_TYPE%DUMMY3DNZP+MG_GRIDS *   þ#  $  a   GRID_TYPE%DUMMY3+MG_GRIDS 2   "%  <  a   GRID_TYPE%GATHERBUFFER2D+MG_GRIDS 0   ^&  T  a   GRID_TYPE%GATHERBUFFER+MG_GRIDS &   ²'  $  a   GRID_TYPE%DU+MG_GRIDS &   Ö(  $  a   GRID_TYPE%DV+MG_GRIDS &   ú)  $  a   GRID_TYPE%DW+MG_GRIDS    +  q       IP+MG_CST    +  q       RP+MG_CST     ,         GRID+MG_GRIDS 3   ,  z       MPI_DOUBLE_PRECISION+MPI_CONSTANTS &   -  b       GROUPEQ+MPI_CONSTANTS *   w-  W   a   GROUPEQ%LHS+MPI_CONSTANTS *   Î-  W   a   GROUPEQ%RHS+MPI_CONSTANTS %   %.  b       INFOEQ+MPI_CONSTANTS )   .  V   a   INFOEQ%LHS+MPI_CONSTANTS )   Ý.  V   a   INFOEQ%RHS+MPI_CONSTANTS (   3/  b       MESSAGEEQ+MPI_CONSTANTS ,   /  Y   a   MESSAGEEQ%LHS+MPI_CONSTANTS ,   î/  Y   a   MESSAGEEQ%RHS+MPI_CONSTANTS $   G0  b       WINEQ+MPI_CONSTANTS (   ©0  U   a   WINEQ%LHS+MPI_CONSTANTS (   þ0  U   a   WINEQ%RHS+MPI_CONSTANTS )   S1  b       DATATYPEEQ+MPI_CONSTANTS -   µ1  Z   a   DATATYPEEQ%LHS+MPI_CONSTANTS -   2  Z   a   DATATYPEEQ%RHS+MPI_CONSTANTS (   i2  b       REQUESTEQ+MPI_CONSTANTS ,   Ë2  Y   a   REQUESTEQ%LHS+MPI_CONSTANTS ,   $3  Y   a   REQUESTEQ%RHS+MPI_CONSTANTS %   }3  b       FILEEQ+MPI_CONSTANTS )   ß3  V   a   FILEEQ%LHS+MPI_CONSTANTS )   54  V   a   FILEEQ%RHS+MPI_CONSTANTS +   4  b       ERRHANDLEREQ+MPI_CONSTANTS /   í4  \   a   ERRHANDLEREQ%LHS+MPI_CONSTANTS /   I5  \   a   ERRHANDLEREQ%RHS+MPI_CONSTANTS %   ¥5  b       COMMEQ+MPI_CONSTANTS )   6  V   a   COMMEQ%LHS+MPI_CONSTANTS )   ]6  V   a   COMMEQ%RHS+MPI_CONSTANTS #   ³6  b       OPEQ+MPI_CONSTANTS '   7  T   a   OPEQ%LHS+MPI_CONSTANTS '   i7  T   a   OPEQ%RHS+MPI_CONSTANTS '   ½7  b       GROUPNEQ+MPI_CONSTANTS +   8  W   a   GROUPNEQ%LHS+MPI_CONSTANTS +   v8  W   a   GROUPNEQ%RHS+MPI_CONSTANTS &   Í8  b       INFONEQ+MPI_CONSTANTS *   /9  V   a   INFONEQ%LHS+MPI_CONSTANTS *   9  V   a   INFONEQ%RHS+MPI_CONSTANTS )   Û9  b       MESSAGENEQ+MPI_CONSTANTS -   =:  Y   a   MESSAGENEQ%LHS+MPI_CONSTANTS -   :  Y   a   MESSAGENEQ%RHS+MPI_CONSTANTS %   ï:  b       WINNEQ+MPI_CONSTANTS )   Q;  U   a   WINNEQ%LHS+MPI_CONSTANTS )   ¦;  U   a   WINNEQ%RHS+MPI_CONSTANTS *   û;  b       DATATYPENEQ+MPI_CONSTANTS .   ]<  Z   a   DATATYPENEQ%LHS+MPI_CONSTANTS .   ·<  Z   a   DATATYPENEQ%RHS+MPI_CONSTANTS )   =  b       REQUESTNEQ+MPI_CONSTANTS -   s=  Y   a   REQUESTNEQ%LHS+MPI_CONSTANTS -   Ì=  Y   a   REQUESTNEQ%RHS+MPI_CONSTANTS &   %>  b       FILENEQ+MPI_CONSTANTS *   >  V   a   FILENEQ%LHS+MPI_CONSTANTS *   Ý>  V   a   FILENEQ%RHS+MPI_CONSTANTS ,   3?  b       ERRHANDLERNEQ+MPI_CONSTANTS 0   ?  \   a   ERRHANDLERNEQ%LHS+MPI_CONSTANTS 0   ñ?  \   a   ERRHANDLERNEQ%RHS+MPI_CONSTANTS &   M@  b       COMMNEQ+MPI_CONSTANTS *   ¯@  V   a   COMMNEQ%LHS+MPI_CONSTANTS *   A  V   a   COMMNEQ%RHS+MPI_CONSTANTS $   [A  b       OPNEQ+MPI_CONSTANTS (   ½A  T   a   OPNEQ%LHS+MPI_CONSTANTS (   B  T   a   OPNEQ%RHS+MPI_CONSTANTS    eB  _       GATHER_2D    ÄB  @   a   GATHER_2D%LEV    C  ¤   a   GATHER_2D%X    ¨C  ¤   a   GATHER_2D%Y    LD  _       GATHER_3D    «D  @   a   GATHER_3D%LEV    ëD  ¼   a   GATHER_3D%X    §E  ¼   a   GATHER_3D%Y    cF  n       SPLIT    ÑF  <      SPLIT%MOD    G  @   a   SPLIT%LEV    MG  ¼   a   SPLIT%X    	H  ¼   a   SPLIT%Y %   ÅH  ]       MPI_OP+MPI_CONSTANTS -   "I  H   a   MPI_OP%MPI_VAL+MPI_CONSTANTS '   jI  ]       MPI_COMM+MPI_CONSTANTS /   ÇI  H   a   MPI_COMM%MPI_VAL+MPI_CONSTANTS -   J  ]       MPI_ERRHANDLER+MPI_CONSTANTS 5   lJ  H   a   MPI_ERRHANDLER%MPI_VAL+MPI_CONSTANTS '   ´J  ]       MPI_FILE+MPI_CONSTANTS /   K  H   a   MPI_FILE%MPI_VAL+MPI_CONSTANTS *   YK  ]       MPI_REQUEST+MPI_CONSTANTS 2   ¶K  H   a   MPI_REQUEST%MPI_VAL+MPI_CONSTANTS +   þK  ]       MPI_DATATYPE+MPI_CONSTANTS 3   [L  H   a   MPI_DATATYPE%MPI_VAL+MPI_CONSTANTS &   £L  ]       MPI_WIN+MPI_CONSTANTS .    M  H   a   MPI_WIN%MPI_VAL+MPI_CONSTANTS *   HM  ]       MPI_MESSAGE+MPI_CONSTANTS 2   ¥M  H   a   MPI_MESSAGE%MPI_VAL+MPI_CONSTANTS '   íM  ]       MPI_INFO+MPI_CONSTANTS /   JN  H   a   MPI_INFO%MPI_VAL+MPI_CONSTANTS (   N  ]       MPI_GROUP+MPI_CONSTANTS 0   ïN  H   a   MPI_GROUP%MPI_VAL+MPI_CONSTANTS 