	  -Z  Á   k820309    h          14.0        !]                                                                                                           
       mg_solvers.f90 MG_SOLVERS                                                    
                                                          
                                                          
                                                          
                                                          
                                                          
                                                          
                                                          
                         @                          	     '             .      #NX 
   #NY    #NZ    #NPX    #NPY    #INCX    #INCY    #NG2D    #NG    #NGP    #GATHER    #NGX    #NGY    #KEY    #LOCALCOMM    #NEIGHB    #CA    #P    #B    #R    #PX    #PY    #PZ     #DX !   #DY "   #DXU #   #DYV $   #DZ %   #DZW &   #ARX '   #ARY (   #ARZ )   #BETA *   #GAMU +   #GAMV ,   #ZXDY -   #ZYDX .   #ALPHA /   #DUMMY3DNZ 0   #DUMMY3DNZP 1   #DUMMY3 2   #GATHERBUFFER2D 3   #GATHERBUFFER 4   #DU 5   #DV 6   #DW 7                                              
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      	                                                      $       
                                                      (                                                             ,                                                             0                                                             4                                                             8                                                                    <                   p          p            p                                                                            `                
            &                   &                   &                   &                                                                   î              y
                                                                                                  ð                
            &                   &                   &                                                                   î              y
                                                                                                  h               
            &                   &                   &                                                                   î              y
                                                                                                  à               
            &                   &                   &                                                                   î              y
                                                                                                  X               
            &                   &                   &                                                                   î              y
                                                                                                  Ð               
            &                   &                   &                                                                   î              y
                                                                                                   H               
            &                   &                   &                                                                   î              y
                                                                                      !            À               
            &                   &                                                                   î              y
                                                                                      "                            
            &                   &                                                                   î              y
                                                                                      #                           
            &                   &                                                                   î              y
                                                                                      $            à               
            &                   &                                                                   î              y
                                                                                      %            @               
            &                   &                   &                                                                   î              y
                                                                                      &            ¸               
            &                   &                   &                                                                   î              y
                                                                                      '            0               
            &                   &                   &                                                                   î              y
                                                                                      (            ¨               
            &                   &                   &                                                                   î              y
                                                                                      )                             
            &                   &                                                                   î              y
                                                                                      *                         !  
            &                   &                                                                   î              y
                                                                                      +            à             "  
            &                   &                                                                   î              y
                                                                                      ,            @             #  
            &                   &                                                                   î              y
                                                                                      -                          $  
            &                   &                   &                                                                   î              y
                                                                                      .            	             %  
            &                   &                   &                                                                   î              y
                                                                                      /            	             &  
            &                   &                   &                                                                   î              y
                                                                                      0            
             '  
            &                   &                   &                                                                   î              y
                                                                                      1            
             (  
            &                   &                   &                                                                   î              y
                                                                                      2            ø
             )  
            &                   &                   &                                                                   î              y
                                                                                      3            p             *  
            &                   &                   &                   &                                                                   î              y
                                                                                      4                          +  
            &                   &                   &                   &                   &                                                                   î              y
                                                                                      5            ¨             ,  
            &                   &                   &                                                                   î              y
                                                                                      6                          -  
            &                   &                   &                                                                   î              y
                                                                                      7                         .  
            &                   &                   &                                                                   î              y
                                                                                           8                                                      4                                            9                                                      8                                            :                                                      8          D  @                              ;                      D                                  <                                                        =                      @                                >                                   &                                           #GRID_TYPE 	   #         @                                  ?                   #TIC%TRIM @   #LEV A   #STRING B                                               @     TRIM           
                                 A                     
                                 B                    1 #         @                                  C                    #LEV D   #SUMLOC E   #SUMGLO F             
                                 D                     
                                 E     
                                                F     
       #         @                                  G                   #COMPUTE_RESIDUAL%SIZE H   #COMPUTE_RESIDUAL%SQRT I   #LEV J   #RES K                                               H     SIZE                                             I     SQRT           
                                 J                                                     K     
                 D                                 L                      D                                 M     
       #         @                                  N                   #TOC%REAL O   #TOC%TRIM P   #LEV Q   #STRING R                                @              O     REAL                                             P     TRIM           
                                 Q                     
                                 R                    1              @                              S            #         @                                  T                    #LEV U             
                                 U           #         @                                  V                   #RELAX%SIZE W   #RELAX%TRIM X   #LEV Y   #NSWEEPS Z                                               W     SIZE                                             X     TRIM           
                                 Y                     
                                 Z                     D  @                              [            #         @                                  \                    #LEV ]             
                                 ]                     D  @                              ^                      D  @                              _            %         @                                `                           #LHS a   #RHS c             
                                  a                   #MPI_GROUP b             
                                  c                   #MPI_GROUP b   %         @                                d                           #LHS e   #RHS g             
                                  e                   #MPI_INFO f             
                                  g                   #MPI_INFO f   %         @                                h                           #LHS i   #RHS k             
                                  i                   #MPI_MESSAGE j             
                                  k                   #MPI_MESSAGE j   %         @                                l                           #LHS m   #RHS o             
                                  m                   #MPI_WIN n             
                                  o                   #MPI_WIN n   %         @                                p                           #LHS q   #RHS s             
                                  q                   #MPI_DATATYPE r             
                                  s                   #MPI_DATATYPE r   %         @                                t                           #LHS u   #RHS w             
                                  u                   #MPI_REQUEST v             
                                  w                   #MPI_REQUEST v   %         @                                x                           #LHS y   #RHS {             
                                  y                   #MPI_FILE z             
                                  {                   #MPI_FILE z   %         @                                |                           #LHS }   #RHS              
                                  }                   #MPI_ERRHANDLER ~             
                                                     #MPI_ERRHANDLER ~   %         @                                                           #LHS    #RHS              
                                                     #MPI_COMM              
                                                     #MPI_COMM    %         @                                                           #LHS    #RHS              
                                                     #MPI_OP              
                                                     #MPI_OP    %         @                                                           #LHS    #RHS              
                                                     #MPI_GROUP b             
                                                     #MPI_GROUP b   %         @                                                           #LHS    #RHS              
                                                     #MPI_INFO f             
                                                     #MPI_INFO f   %         @                                                           #LHS    #RHS              
                                                     #MPI_MESSAGE j             
                                                     #MPI_MESSAGE j   %         @                                                           #LHS    #RHS              
                                                     #MPI_WIN n             
                                                     #MPI_WIN n   %         @                                                           #LHS    #RHS              
                                                     #MPI_DATATYPE r             
                                                     #MPI_DATATYPE r   %         @                                                           #LHS    #RHS              
                                                     #MPI_REQUEST v             
                                                     #MPI_REQUEST v   %         @                                                           #LHS    #RHS              
                                                     #MPI_FILE z             
                                                     #MPI_FILE z   %         @                                                           #LHS    #RHS              
                                                     #MPI_ERRHANDLER ~             
                                                     #MPI_ERRHANDLER ~   %         @                                                            #LHS ¡   #RHS ¢             
                                  ¡                   #MPI_COMM              
                                  ¢                   #MPI_COMM    %         @                                £                           #LHS ¤   #RHS ¥             
                                  ¤                   #MPI_OP              
                                  ¥                   #MPI_OP    #         @                                   ¦                    #LEV §   #X ¨   #Y ©             
                                 §                    
                                ¨                   
              &                   &                                                                                   ©                   
               &                   &                                           #         @                                   ª                    #LEV «   #X ¬   #Y ­             
                                 «                    
                                ¬                   
              &                   &                   &                                                                                   ­                   
               &                   &                   &                                           #         @                                   ®                    #SOLVE_P%LOG ¯   #SOLVE_P%SQRT °   #SOLVE_P%SUM ±   #SOLVE_P%MOD ²   #SOLVE_P%REAL ³                                                                    ¯     LOG                                            °     SQRT                                            ±     SUM                                            ²     MOD                             @              ³     REAL #         @                                  ´                     #         @                                  µ                    #LEV1 ¶             
                                 ¶                            @                                '                    #MPI_VAL ·                                               ·                                    @                                '                    #MPI_VAL ¸                                               ¸                                    @                           ~     '                    #MPI_VAL ¹                                               ¹                                    @                           z     '                    #MPI_VAL º                                               º                                    @                           v     '                    #MPI_VAL »                                               »                                    @                           r     '                    #MPI_VAL ¼                                               ¼                                    @                           n     '                    #MPI_VAL ½                                               ½                                    @                           j     '                    #MPI_VAL ¾                                               ¾                                    @                           f     '                    #MPI_VAL ¿                                               ¿                                    @                           b     '                    #MPI_VAL À                                               À                          "      fn#fn    Â   @   J   MG_CST      @   J   MG_MPI    B  @   J   MG_TICTOC      @   J   MG_NAMELIST    Â  @   J   MG_GRIDS      @   J   MG_INTERGRIDS    B  @   J   MG_RELAX      @   J   MG_NETCDF_OUT #   Â        GRID_TYPE+MG_GRIDS &   Ö  H   a   GRID_TYPE%NX+MG_GRIDS &     H   a   GRID_TYPE%NY+MG_GRIDS &   f  H   a   GRID_TYPE%NZ+MG_GRIDS '   ®  H   a   GRID_TYPE%NPX+MG_GRIDS '   ö  H   a   GRID_TYPE%NPY+MG_GRIDS (   >  H   a   GRID_TYPE%INCX+MG_GRIDS (     H   a   GRID_TYPE%INCY+MG_GRIDS (   Î  H   a   GRID_TYPE%NG2D+MG_GRIDS &     H   a   GRID_TYPE%NG+MG_GRIDS '   ^  H   a   GRID_TYPE%NGP+MG_GRIDS *   ¦  H   a   GRID_TYPE%GATHER+MG_GRIDS '   î  H   a   GRID_TYPE%NGX+MG_GRIDS '   6  H   a   GRID_TYPE%NGY+MG_GRIDS '   ~  H   a   GRID_TYPE%KEY+MG_GRIDS -   Æ  H   a   GRID_TYPE%LOCALCOMM+MG_GRIDS *   	     a   GRID_TYPE%NEIGHB+MG_GRIDS &   ª	  <  a   GRID_TYPE%CA+MG_GRIDS %   æ
  $  a   GRID_TYPE%P+MG_GRIDS %   
  $  a   GRID_TYPE%B+MG_GRIDS %   .  $  a   GRID_TYPE%R+MG_GRIDS &   R  $  a   GRID_TYPE%PX+MG_GRIDS &   v  $  a   GRID_TYPE%PY+MG_GRIDS &     $  a   GRID_TYPE%PZ+MG_GRIDS &   ¾    a   GRID_TYPE%DX+MG_GRIDS &   Ê    a   GRID_TYPE%DY+MG_GRIDS '   Ö    a   GRID_TYPE%DXU+MG_GRIDS '   â    a   GRID_TYPE%DYV+MG_GRIDS &   î  $  a   GRID_TYPE%DZ+MG_GRIDS '     $  a   GRID_TYPE%DZW+MG_GRIDS '   6  $  a   GRID_TYPE%ARX+MG_GRIDS '   Z  $  a   GRID_TYPE%ARY+MG_GRIDS '   ~    a   GRID_TYPE%ARZ+MG_GRIDS (       a   GRID_TYPE%BETA+MG_GRIDS (       a   GRID_TYPE%GAMU+MG_GRIDS (   ¢    a   GRID_TYPE%GAMV+MG_GRIDS (   ®  $  a   GRID_TYPE%ZXDY+MG_GRIDS (   Ò  $  a   GRID_TYPE%ZYDX+MG_GRIDS )   ö   $  a   GRID_TYPE%ALPHA+MG_GRIDS -   "  $  a   GRID_TYPE%DUMMY3DNZ+MG_GRIDS .   >#  $  a   GRID_TYPE%DUMMY3DNZP+MG_GRIDS *   b$  $  a   GRID_TYPE%DUMMY3+MG_GRIDS 2   %  <  a   GRID_TYPE%GATHERBUFFER2D+MG_GRIDS 0   Â&  T  a   GRID_TYPE%GATHERBUFFER+MG_GRIDS &   (  $  a   GRID_TYPE%DU+MG_GRIDS &   :)  $  a   GRID_TYPE%DV+MG_GRIDS &   ^*  $  a   GRID_TYPE%DW+MG_GRIDS    +  q       IP+MG_CST    ó+  q       RP+MG_CST    d,  q       LG+MG_TICTOC (   Õ,  @       OUTPUT_FREQ+MG_NAMELIST %   -  @       AUTOTUNE+MG_NAMELIST    U-  @       MYRANK+MG_MPI    -         GRID+MG_GRIDS    0.  k       TIC+MG_TICTOC #   .  =      TIC%TRIM+MG_TICTOC "   Ø.  @   a   TIC%LEV+MG_TICTOC %   /  L   a   TIC%STRING+MG_TICTOC +   d/  i       GLOBAL_SUM+MG_MPI_EXCHANGE /   Í/  @   a   GLOBAL_SUM%LEV+MG_MPI_EXCHANGE 2   0  @   a   GLOBAL_SUM%SUMLOC+MG_MPI_EXCHANGE 2   M0  @   a   GLOBAL_SUM%SUMGLO+MG_MPI_EXCHANGE *   0         COMPUTE_RESIDUAL+MG_RELAX /   1  =      COMPUTE_RESIDUAL%SIZE+MG_RELAX /   Z1  =      COMPUTE_RESIDUAL%SQRT+MG_RELAX .   1  @   a   COMPUTE_RESIDUAL%LEV+MG_RELAX .   ×1  @   a   COMPUTE_RESIDUAL%RES+MG_RELAX +   2  @       SOLVER_MAXITER+MG_NAMELIST (   W2  @       SOLVER_PREC+MG_NAMELIST    2  y       TOC+MG_TICTOC #   3  =      TOC%REAL+MG_TICTOC #   M3  =      TOC%TRIM+MG_TICTOC "   3  @   a   TOC%LEV+MG_TICTOC %   Ê3  L   a   TOC%STRING+MG_TICTOC    4  @       NLEVS+MG_GRIDS *   V4  Q       FINE2COARSE+MG_INTERGRIDS .   §4  @   a   FINE2COARSE%LEV+MG_INTERGRIDS    ç4  ~       RELAX+MG_RELAX $   e5  =      RELAX%SIZE+MG_RELAX $   ¢5  =      RELAX%TRIM+MG_RELAX #   ß5  @   a   RELAX%LEV+MG_RELAX '   6  @   a   RELAX%NSWEEPS+MG_RELAX (   _6  @       NS_COARSEST+MG_NAMELIST *   6  Q       COARSE2FINE+MG_INTERGRIDS .   ð6  @   a   COARSE2FINE%LEV+MG_INTERGRIDS #   07  @       NS_PRE+MG_NAMELIST $   p7  @       NS_POST+MG_NAMELIST &   °7  b       GROUPEQ+MPI_CONSTANTS *   8  W   a   GROUPEQ%LHS+MPI_CONSTANTS *   i8  W   a   GROUPEQ%RHS+MPI_CONSTANTS %   À8  b       INFOEQ+MPI_CONSTANTS )   "9  V   a   INFOEQ%LHS+MPI_CONSTANTS )   x9  V   a   INFOEQ%RHS+MPI_CONSTANTS (   Î9  b       MESSAGEEQ+MPI_CONSTANTS ,   0:  Y   a   MESSAGEEQ%LHS+MPI_CONSTANTS ,   :  Y   a   MESSAGEEQ%RHS+MPI_CONSTANTS $   â:  b       WINEQ+MPI_CONSTANTS (   D;  U   a   WINEQ%LHS+MPI_CONSTANTS (   ;  U   a   WINEQ%RHS+MPI_CONSTANTS )   î;  b       DATATYPEEQ+MPI_CONSTANTS -   P<  Z   a   DATATYPEEQ%LHS+MPI_CONSTANTS -   ª<  Z   a   DATATYPEEQ%RHS+MPI_CONSTANTS (   =  b       REQUESTEQ+MPI_CONSTANTS ,   f=  Y   a   REQUESTEQ%LHS+MPI_CONSTANTS ,   ¿=  Y   a   REQUESTEQ%RHS+MPI_CONSTANTS %   >  b       FILEEQ+MPI_CONSTANTS )   z>  V   a   FILEEQ%LHS+MPI_CONSTANTS )   Ð>  V   a   FILEEQ%RHS+MPI_CONSTANTS +   &?  b       ERRHANDLEREQ+MPI_CONSTANTS /   ?  \   a   ERRHANDLEREQ%LHS+MPI_CONSTANTS /   ä?  \   a   ERRHANDLEREQ%RHS+MPI_CONSTANTS %   @@  b       COMMEQ+MPI_CONSTANTS )   ¢@  V   a   COMMEQ%LHS+MPI_CONSTANTS )   ø@  V   a   COMMEQ%RHS+MPI_CONSTANTS #   NA  b       OPEQ+MPI_CONSTANTS '   °A  T   a   OPEQ%LHS+MPI_CONSTANTS '   B  T   a   OPEQ%RHS+MPI_CONSTANTS '   XB  b       GROUPNEQ+MPI_CONSTANTS +   ºB  W   a   GROUPNEQ%LHS+MPI_CONSTANTS +   C  W   a   GROUPNEQ%RHS+MPI_CONSTANTS &   hC  b       INFONEQ+MPI_CONSTANTS *   ÊC  V   a   INFONEQ%LHS+MPI_CONSTANTS *    D  V   a   INFONEQ%RHS+MPI_CONSTANTS )   vD  b       MESSAGENEQ+MPI_CONSTANTS -   ØD  Y   a   MESSAGENEQ%LHS+MPI_CONSTANTS -   1E  Y   a   MESSAGENEQ%RHS+MPI_CONSTANTS %   E  b       WINNEQ+MPI_CONSTANTS )   ìE  U   a   WINNEQ%LHS+MPI_CONSTANTS )   AF  U   a   WINNEQ%RHS+MPI_CONSTANTS *   F  b       DATATYPENEQ+MPI_CONSTANTS .   øF  Z   a   DATATYPENEQ%LHS+MPI_CONSTANTS .   RG  Z   a   DATATYPENEQ%RHS+MPI_CONSTANTS )   ¬G  b       REQUESTNEQ+MPI_CONSTANTS -   H  Y   a   REQUESTNEQ%LHS+MPI_CONSTANTS -   gH  Y   a   REQUESTNEQ%RHS+MPI_CONSTANTS &   ÀH  b       FILENEQ+MPI_CONSTANTS *   "I  V   a   FILENEQ%LHS+MPI_CONSTANTS *   xI  V   a   FILENEQ%RHS+MPI_CONSTANTS ,   ÎI  b       ERRHANDLERNEQ+MPI_CONSTANTS 0   0J  \   a   ERRHANDLERNEQ%LHS+MPI_CONSTANTS 0   J  \   a   ERRHANDLERNEQ%RHS+MPI_CONSTANTS &   èJ  b       COMMNEQ+MPI_CONSTANTS *   JK  V   a   COMMNEQ%LHS+MPI_CONSTANTS *    K  V   a   COMMNEQ%RHS+MPI_CONSTANTS $   öK  b       OPNEQ+MPI_CONSTANTS (   XL  T   a   OPNEQ%LHS+MPI_CONSTANTS (   ¬L  T   a   OPNEQ%RHS+MPI_CONSTANTS $    M  _       GATHER_2D+MG_GATHER (   _M  @   a   GATHER_2D%LEV+MG_GATHER &   M  ¤   a   GATHER_2D%X+MG_GATHER &   CN  ¤   a   GATHER_2D%Y+MG_GATHER $   çN  _       GATHER_3D+MG_GATHER (   FO  @   a   GATHER_3D%LEV+MG_GATHER &   O  ¼   a   GATHER_3D%X+MG_GATHER &   BP  ¼   a   GATHER_3D%Y+MG_GATHER    þP  µ       SOLVE_P    ³Q  <      SOLVE_P%LOG    ïQ  =      SOLVE_P%SQRT    ,R  <      SOLVE_P%SUM    hR  <      SOLVE_P%MOD    ¤R  =      SOLVE_P%REAL    áR  H       FCYCLE    )S  R       VCYCLE    {S  @   a   VCYCLE%LEV1 %   »S  ]       MPI_OP+MPI_CONSTANTS -   T  H   a   MPI_OP%MPI_VAL+MPI_CONSTANTS '   `T  ]       MPI_COMM+MPI_CONSTANTS /   ½T  H   a   MPI_COMM%MPI_VAL+MPI_CONSTANTS -   U  ]       MPI_ERRHANDLER+MPI_CONSTANTS 5   bU  H   a   MPI_ERRHANDLER%MPI_VAL+MPI_CONSTANTS '   ªU  ]       MPI_FILE+MPI_CONSTANTS /   V  H   a   MPI_FILE%MPI_VAL+MPI_CONSTANTS *   OV  ]       MPI_REQUEST+MPI_CONSTANTS 2   ¬V  H   a   MPI_REQUEST%MPI_VAL+MPI_CONSTANTS +   ôV  ]       MPI_DATATYPE+MPI_CONSTANTS 3   QW  H   a   MPI_DATATYPE%MPI_VAL+MPI_CONSTANTS &   W  ]       MPI_WIN+MPI_CONSTANTS .   öW  H   a   MPI_WIN%MPI_VAL+MPI_CONSTANTS *   >X  ]       MPI_MESSAGE+MPI_CONSTANTS 2   X  H   a   MPI_MESSAGE%MPI_VAL+MPI_CONSTANTS '   ãX  ]       MPI_INFO+MPI_CONSTANTS /   @Y  H   a   MPI_INFO%MPI_VAL+MPI_CONSTANTS (   Y  ]       MPI_GROUP+MPI_CONSTANTS 0   åY  H   a   MPI_GROUP%MPI_VAL+MPI_CONSTANTS 