	  pg  Ó   k820309    ?          14.0        +qZ                                                                                                           
       mg_projection.f90 MG_PROJECTION                                                    
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
            &                   &                   &                   &                                                                                 y
                                                                                                  ð                
            &                   &                   &                                                                                 y
                                                                                                  h               
            &                   &                   &                                                                                 y
                                                                                                  à               
            &                   &                   &                                                                                 y
                                                                                                  X               
            &                   &                   &                                                                                 y
                                                                                                  Ð               
            &                   &                   &                                                                                 y
                                                                                                   H               
            &                   &                   &                                                                                 y
                                                                                      !            À               
            &                   &                                                                                 y
                                                                                      "                            
            &                   &                                                                                 y
                                                                                      #                           
            &                   &                                                                                 y
                                                                                      $            à               
            &                   &                                                                                 y
                                                                                      %            @               
            &                   &                   &                                                                                 y
                                                                                      &            ¸               
            &                   &                   &                                                                                 y
                                                                                      '            0               
            &                   &                   &                                                                                 y
                                                                                      (            ¨               
            &                   &                   &                                                                                 y
                                                                                      )                             
            &                   &                                                                                 y
                                                                                      *                         !  
            &                   &                                                                                 y
                                                                                      +            à             "  
            &                   &                                                                                 y
                                                                                      ,            @             #  
            &                   &                                                                                 y
                                                                                      -                          $  
            &                   &                   &                                                                                 y
                                                                                      .            	             %  
            &                   &                   &                                                                                 y
                                                                                      /            	             &  
            &                   &                   &                                                                                 y
                                                                                      0            
             '  
            &                   &                   &                                                                                 y
                                                                                      1            
             (  
            &                   &                   &                                                                                 y
                                                                                      2            ø
             )  
            &                   &                   &                                                                                 y
                                                                                      3            p             *  
            &                   &                   &                   &                                                                                 y
                                                                                      4                          +  
            &                   &                   &                   &                   &                                                                                 y
                                                                                      5            ¨             ,  
            &                   &                   &                                                                                 y
                                                                                      6                          -  
            &                   &                   &                                                                                 y
                                                                                      7                         .  
            &                   &                   &                                                                                 y
                                                                                           8                                                      4                                            9                                                      8          D                                  :                                                        ;                      @ P                               <                                   &                                           #GRID_TYPE 	                                               =     
                
                       Ð?        0.25                                            >     
                
                       à?        0.5          D                                  ?                         @                              @                                                        A     
                
                       ð?        1.                                            B     
                
                        @        2.%         @                                C                           #LHS D   #RHS F             
                                  D                   #MPI_GROUP E             
                                  F                   #MPI_GROUP E   %         @                                G                           #LHS H   #RHS J             
                                  H                   #MPI_INFO I             
                                  J                   #MPI_INFO I   %         @                                K                           #LHS L   #RHS N             
                                  L                   #MPI_MESSAGE M             
                                  N                   #MPI_MESSAGE M   %         @                                O                           #LHS P   #RHS R             
                                  P                   #MPI_WIN Q             
                                  R                   #MPI_WIN Q   %         @                                S                           #LHS T   #RHS V             
                                  T                   #MPI_DATATYPE U             
                                  V                   #MPI_DATATYPE U   %         @                                W                           #LHS X   #RHS Z             
                                  X                   #MPI_REQUEST Y             
                                  Z                   #MPI_REQUEST Y   %         @                                [                           #LHS \   #RHS ^             
                                  \                   #MPI_FILE ]             
                                  ^                   #MPI_FILE ]   %         @                                _                           #LHS `   #RHS b             
                                  `                   #MPI_ERRHANDLER a             
                                  b                   #MPI_ERRHANDLER a   %         @                                c                           #LHS d   #RHS f             
                                  d                   #MPI_COMM e             
                                  f                   #MPI_COMM e   %         @                                g                           #LHS h   #RHS j             
                                  h                   #MPI_OP i             
                                  j                   #MPI_OP i   %         @                                k                           #LHS l   #RHS m             
                                  l                   #MPI_GROUP E             
                                  m                   #MPI_GROUP E   %         @                                n                           #LHS o   #RHS p             
                                  o                   #MPI_INFO I             
                                  p                   #MPI_INFO I   %         @                                q                           #LHS r   #RHS s             
                                  r                   #MPI_MESSAGE M             
                                  s                   #MPI_MESSAGE M   %         @                                t                           #LHS u   #RHS v             
                                  u                   #MPI_WIN Q             
                                  v                   #MPI_WIN Q   %         @                                w                           #LHS x   #RHS y             
                                  x                   #MPI_DATATYPE U             
                                  y                   #MPI_DATATYPE U   %         @                                z                           #LHS {   #RHS |             
                                  {                   #MPI_REQUEST Y             
                                  |                   #MPI_REQUEST Y   %         @                                }                           #LHS ~   #RHS              
                                  ~                   #MPI_FILE ]             
                                                     #MPI_FILE ]   %         @                                                           #LHS    #RHS              
                                                     #MPI_ERRHANDLER a             
                                                     #MPI_ERRHANDLER a   %         @                                                           #LHS    #RHS              
                                                     #MPI_COMM e             
                                                     #MPI_COMM e   %         @                                                           #LHS    #RHS              
                                                     #MPI_OP i             
                                                     #MPI_OP i   #         @                                                      #FILL_HALO_2D%PRESENT    #FILL_HALO_2D%TRIM    #LEV    #A2D    #LBC_NULL                                                     PRESENT                                                  TRIM           
                                                     
                                                  
               &                   &                                                     
                                                           #         @                                                      #FILL_HALO_3D%PRESENT    #FILL_HALO_3D%TRIM    #FILL_HALO_3D%SIZE    #LEV    #P    #LBC_NULL                                                     PRESENT                                                  TRIM                                                  SIZE           
                                                     
                                                 
 )              &                   &                   &                                                     
                                                           #         @                                                      #FILL_HALO_3D_RELAX%SIZE    #LEV    #P    #NX    #NY    #NZ                                                     SIZE           
                                                     
                                                   
         p            5 O p        n                                           1p         p        p         5 O p        p          5 O p          & p          5 O p        n                                      1  & p          5 O p        n                                          1    5 O p             5 O p        n                                      1p         p             5 O p        n                                          1p         p                                                                      @                                                     @                                                     @                     #         @                                                     #FILL_HALO_4D%SIZE    #LEV    #CA                                                      SIZE           
                                                     
                                                  
 =              &                   &                   &                   &                                           #         @                                   ¡                    #LEV ¢   #X £   #Y ¤             
                                 ¢                    
                                £                   
              &                   &                                                                                   ¤                   
               &                   &                                           #         @                                   ¥                    #LEV ¦   #X §   #Y ¨             
                                 ¦                    
                                §                   
              &                   &                   &                                                                                   ¨                   
               &                   &                   &                                           #         @                                   ©                   #SUB_NETCDF_WRITE_FAST_R2D%PRESENT ª   #SUB_NETCDF_WRITE_FAST_R2D%TRIM «   #SUB_NETCDF_WRITE_FAST_R2D%LEN_TRIM ¬   #SUB_NETCDF_WRITE_FAST_R2D%SIZE ­   #ARRAY_R2D ®   #VNAME ¯   #NETCDF_FILE_NAME °   #RANK ±   #ITER ²                                               ª     PRESENT                                             «     TRIM                                             ¬     LEN_TRIM                                             ­     SIZE           
                                ®                   
              &                   &                                                     
                                ¯                    1           
                                °                    1           
                                ±                     
                                ²           #         @                                   ³                   #SUB_NETCDF_WRITE_FAST_R3D%PRESENT ´   #SUB_NETCDF_WRITE_FAST_R3D%TRIM µ   #SUB_NETCDF_WRITE_FAST_R3D%LEN_TRIM ¶   #SUB_NETCDF_WRITE_FAST_R3D%SIZE ·   #ARRAY_R3D ¸   #VNAME ¹   #NETCDF_FILE_NAME º   #RANK »   #ITER ¼                                               ´     PRESENT                                             µ     TRIM                                             ¶     LEN_TRIM                                             ·     SIZE           
                                ¸                   
              &                   &                   &                                                     
                                ¹                    1           
                                º                    1           
                                »                     
                                ¼           #         @                                  ½                   #SUB_NETCDF_WRITE_FAST_R4D%PRESENT ¾   #SUB_NETCDF_WRITE_FAST_R4D%TRIM ¿   #SUB_NETCDF_WRITE_FAST_R4D%LEN_TRIM À   #SUB_NETCDF_WRITE_FAST_R4D%SIZE Á   #ARRAY_R4D Â   #VNAME Ã   #NETCDF_FILE_NAME Ä   #RANK Å   #ITER Æ                                               ¾     PRESENT                                             ¿     TRIM                                             À     LEN_TRIM                                             Á     SIZE          
                               Â                   
              &                   &                   &                   &                                                     
                                Ã                    1           
                                Ä                    1           
                                Å                     
                                Æ           #         @                                   Ç                     #         @                                   È                                      @                           i     '                    #MPI_VAL É                                               É                                    @                           e     '                    #MPI_VAL Ê                                               Ê                                    @                           a     '                    #MPI_VAL Ë                                               Ë                                    @                           ]     '                    #MPI_VAL Ì                                               Ì                                    @                           Y     '                    #MPI_VAL Í                                               Í                                    @                           U     '                    #MPI_VAL Î                                               Î                                    @                           Q     '                    #MPI_VAL Ï                                               Ï                                    @                           M     '                    #MPI_VAL Ð                                               Ð                                    @                           I     '                    #MPI_VAL Ñ                                               Ñ                                    @                           E     '                    #MPI_VAL Ò                                               Ò                          (      fn#fn    È   @   J   MG_CST      @   J   MG_MPI    H  @   J   MG_TICTOC      @   J   MG_NAMELIST    È  @   J   MG_GRIDS       @   J   MG_MPI_EXCHANGE    H  @   J   MG_GATHER      @   J   MG_NETCDF_OUT #   È        GRID_TYPE+MG_GRIDS &   Ü  H   a   GRID_TYPE%NX+MG_GRIDS &   $  H   a   GRID_TYPE%NY+MG_GRIDS &   l  H   a   GRID_TYPE%NZ+MG_GRIDS '   ´  H   a   GRID_TYPE%NPX+MG_GRIDS '   ü  H   a   GRID_TYPE%NPY+MG_GRIDS (   D  H   a   GRID_TYPE%INCX+MG_GRIDS (     H   a   GRID_TYPE%INCY+MG_GRIDS (   Ô  H   a   GRID_TYPE%NG2D+MG_GRIDS &     H   a   GRID_TYPE%NG+MG_GRIDS '   d  H   a   GRID_TYPE%NGP+MG_GRIDS *   ¬  H   a   GRID_TYPE%GATHER+MG_GRIDS '   ô  H   a   GRID_TYPE%NGX+MG_GRIDS '   <  H   a   GRID_TYPE%NGY+MG_GRIDS '     H   a   GRID_TYPE%KEY+MG_GRIDS -   Ì  H   a   GRID_TYPE%LOCALCOMM+MG_GRIDS *   	     a   GRID_TYPE%NEIGHB+MG_GRIDS &   °	  <  a   GRID_TYPE%CA+MG_GRIDS %   ì
  $  a   GRID_TYPE%P+MG_GRIDS %     $  a   GRID_TYPE%B+MG_GRIDS %   4  $  a   GRID_TYPE%R+MG_GRIDS &   X  $  a   GRID_TYPE%PX+MG_GRIDS &   |  $  a   GRID_TYPE%PY+MG_GRIDS &      $  a   GRID_TYPE%PZ+MG_GRIDS &   Ä    a   GRID_TYPE%DX+MG_GRIDS &   Ð    a   GRID_TYPE%DY+MG_GRIDS '   Ü    a   GRID_TYPE%DXU+MG_GRIDS '   è    a   GRID_TYPE%DYV+MG_GRIDS &   ô  $  a   GRID_TYPE%DZ+MG_GRIDS '     $  a   GRID_TYPE%DZW+MG_GRIDS '   <  $  a   GRID_TYPE%ARX+MG_GRIDS '   `  $  a   GRID_TYPE%ARY+MG_GRIDS '       a   GRID_TYPE%ARZ+MG_GRIDS (       a   GRID_TYPE%BETA+MG_GRIDS (       a   GRID_TYPE%GAMU+MG_GRIDS (   ¨    a   GRID_TYPE%GAMV+MG_GRIDS (   ´  $  a   GRID_TYPE%ZXDY+MG_GRIDS (   Ø  $  a   GRID_TYPE%ZYDX+MG_GRIDS )   ü   $  a   GRID_TYPE%ALPHA+MG_GRIDS -    "  $  a   GRID_TYPE%DUMMY3DNZ+MG_GRIDS .   D#  $  a   GRID_TYPE%DUMMY3DNZP+MG_GRIDS *   h$  $  a   GRID_TYPE%DUMMY3+MG_GRIDS 2   %  <  a   GRID_TYPE%GATHERBUFFER2D+MG_GRIDS 0   È&  T  a   GRID_TYPE%GATHERBUFFER+MG_GRIDS &   (  $  a   GRID_TYPE%DU+MG_GRIDS &   @)  $  a   GRID_TYPE%DV+MG_GRIDS &   d*  $  a   GRID_TYPE%DW+MG_GRIDS    +  q       IP+MG_CST    ù+  q       RP+MG_CST ,   j,  @       SURFACE_NEUMANN+MG_NAMELIST    ª,  @       NLEVS+MG_GRIDS    ê,         GRID+MG_GRIDS    -  t       QRT+MG_CST    ù-  s       HLF+MG_CST *   l.  @       NETCDF_OUTPUT+MG_NAMELIST    ¬.  @       MYRANK+MG_MPI    ì.  r       ONE+MG_CST    ^/  r       TWO+MG_CST &   Ð/  b       GROUPEQ+MPI_CONSTANTS *   20  W   a   GROUPEQ%LHS+MPI_CONSTANTS *   0  W   a   GROUPEQ%RHS+MPI_CONSTANTS %   à0  b       INFOEQ+MPI_CONSTANTS )   B1  V   a   INFOEQ%LHS+MPI_CONSTANTS )   1  V   a   INFOEQ%RHS+MPI_CONSTANTS (   î1  b       MESSAGEEQ+MPI_CONSTANTS ,   P2  Y   a   MESSAGEEQ%LHS+MPI_CONSTANTS ,   ©2  Y   a   MESSAGEEQ%RHS+MPI_CONSTANTS $   3  b       WINEQ+MPI_CONSTANTS (   d3  U   a   WINEQ%LHS+MPI_CONSTANTS (   ¹3  U   a   WINEQ%RHS+MPI_CONSTANTS )   4  b       DATATYPEEQ+MPI_CONSTANTS -   p4  Z   a   DATATYPEEQ%LHS+MPI_CONSTANTS -   Ê4  Z   a   DATATYPEEQ%RHS+MPI_CONSTANTS (   $5  b       REQUESTEQ+MPI_CONSTANTS ,   5  Y   a   REQUESTEQ%LHS+MPI_CONSTANTS ,   ß5  Y   a   REQUESTEQ%RHS+MPI_CONSTANTS %   86  b       FILEEQ+MPI_CONSTANTS )   6  V   a   FILEEQ%LHS+MPI_CONSTANTS )   ð6  V   a   FILEEQ%RHS+MPI_CONSTANTS +   F7  b       ERRHANDLEREQ+MPI_CONSTANTS /   ¨7  \   a   ERRHANDLEREQ%LHS+MPI_CONSTANTS /   8  \   a   ERRHANDLEREQ%RHS+MPI_CONSTANTS %   `8  b       COMMEQ+MPI_CONSTANTS )   Â8  V   a   COMMEQ%LHS+MPI_CONSTANTS )   9  V   a   COMMEQ%RHS+MPI_CONSTANTS #   n9  b       OPEQ+MPI_CONSTANTS '   Ð9  T   a   OPEQ%LHS+MPI_CONSTANTS '   $:  T   a   OPEQ%RHS+MPI_CONSTANTS '   x:  b       GROUPNEQ+MPI_CONSTANTS +   Ú:  W   a   GROUPNEQ%LHS+MPI_CONSTANTS +   1;  W   a   GROUPNEQ%RHS+MPI_CONSTANTS &   ;  b       INFONEQ+MPI_CONSTANTS *   ê;  V   a   INFONEQ%LHS+MPI_CONSTANTS *   @<  V   a   INFONEQ%RHS+MPI_CONSTANTS )   <  b       MESSAGENEQ+MPI_CONSTANTS -   ø<  Y   a   MESSAGENEQ%LHS+MPI_CONSTANTS -   Q=  Y   a   MESSAGENEQ%RHS+MPI_CONSTANTS %   ª=  b       WINNEQ+MPI_CONSTANTS )   >  U   a   WINNEQ%LHS+MPI_CONSTANTS )   a>  U   a   WINNEQ%RHS+MPI_CONSTANTS *   ¶>  b       DATATYPENEQ+MPI_CONSTANTS .   ?  Z   a   DATATYPENEQ%LHS+MPI_CONSTANTS .   r?  Z   a   DATATYPENEQ%RHS+MPI_CONSTANTS )   Ì?  b       REQUESTNEQ+MPI_CONSTANTS -   .@  Y   a   REQUESTNEQ%LHS+MPI_CONSTANTS -   @  Y   a   REQUESTNEQ%RHS+MPI_CONSTANTS &   à@  b       FILENEQ+MPI_CONSTANTS *   BA  V   a   FILENEQ%LHS+MPI_CONSTANTS *   A  V   a   FILENEQ%RHS+MPI_CONSTANTS ,   îA  b       ERRHANDLERNEQ+MPI_CONSTANTS 0   PB  \   a   ERRHANDLERNEQ%LHS+MPI_CONSTANTS 0   ¬B  \   a   ERRHANDLERNEQ%RHS+MPI_CONSTANTS &   C  b       COMMNEQ+MPI_CONSTANTS *   jC  V   a   COMMNEQ%LHS+MPI_CONSTANTS *   ÀC  V   a   COMMNEQ%RHS+MPI_CONSTANTS $   D  b       OPNEQ+MPI_CONSTANTS (   xD  T   a   OPNEQ%LHS+MPI_CONSTANTS (   ÌD  T   a   OPNEQ%RHS+MPI_CONSTANTS -    E         FILL_HALO_2D+MG_MPI_EXCHANGE 5   ¹E  @      FILL_HALO_2D%PRESENT+MG_MPI_EXCHANGE 2   ùE  =      FILL_HALO_2D%TRIM+MG_MPI_EXCHANGE 1   6F  @   a   FILL_HALO_2D%LEV+MG_MPI_EXCHANGE 1   vF  ¤   a   FILL_HALO_2D%A2D+MG_MPI_EXCHANGE 6   G  P   a   FILL_HALO_2D%LBC_NULL+MG_MPI_EXCHANGE -   jG  ®       FILL_HALO_3D+MG_MPI_EXCHANGE 5   H  @      FILL_HALO_3D%PRESENT+MG_MPI_EXCHANGE 2   XH  =      FILL_HALO_3D%TRIM+MG_MPI_EXCHANGE 2   H  =      FILL_HALO_3D%SIZE+MG_MPI_EXCHANGE 1   ÒH  @   a   FILL_HALO_3D%LEV+MG_MPI_EXCHANGE /   I  ¼   a   FILL_HALO_3D%P+MG_MPI_EXCHANGE 6   ÎI  P   a   FILL_HALO_3D%LBC_NULL+MG_MPI_EXCHANGE 3   J         FILL_HALO_3D_RELAX+MG_MPI_EXCHANGE 8   «J  =      FILL_HALO_3D_RELAX%SIZE+MG_MPI_EXCHANGE 7   èJ  @   a   FILL_HALO_3D_RELAX%LEV+MG_MPI_EXCHANGE 5   (K  ñ  a   FILL_HALO_3D_RELAX%P+MG_MPI_EXCHANGE 6   N  @   a   FILL_HALO_3D_RELAX%NX+MG_MPI_EXCHANGE 6   YN  @   a   FILL_HALO_3D_RELAX%NY+MG_MPI_EXCHANGE 6   N  @   a   FILL_HALO_3D_RELAX%NZ+MG_MPI_EXCHANGE -   ÙN  p       FILL_HALO_4D+MG_MPI_EXCHANGE 2   IO  =      FILL_HALO_4D%SIZE+MG_MPI_EXCHANGE 1   O  @   a   FILL_HALO_4D%LEV+MG_MPI_EXCHANGE 0   ÆO  Ô   a   FILL_HALO_4D%CA+MG_MPI_EXCHANGE $   P  _       GATHER_2D+MG_GATHER (   ùP  @   a   GATHER_2D%LEV+MG_GATHER &   9Q  ¤   a   GATHER_2D%X+MG_GATHER &   ÝQ  ¤   a   GATHER_2D%Y+MG_GATHER $   R  _       GATHER_3D+MG_GATHER (   àR  @   a   GATHER_3D%LEV+MG_GATHER &    S  ¼   a   GATHER_3D%X+MG_GATHER &   ÜS  ¼   a   GATHER_3D%Y+MG_GATHER 8   T  #      SUB_NETCDF_WRITE_FAST_R2D+MG_NETCDF_OUT @   »U  @      SUB_NETCDF_WRITE_FAST_R2D%PRESENT+MG_NETCDF_OUT =   ûU  =      SUB_NETCDF_WRITE_FAST_R2D%TRIM+MG_NETCDF_OUT A   8V  A      SUB_NETCDF_WRITE_FAST_R2D%LEN_TRIM+MG_NETCDF_OUT =   yV  =      SUB_NETCDF_WRITE_FAST_R2D%SIZE+MG_NETCDF_OUT B   ¶V  ¤   a   SUB_NETCDF_WRITE_FAST_R2D%ARRAY_R2D+MG_NETCDF_OUT >   ZW  L   a   SUB_NETCDF_WRITE_FAST_R2D%VNAME+MG_NETCDF_OUT I   ¦W  L   a   SUB_NETCDF_WRITE_FAST_R2D%NETCDF_FILE_NAME+MG_NETCDF_OUT =   òW  @   a   SUB_NETCDF_WRITE_FAST_R2D%RANK+MG_NETCDF_OUT =   2X  @   a   SUB_NETCDF_WRITE_FAST_R2D%ITER+MG_NETCDF_OUT 8   rX  #      SUB_NETCDF_WRITE_FAST_R3D+MG_NETCDF_OUT @   Y  @      SUB_NETCDF_WRITE_FAST_R3D%PRESENT+MG_NETCDF_OUT =   ÕY  =      SUB_NETCDF_WRITE_FAST_R3D%TRIM+MG_NETCDF_OUT A   Z  A      SUB_NETCDF_WRITE_FAST_R3D%LEN_TRIM+MG_NETCDF_OUT =   SZ  =      SUB_NETCDF_WRITE_FAST_R3D%SIZE+MG_NETCDF_OUT B   Z  ¼   a   SUB_NETCDF_WRITE_FAST_R3D%ARRAY_R3D+MG_NETCDF_OUT >   L[  L   a   SUB_NETCDF_WRITE_FAST_R3D%VNAME+MG_NETCDF_OUT I   [  L   a   SUB_NETCDF_WRITE_FAST_R3D%NETCDF_FILE_NAME+MG_NETCDF_OUT =   ä[  @   a   SUB_NETCDF_WRITE_FAST_R3D%RANK+MG_NETCDF_OUT =   $\  @   a   SUB_NETCDF_WRITE_FAST_R3D%ITER+MG_NETCDF_OUT 8   d\  #      SUB_NETCDF_WRITE_FAST_R4D+MG_NETCDF_OUT @   ]  @      SUB_NETCDF_WRITE_FAST_R4D%PRESENT+MG_NETCDF_OUT =   Ç]  =      SUB_NETCDF_WRITE_FAST_R4D%TRIM+MG_NETCDF_OUT A   ^  A      SUB_NETCDF_WRITE_FAST_R4D%LEN_TRIM+MG_NETCDF_OUT =   E^  =      SUB_NETCDF_WRITE_FAST_R4D%SIZE+MG_NETCDF_OUT B   ^  Ô   a   SUB_NETCDF_WRITE_FAST_R4D%ARRAY_R4D+MG_NETCDF_OUT >   V_  L   a   SUB_NETCDF_WRITE_FAST_R4D%VNAME+MG_NETCDF_OUT I   ¢_  L   a   SUB_NETCDF_WRITE_FAST_R4D%NETCDF_FILE_NAME+MG_NETCDF_OUT =   î_  @   a   SUB_NETCDF_WRITE_FAST_R4D%RANK+MG_NETCDF_OUT =   .`  @   a   SUB_NETCDF_WRITE_FAST_R4D%ITER+MG_NETCDF_OUT    n`  H       SET_MATRICES    ¶`  H       CORRECTION_UVW %   þ`  ]       MPI_OP+MPI_CONSTANTS -   [a  H   a   MPI_OP%MPI_VAL+MPI_CONSTANTS '   £a  ]       MPI_COMM+MPI_CONSTANTS /    b  H   a   MPI_COMM%MPI_VAL+MPI_CONSTANTS -   Hb  ]       MPI_ERRHANDLER+MPI_CONSTANTS 5   ¥b  H   a   MPI_ERRHANDLER%MPI_VAL+MPI_CONSTANTS '   íb  ]       MPI_FILE+MPI_CONSTANTS /   Jc  H   a   MPI_FILE%MPI_VAL+MPI_CONSTANTS *   c  ]       MPI_REQUEST+MPI_CONSTANTS 2   ïc  H   a   MPI_REQUEST%MPI_VAL+MPI_CONSTANTS +   7d  ]       MPI_DATATYPE+MPI_CONSTANTS 3   d  H   a   MPI_DATATYPE%MPI_VAL+MPI_CONSTANTS &   Üd  ]       MPI_WIN+MPI_CONSTANTS .   9e  H   a   MPI_WIN%MPI_VAL+MPI_CONSTANTS *   e  ]       MPI_MESSAGE+MPI_CONSTANTS 2   Þe  H   a   MPI_MESSAGE%MPI_VAL+MPI_CONSTANTS '   &f  ]       MPI_INFO+MPI_CONSTANTS /   f  H   a   MPI_INFO%MPI_VAL+MPI_CONSTANTS (   Ëf  ]       MPI_GROUP+MPI_CONSTANTS 0   (g  H   a   MPI_GROUP%MPI_VAL+MPI_CONSTANTS 