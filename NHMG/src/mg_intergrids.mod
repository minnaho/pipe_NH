	  j  �   k820309    h          14.0        ��]                                                                                                           
       mg_intergrids.f90 MG_INTERGRIDS                                                    
                                                          
                                                          
                                                          
                                                          
                                                          
                                                          
                         @                               '             .      #NX 	   #NY 
   #NZ    #NPX    #NPY    #INCX    #INCY    #NG2D    #NG    #NGP    #GATHER    #NGX    #NGY    #KEY    #LOCALCOMM    #NEIGHB    #CA    #P    #B    #R    #PX    #PY    #PZ    #DX     #DY !   #DXU "   #DYV #   #DZ $   #DZW %   #ARX &   #ARY '   #ARZ (   #BETA )   #GAMU *   #GAMV +   #ZXDY ,   #ZYDX -   #ALPHA .   #DUMMY3DNZ /   #DUMMY3DNZP 0   #DUMMY3 1   #GATHERBUFFER2D 2   #GATHERBUFFER 3   #DU 4   #DV 5   #DW 6                �                              	                                �                              
                               �                                                             �                                                             �                                                             �                                                             �                                                             �                                                             �                                           	                   �                                   $       
                   �                                   (                          �                                   ,                          �                                   0                          �                                   4                          �                                   8                          �                                          <                   p          p            p                                    �                                        `                
            &                   &                   &                   &                                                                   �              y
                                                          �                                        �                
            &                   &                   &                                                                   �              y
                                                          �                                        h               
            &                   &                   &                                                                   �              y
                                                          �                                        �               
            &                   &                   &                                                                   �              y
                                                          �                                        X               
            &                   &                   &                                                                   �              y
                                                          �                                        �               
            &                   &                   &                                                                   �              y
                                                          �                                        H               
            &                   &                   &                                                                   �              y
                                                          �                                         �               
            &                   &                                                                   �              y
                                                          �                            !                            
            &                   &                                                                   �              y
                                                          �                            "            �               
            &                   &                                                                   �              y
                                                          �                            #            �               
            &                   &                                                                   �              y
                                                          �                            $            @               
            &                   &                   &                                                                   �              y
                                                          �                            %            �               
            &                   &                   &                                                                   �              y
                                                          �                            &            0               
            &                   &                   &                                                                   �              y
                                                          �                            '            �               
            &                   &                   &                                                                   �              y
                                                          �                            (                             
            &                   &                                                                   �              y
                                                          �                            )            �             !  
            &                   &                                                                   �              y
                                                          �                            *            �             "  
            &                   &                                                                   �              y
                                                          �                            +            @             #  
            &                   &                                                                   �              y
                                                          �                            ,            �             $  
            &                   &                   &                                                                   �              y
                                                          �                            -            	             %  
            &                   &                   &                                                                   �              y
                                                          �                            .            �	             &  
            &                   &                   &                                                                   �              y
                                                          �                            /            
             '  
            &                   &                   &                                                                   �              y
                                                          �                            0            �
             (  
            &                   &                   &                                                                   �              y
                                                          �                            1            �
             )  
            &                   &                   &                                                                   �              y
                                                          �                            2            p             *  
            &                   &                   &                   &                                                                   �              y
                                                          �                            3                          +  
            &                   &                   &                   &                   &                                                                   �              y
                                                          �                            4            �             ,  
            &                   &                   &                                                                   �              y
                                                          �                            5                          -  
            &                   &                   &                                                                   �              y
                                                          �                            6            �             .  
            &                   &                   &                                                                   �              y
                                                                                           7                                                      4                                            8                                                      8          @ P                               9                                   &                                           #GRID_TYPE    #         @                                  :                   #TIC%TRIM ;   #LEV <   #STRING =                                               ;     TRIM           
                                 <                     
                                 =                    1 #         @                                  >                   #TOC%REAL ?   #TOC%TRIM @   #LEV A   #STRING B                                @              ?     REAL                                             @     TRIM           
                                 A                     
                                 B                    1                                             C     
                
                                 0.#         @                                  D                   #SPLIT%MOD E   #LEV F   #X G   #Y H                                               E     MOD           
                                 F                    
                                G                   
              &                   &                   &                                                                                   H                   
               &                   &                   &                                                     D                                  I                                                        J     
                
                       �?        0.5%         @                                K                           #LHS L   #RHS N             
                                  L                   #MPI_GROUP M             
                                  N                   #MPI_GROUP M   %         @                                O                           #LHS P   #RHS R             
                                  P                   #MPI_INFO Q             
                                  R                   #MPI_INFO Q   %         @                                S                           #LHS T   #RHS V             
                                  T                   #MPI_MESSAGE U             
                                  V                   #MPI_MESSAGE U   %         @                                W                           #LHS X   #RHS Z             
                                  X                   #MPI_WIN Y             
                                  Z                   #MPI_WIN Y   %         @                                [                           #LHS \   #RHS ^             
                                  \                   #MPI_DATATYPE ]             
                                  ^                   #MPI_DATATYPE ]   %         @                                _                           #LHS `   #RHS b             
                                  `                   #MPI_REQUEST a             
                                  b                   #MPI_REQUEST a   %         @                                c                           #LHS d   #RHS f             
                                  d                   #MPI_FILE e             
                                  f                   #MPI_FILE e   %         @                                g                           #LHS h   #RHS j             
                                  h                   #MPI_ERRHANDLER i             
                                  j                   #MPI_ERRHANDLER i   %         @                                k                           #LHS l   #RHS n             
                                  l                   #MPI_COMM m             
                                  n                   #MPI_COMM m   %         @                                o                           #LHS p   #RHS r             
                                  p                   #MPI_OP q             
                                  r                   #MPI_OP q   %         @                                s                           #LHS t   #RHS u             
                                  t                   #MPI_GROUP M             
                                  u                   #MPI_GROUP M   %         @                                v                           #LHS w   #RHS x             
                                  w                   #MPI_INFO Q             
                                  x                   #MPI_INFO Q   %         @                                y                           #LHS z   #RHS {             
                                  z                   #MPI_MESSAGE U             
                                  {                   #MPI_MESSAGE U   %         @                                |                           #LHS }   #RHS ~             
                                  }                   #MPI_WIN Y             
                                  ~                   #MPI_WIN Y   %         @                                                           #LHS �   #RHS �             
                                  �                   #MPI_DATATYPE ]             
                                  �                   #MPI_DATATYPE ]   %         @                                �                           #LHS �   #RHS �             
                                  �                   #MPI_REQUEST a             
                                  �                   #MPI_REQUEST a   %         @                                �                           #LHS �   #RHS �             
                                  �                   #MPI_FILE e             
                                  �                   #MPI_FILE e   %         @                                �                           #LHS �   #RHS �             
                                  �                   #MPI_ERRHANDLER i             
                                  �                   #MPI_ERRHANDLER i   %         @                                �                           #LHS �   #RHS �             
                                  �                   #MPI_COMM m             
                                  �                   #MPI_COMM m   %         @                                �                           #LHS �   #RHS �             
                                  �                   #MPI_OP q             
                                  �                   #MPI_OP q   #         @                                   �                   #FILL_HALO_2D%PRESENT �   #FILL_HALO_2D%TRIM �   #LEV �   #A2D �   #LBC_NULL �                                               �     PRESENT                                             �     TRIM           
                                 �                    
                               �                   
               &                   &                                                     
                                �                           #         @                                  �                   #FILL_HALO_3D%PRESENT �   #FILL_HALO_3D%TRIM �   #FILL_HALO_3D%SIZE �   #LEV �   #P �   #LBC_NULL �                                               �     PRESENT                                             �     TRIM                                             �     SIZE           
                                 �                    
                              �                   
 )              &                   &                   &                                                     
                                �                           #         @                                   �                   #FILL_HALO_3D_RELAX%SIZE �   #LEV �   #P �   #NX �   #NY �   #NZ �                                               �     SIZE           
                                 �                    
    �                           �                    
         p            5 O p        n                                           1p         p        p         5 O p        p          5 O p          & p          5 O p        n                                      1  & p          5 O p        n                                          1    5 O p             5 O p        n                                      1p         p             5 O p        n                                          1p         p                                                                      @         �                                            @         �                                            @         �            #         @                                   �                   #FILL_HALO_4D%SIZE �   #LEV �   #CA �                                               �     SIZE           
                                 �                    
                              �                   
 =              &                   &                   &                   &                                           #         @                                   �                    #LEV �   #X �   #Y �             
                                 �                    
                                �                   
              &                   &                                                                                   �                   
               &                   &                                           #         @                                  �                    #LEV �   #X �   #Y �             
                                 �                    
                                �                   
              &                   &                   &                                                                                   �                   
               &                   &                   &                                           #         @                                   �                    #LEV �             
  @                              �           #         @                                  �                   #FINE2COARSE_2D%SIZE �   #X �   #Y �   #NX �   #NY �                                              �     SIZE          
 P                              �                   
              &                   &                   &                                                    D                               �                   
               &                   &                   &                                                     
                                 �                     
                                 �           #         @                                  �                    #X �   #Y �   #NX �   #NY �   #NZ �            
                                �                   
              &                   &                   &                                                    D                               �                   
               &                   &                   &                                                     
                                 �                     
                                 �                     
                                 �           #         @                                   �                    #LEV �             
  @                              �           #         @                                  �                    #XF �   #XC �   #NX �   #NY �            D                               �                   
 	              &                   &                   &                                                    
                                �                   
 
             &                   &                   &                                                     
                                 �                     
                                 �           #         @                                  �                   #COARSE2FINE_3D_LINEAR%MOD �   #XF �   #XC �   #NX �   #NY �   #NZ �                                              �     MOD          D                               �                   
               &                   &                   &                                                    
                                �                   
              &                   &                   &                                                     
                                 �                     
                                 �                     
                                 �                            @                           q     '                    #MPI_VAL �                �                               �                                    @                           m     '                    #MPI_VAL �                �                               �                                    @                           i     '                    #MPI_VAL �                �                               �                                    @                           e     '                    #MPI_VAL �                �                               �                                    @                           a     '                    #MPI_VAL �                �                               �                                    @                           ]     '                    #MPI_VAL �                �                               �                                    @                           Y     '                    #MPI_VAL �                �                               �                                    @                           U     '                    #MPI_VAL �                �                               �                                    @                           Q     '                    #MPI_VAL �                �                               �                                    @                           M     '                    #MPI_VAL �                �                               �                      �   (      fn#fn    �   @   J   MG_CST      @   J   MG_MPI    H  @   J   MG_TICTOC    �  @   J   MG_NAMELIST    �  @   J   MG_GRIDS       @   J   MG_MPI_EXCHANGE    H  @   J   MG_GATHER #   �        GRID_TYPE+MG_GRIDS &   �  H   a   GRID_TYPE%NX+MG_GRIDS &   �  H   a   GRID_TYPE%NY+MG_GRIDS &   ,  H   a   GRID_TYPE%NZ+MG_GRIDS '   t  H   a   GRID_TYPE%NPX+MG_GRIDS '   �  H   a   GRID_TYPE%NPY+MG_GRIDS (     H   a   GRID_TYPE%INCX+MG_GRIDS (   L  H   a   GRID_TYPE%INCY+MG_GRIDS (   �  H   a   GRID_TYPE%NG2D+MG_GRIDS &   �  H   a   GRID_TYPE%NG+MG_GRIDS '   $  H   a   GRID_TYPE%NGP+MG_GRIDS *   l  H   a   GRID_TYPE%GATHER+MG_GRIDS '   �  H   a   GRID_TYPE%NGX+MG_GRIDS '   �  H   a   GRID_TYPE%NGY+MG_GRIDS '   D  H   a   GRID_TYPE%KEY+MG_GRIDS -   �  H   a   GRID_TYPE%LOCALCOMM+MG_GRIDS *   �  �   a   GRID_TYPE%NEIGHB+MG_GRIDS &   p	  <  a   GRID_TYPE%CA+MG_GRIDS %   �
  $  a   GRID_TYPE%P+MG_GRIDS %   �  $  a   GRID_TYPE%B+MG_GRIDS %   �  $  a   GRID_TYPE%R+MG_GRIDS &     $  a   GRID_TYPE%PX+MG_GRIDS &   <  $  a   GRID_TYPE%PY+MG_GRIDS &   `  $  a   GRID_TYPE%PZ+MG_GRIDS &   �    a   GRID_TYPE%DX+MG_GRIDS &   �    a   GRID_TYPE%DY+MG_GRIDS '   �    a   GRID_TYPE%DXU+MG_GRIDS '   �    a   GRID_TYPE%DYV+MG_GRIDS &   �  $  a   GRID_TYPE%DZ+MG_GRIDS '   �  $  a   GRID_TYPE%DZW+MG_GRIDS '   �  $  a   GRID_TYPE%ARX+MG_GRIDS '      $  a   GRID_TYPE%ARY+MG_GRIDS '   D    a   GRID_TYPE%ARZ+MG_GRIDS (   P    a   GRID_TYPE%BETA+MG_GRIDS (   \    a   GRID_TYPE%GAMU+MG_GRIDS (   h    a   GRID_TYPE%GAMV+MG_GRIDS (   t  $  a   GRID_TYPE%ZXDY+MG_GRIDS (   �  $  a   GRID_TYPE%ZYDX+MG_GRIDS )   �   $  a   GRID_TYPE%ALPHA+MG_GRIDS -   �!  $  a   GRID_TYPE%DUMMY3DNZ+MG_GRIDS .   #  $  a   GRID_TYPE%DUMMY3DNZP+MG_GRIDS *   ($  $  a   GRID_TYPE%DUMMY3+MG_GRIDS 2   L%  <  a   GRID_TYPE%GATHERBUFFER2D+MG_GRIDS 0   �&  T  a   GRID_TYPE%GATHERBUFFER+MG_GRIDS &   �'  $  a   GRID_TYPE%DU+MG_GRIDS &    )  $  a   GRID_TYPE%DV+MG_GRIDS &   $*  $  a   GRID_TYPE%DW+MG_GRIDS    H+  q       IP+MG_CST    �+  q       RP+MG_CST    *,  �       GRID+MG_GRIDS    �,  k       TIC+MG_TICTOC #   0-  =      TIC%TRIM+MG_TICTOC "   m-  @   a   TIC%LEV+MG_TICTOC %   �-  L   a   TIC%STRING+MG_TICTOC    �-  y       TOC+MG_TICTOC #   r.  =      TOC%REAL+MG_TICTOC #   �.  =      TOC%TRIM+MG_TICTOC "   �.  @   a   TOC%LEV+MG_TICTOC %   ,/  L   a   TOC%STRING+MG_TICTOC    x/  r       ZERO+MG_CST     �/  n       SPLIT+MG_GATHER $   X0  <      SPLIT%MOD+MG_GATHER $   �0  @   a   SPLIT%LEV+MG_GATHER "   �0  �   a   SPLIT%X+MG_GATHER "   �1  �   a   SPLIT%Y+MG_GATHER ,   L2  @       SURFACE_NEUMANN+MG_NAMELIST    �2  s       HLF+MG_CST &   �2  b       GROUPEQ+MPI_CONSTANTS *   a3  W   a   GROUPEQ%LHS+MPI_CONSTANTS *   �3  W   a   GROUPEQ%RHS+MPI_CONSTANTS %   4  b       INFOEQ+MPI_CONSTANTS )   q4  V   a   INFOEQ%LHS+MPI_CONSTANTS )   �4  V   a   INFOEQ%RHS+MPI_CONSTANTS (   5  b       MESSAGEEQ+MPI_CONSTANTS ,   5  Y   a   MESSAGEEQ%LHS+MPI_CONSTANTS ,   �5  Y   a   MESSAGEEQ%RHS+MPI_CONSTANTS $   16  b       WINEQ+MPI_CONSTANTS (   �6  U   a   WINEQ%LHS+MPI_CONSTANTS (   �6  U   a   WINEQ%RHS+MPI_CONSTANTS )   =7  b       DATATYPEEQ+MPI_CONSTANTS -   �7  Z   a   DATATYPEEQ%LHS+MPI_CONSTANTS -   �7  Z   a   DATATYPEEQ%RHS+MPI_CONSTANTS (   S8  b       REQUESTEQ+MPI_CONSTANTS ,   �8  Y   a   REQUESTEQ%LHS+MPI_CONSTANTS ,   9  Y   a   REQUESTEQ%RHS+MPI_CONSTANTS %   g9  b       FILEEQ+MPI_CONSTANTS )   �9  V   a   FILEEQ%LHS+MPI_CONSTANTS )   :  V   a   FILEEQ%RHS+MPI_CONSTANTS +   u:  b       ERRHANDLEREQ+MPI_CONSTANTS /   �:  \   a   ERRHANDLEREQ%LHS+MPI_CONSTANTS /   3;  \   a   ERRHANDLEREQ%RHS+MPI_CONSTANTS %   �;  b       COMMEQ+MPI_CONSTANTS )   �;  V   a   COMMEQ%LHS+MPI_CONSTANTS )   G<  V   a   COMMEQ%RHS+MPI_CONSTANTS #   �<  b       OPEQ+MPI_CONSTANTS '   �<  T   a   OPEQ%LHS+MPI_CONSTANTS '   S=  T   a   OPEQ%RHS+MPI_CONSTANTS '   �=  b       GROUPNEQ+MPI_CONSTANTS +   	>  W   a   GROUPNEQ%LHS+MPI_CONSTANTS +   `>  W   a   GROUPNEQ%RHS+MPI_CONSTANTS &   �>  b       INFONEQ+MPI_CONSTANTS *   ?  V   a   INFONEQ%LHS+MPI_CONSTANTS *   o?  V   a   INFONEQ%RHS+MPI_CONSTANTS )   �?  b       MESSAGENEQ+MPI_CONSTANTS -   '@  Y   a   MESSAGENEQ%LHS+MPI_CONSTANTS -   �@  Y   a   MESSAGENEQ%RHS+MPI_CONSTANTS %   �@  b       WINNEQ+MPI_CONSTANTS )   ;A  U   a   WINNEQ%LHS+MPI_CONSTANTS )   �A  U   a   WINNEQ%RHS+MPI_CONSTANTS *   �A  b       DATATYPENEQ+MPI_CONSTANTS .   GB  Z   a   DATATYPENEQ%LHS+MPI_CONSTANTS .   �B  Z   a   DATATYPENEQ%RHS+MPI_CONSTANTS )   �B  b       REQUESTNEQ+MPI_CONSTANTS -   ]C  Y   a   REQUESTNEQ%LHS+MPI_CONSTANTS -   �C  Y   a   REQUESTNEQ%RHS+MPI_CONSTANTS &   D  b       FILENEQ+MPI_CONSTANTS *   qD  V   a   FILENEQ%LHS+MPI_CONSTANTS *   �D  V   a   FILENEQ%RHS+MPI_CONSTANTS ,   E  b       ERRHANDLERNEQ+MPI_CONSTANTS 0   E  \   a   ERRHANDLERNEQ%LHS+MPI_CONSTANTS 0   �E  \   a   ERRHANDLERNEQ%RHS+MPI_CONSTANTS &   7F  b       COMMNEQ+MPI_CONSTANTS *   �F  V   a   COMMNEQ%LHS+MPI_CONSTANTS *   �F  V   a   COMMNEQ%RHS+MPI_CONSTANTS $   EG  b       OPNEQ+MPI_CONSTANTS (   �G  T   a   OPNEQ%LHS+MPI_CONSTANTS (   �G  T   a   OPNEQ%RHS+MPI_CONSTANTS -   OH  �       FILL_HALO_2D+MG_MPI_EXCHANGE 5   �H  @      FILL_HALO_2D%PRESENT+MG_MPI_EXCHANGE 2   (I  =      FILL_HALO_2D%TRIM+MG_MPI_EXCHANGE 1   eI  @   a   FILL_HALO_2D%LEV+MG_MPI_EXCHANGE 1   �I  �   a   FILL_HALO_2D%A2D+MG_MPI_EXCHANGE 6   IJ  P   a   FILL_HALO_2D%LBC_NULL+MG_MPI_EXCHANGE -   �J  �       FILL_HALO_3D+MG_MPI_EXCHANGE 5   GK  @      FILL_HALO_3D%PRESENT+MG_MPI_EXCHANGE 2   �K  =      FILL_HALO_3D%TRIM+MG_MPI_EXCHANGE 2   �K  =      FILL_HALO_3D%SIZE+MG_MPI_EXCHANGE 1   L  @   a   FILL_HALO_3D%LEV+MG_MPI_EXCHANGE /   AL  �   a   FILL_HALO_3D%P+MG_MPI_EXCHANGE 6   �L  P   a   FILL_HALO_3D%LBC_NULL+MG_MPI_EXCHANGE 3   MM  �       FILL_HALO_3D_RELAX+MG_MPI_EXCHANGE 8   �M  =      FILL_HALO_3D_RELAX%SIZE+MG_MPI_EXCHANGE 7   N  @   a   FILL_HALO_3D_RELAX%LEV+MG_MPI_EXCHANGE 5   WN  �  a   FILL_HALO_3D_RELAX%P+MG_MPI_EXCHANGE 6   HQ  @   a   FILL_HALO_3D_RELAX%NX+MG_MPI_EXCHANGE 6   �Q  @   a   FILL_HALO_3D_RELAX%NY+MG_MPI_EXCHANGE 6   �Q  @   a   FILL_HALO_3D_RELAX%NZ+MG_MPI_EXCHANGE -   R  p       FILL_HALO_4D+MG_MPI_EXCHANGE 2   xR  =      FILL_HALO_4D%SIZE+MG_MPI_EXCHANGE 1   �R  @   a   FILL_HALO_4D%LEV+MG_MPI_EXCHANGE 0   �R  �   a   FILL_HALO_4D%CA+MG_MPI_EXCHANGE $   �S  _       GATHER_2D+MG_GATHER (   (T  @   a   GATHER_2D%LEV+MG_GATHER &   hT  �   a   GATHER_2D%X+MG_GATHER &   U  �   a   GATHER_2D%Y+MG_GATHER $   �U  _       GATHER_3D+MG_GATHER (   V  @   a   GATHER_3D%LEV+MG_GATHER &   OV  �   a   GATHER_3D%X+MG_GATHER &   W  �   a   GATHER_3D%Y+MG_GATHER    �W  Q       FINE2COARSE     X  @   a   FINE2COARSE%LEV    XX         FINE2COARSE_2D $   �X  =      FINE2COARSE_2D%SIZE !   Y  �   a   FINE2COARSE_2D%X !   �Y  �   a   FINE2COARSE_2D%Y "   �Z  @   a   FINE2COARSE_2D%NX "   �Z  @   a   FINE2COARSE_2D%NY    [  n       FINE2COARSE_3D !   z[  �   a   FINE2COARSE_3D%X !   6\  �   a   FINE2COARSE_3D%Y "   �\  @   a   FINE2COARSE_3D%NX "   2]  @   a   FINE2COARSE_3D%NY "   r]  @   a   FINE2COARSE_3D%NZ    �]  Q       COARSE2FINE     ^  @   a   COARSE2FINE%LEV &   C^  h       COARSE2FINE_2D_LINEAR )   �^  �   a   COARSE2FINE_2D_LINEAR%XF )   g_  �   a   COARSE2FINE_2D_LINEAR%XC )   #`  @   a   COARSE2FINE_2D_LINEAR%NX )   c`  @   a   COARSE2FINE_2D_LINEAR%NY &   �`  �       COARSE2FINE_3D_LINEAR *   2a  <      COARSE2FINE_3D_LINEAR%MOD )   na  �   a   COARSE2FINE_3D_LINEAR%XF )   *b  �   a   COARSE2FINE_3D_LINEAR%XC )   �b  @   a   COARSE2FINE_3D_LINEAR%NX )   &c  @   a   COARSE2FINE_3D_LINEAR%NY )   fc  @   a   COARSE2FINE_3D_LINEAR%NZ %   �c  ]       MPI_OP+MPI_CONSTANTS -   d  H   a   MPI_OP%MPI_VAL+MPI_CONSTANTS '   Kd  ]       MPI_COMM+MPI_CONSTANTS /   �d  H   a   MPI_COMM%MPI_VAL+MPI_CONSTANTS -   �d  ]       MPI_ERRHANDLER+MPI_CONSTANTS 5   Me  H   a   MPI_ERRHANDLER%MPI_VAL+MPI_CONSTANTS '   �e  ]       MPI_FILE+MPI_CONSTANTS /   �e  H   a   MPI_FILE%MPI_VAL+MPI_CONSTANTS *   :f  ]       MPI_REQUEST+MPI_CONSTANTS 2   �f  H   a   MPI_REQUEST%MPI_VAL+MPI_CONSTANTS +   �f  ]       MPI_DATATYPE+MPI_CONSTANTS 3   <g  H   a   MPI_DATATYPE%MPI_VAL+MPI_CONSTANTS &   �g  ]       MPI_WIN+MPI_CONSTANTS .   �g  H   a   MPI_WIN%MPI_VAL+MPI_CONSTANTS *   )h  ]       MPI_MESSAGE+MPI_CONSTANTS 2   �h  H   a   MPI_MESSAGE%MPI_VAL+MPI_CONSTANTS '   �h  ]       MPI_INFO+MPI_CONSTANTS /   +i  H   a   MPI_INFO%MPI_VAL+MPI_CONSTANTS (   si  ]       MPI_GROUP+MPI_CONSTANTS 0   �i  H   a   MPI_GROUP%MPI_VAL+MPI_CONSTANTS 