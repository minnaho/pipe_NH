	  üX  ­   k820309    h          14.0        k?]                                                                                                           
       mg_horiz_grids.f90 MG_HORIZ_GRIDS                                                    
                                                          
                                                          
                                                          
                                                          
                                                          
                                                         
                         @                               '             .      #NX 	   #NY 
   #NZ    #NPX    #NPY    #INCX    #INCY    #NG2D    #NG    #NGP    #GATHER    #NGX    #NGY    #KEY    #LOCALCOMM    #NEIGHB    #CA    #P    #B    #R    #PX    #PY    #PZ    #DX     #DY !   #DXU "   #DYV #   #DZ $   #DZW %   #ARX &   #ARY '   #ARZ (   #BETA )   #GAMU *   #GAMV +   #ZXDY ,   #ZYDX -   #ALPHA .   #DUMMY3DNZ /   #DUMMY3DNZP 0   #DUMMY3 1   #GATHERBUFFER2D 2   #GATHERBUFFER 3   #DU 4   #DV 5   #DW 6                                              	                                                              
                                                                                                                                                                                                                                                                                                                                                                                                                                                        	                                                      $       
                                                      (                                                             ,                                                             0                                                             4                                                             8                                                                    <                   p          p            p                                                                            `                
            &                   &                   &                   &                                                                                 y
                                                                                                  ð                
            &                   &                   &                                                                                 y
                                                                                                  h               
            &                   &                   &                                                                                 y
                                                                                                  à               
            &                   &                   &                                                                                 y
                                                                                                  X               
            &                   &                   &                                                                                 y
                                                                                                  Ð               
            &                   &                   &                                                                                 y
                                                                                                  H               
            &                   &                   &                                                                                 y
                                                                                                   À               
            &                   &                                                                                 y
                                                                                      !                            
            &                   &                                                                                 y
                                                                                      "                           
            &                   &                                                                                 y
                                                                                      #            à               
            &                   &                                                                                 y
                                                                                      $            @               
            &                   &                   &                                                                                 y
                                                                                      %            ¸               
            &                   &                   &                                                                                 y
                                                                                      &            0               
            &                   &                   &                                                                                 y
                                                                                      '            ¨               
            &                   &                   &                                                                                 y
                                                                                      (                             
            &                   &                                                                                 y
                                                                                      )                         !  
            &                   &                                                                                 y
                                                                                      *            à             "  
            &                   &                                                                                 y
                                                                                      +            @             #  
            &                   &                                                                                 y
                                                                                      ,                          $  
            &                   &                   &                                                                                 y
                                                                                      -            	             %  
            &                   &                   &                                                                                 y
                                                                                      .            	             &  
            &                   &                   &                                                                                 y
                                                                                      /            
             '  
            &                   &                   &                                                                                 y
                                                                                      0            
             (  
            &                   &                   &                                                                                 y
                                                                                      1            ø
             )  
            &                   &                   &                                                                                 y
                                                                                      2            p             *  
            &                   &                   &                   &                                                                                 y
                                                                                      3                          +  
            &                   &                   &                   &                   &                                                                                 y
                                                                                      4            ¨             ,  
            &                   &                   &                                                                                 y
                                                                                      5                          -  
            &                   &                   &                                                                                 y
                                                                                      6                         .  
            &                   &                   &                                                                                 y
                                                                                           7                                                      4                                            8                                                      8                                            9                      @ P                               :                                   &                                           #GRID_TYPE                                                ;     
                
                       à?        0.5%         @                                <                           #LHS =   #RHS ?             
                                  =                   #MPI_GROUP >             
                                  ?                   #MPI_GROUP >   %         @                                @                           #LHS A   #RHS C             
                                  A                   #MPI_INFO B             
                                  C                   #MPI_INFO B   %         @                                D                           #LHS E   #RHS G             
                                  E                   #MPI_MESSAGE F             
                                  G                   #MPI_MESSAGE F   %         @                                H                           #LHS I   #RHS K             
                                  I                   #MPI_WIN J             
                                  K                   #MPI_WIN J   %         @                                L                           #LHS M   #RHS O             
                                  M                   #MPI_DATATYPE N             
                                  O                   #MPI_DATATYPE N   %         @                                P                           #LHS Q   #RHS S             
                                  Q                   #MPI_REQUEST R             
                                  S                   #MPI_REQUEST R   %         @                                T                           #LHS U   #RHS W             
                                  U                   #MPI_FILE V             
                                  W                   #MPI_FILE V   %         @                                X                           #LHS Y   #RHS [             
                                  Y                   #MPI_ERRHANDLER Z             
                                  [                   #MPI_ERRHANDLER Z   %         @                                \                           #LHS ]   #RHS _             
                                  ]                   #MPI_COMM ^             
                                  _                   #MPI_COMM ^   %         @                                `                           #LHS a   #RHS c             
                                  a                   #MPI_OP b             
                                  c                   #MPI_OP b   %         @                                d                           #LHS e   #RHS f             
                                  e                   #MPI_GROUP >             
                                  f                   #MPI_GROUP >   %         @                                g                           #LHS h   #RHS i             
                                  h                   #MPI_INFO B             
                                  i                   #MPI_INFO B   %         @                                j                           #LHS k   #RHS l             
                                  k                   #MPI_MESSAGE F             
                                  l                   #MPI_MESSAGE F   %         @                                m                           #LHS n   #RHS o             
                                  n                   #MPI_WIN J             
                                  o                   #MPI_WIN J   %         @                                p                           #LHS q   #RHS r             
                                  q                   #MPI_DATATYPE N             
                                  r                   #MPI_DATATYPE N   %         @                                s                           #LHS t   #RHS u             
                                  t                   #MPI_REQUEST R             
                                  u                   #MPI_REQUEST R   %         @                                v                           #LHS w   #RHS x             
                                  w                   #MPI_FILE V             
                                  x                   #MPI_FILE V   %         @                                y                           #LHS z   #RHS {             
                                  z                   #MPI_ERRHANDLER Z             
                                  {                   #MPI_ERRHANDLER Z   %         @                                |                           #LHS }   #RHS ~             
                                  }                   #MPI_COMM ^             
                                  ~                   #MPI_COMM ^   %         @                                                           #LHS    #RHS              
                                                     #MPI_OP b             
                                                     #MPI_OP b   #         @                                                     #FILL_HALO_2D%PRESENT    #FILL_HALO_2D%TRIM    #LEV    #A2D    #LBC_NULL                                                     PRESENT                                                  TRIM           
                                                     
                                                  
               &                   &                                                     
                                                           #         @                                                      #FILL_HALO_3D%PRESENT    #FILL_HALO_3D%TRIM    #FILL_HALO_3D%SIZE    #LEV    #P    #LBC_NULL                                                     PRESENT                                                  TRIM                                                  SIZE           
                                                     
                                                 
 )              &                   &                   &                                                     
                                                           #         @                                                      #FILL_HALO_3D_RELAX%SIZE    #LEV    #P    #NX    #NY    #NZ                                                     SIZE           
                                                     
                                                   
         p            5 O p        n                                           1p         p        p         5 O p        p          5 O p          & p          5 O p        n                                      1  & p          5 O p        n                                          1    5 O p             5 O p        n                                      1p         p             5 O p        n                                          1p         p                                                                      @                                                     @                                                     @                     #         @                                                      #FILL_HALO_4D%SIZE    #LEV    #CA                                                     SIZE           
                                                     
                                                 
 =              &                   &                   &                   &                                           #         @                                                      #LEV    #X    #Y              
                                                     
                                                   
              &                   &                                                                                                      
               &                   &                                           #         @                                                       #LEV    #X     #Y ¡             
                                                     
                                                    
              &                   &                   &                                                                                   ¡                   
               &                   &                   &                                           #         @                                   ¢                                      @                           b     '                    #MPI_VAL £                                               £                                    @                           ^     '                    #MPI_VAL ¤                                               ¤                                    @                           Z     '                    #MPI_VAL ¥                                               ¥                                    @                           V     '                    #MPI_VAL ¦                                               ¦                                    @                           R     '                    #MPI_VAL §                                               §                                    @                           N     '                    #MPI_VAL ¨                                               ¨                                    @                           J     '                    #MPI_VAL ©                                               ©                                    @                           F     '                    #MPI_VAL ª                                               ª                                    @                           B     '                    #MPI_VAL «                                               «                                    @                           >     '                    #MPI_VAL ¬                                               ¬                          *      fn#fn    Ê   @   J   MG_CST    
  @   J   MG_MPI    J  @   J   MG_TICTOC      @   J   MG_GRIDS     Ê  @   J   MG_MPI_EXCHANGE    
  @   J   MG_GATHER    J  @   J   MG_NAMELIST #           GRID_TYPE+MG_GRIDS &     H   a   GRID_TYPE%NX+MG_GRIDS &   æ  H   a   GRID_TYPE%NY+MG_GRIDS &   .  H   a   GRID_TYPE%NZ+MG_GRIDS '   v  H   a   GRID_TYPE%NPX+MG_GRIDS '   ¾  H   a   GRID_TYPE%NPY+MG_GRIDS (     H   a   GRID_TYPE%INCX+MG_GRIDS (   N  H   a   GRID_TYPE%INCY+MG_GRIDS (     H   a   GRID_TYPE%NG2D+MG_GRIDS &   Þ  H   a   GRID_TYPE%NG+MG_GRIDS '   &  H   a   GRID_TYPE%NGP+MG_GRIDS *   n  H   a   GRID_TYPE%GATHER+MG_GRIDS '   ¶  H   a   GRID_TYPE%NGX+MG_GRIDS '   þ  H   a   GRID_TYPE%NGY+MG_GRIDS '   F  H   a   GRID_TYPE%KEY+MG_GRIDS -     H   a   GRID_TYPE%LOCALCOMM+MG_GRIDS *   Ö     a   GRID_TYPE%NEIGHB+MG_GRIDS &   r	  <  a   GRID_TYPE%CA+MG_GRIDS %   ®
  $  a   GRID_TYPE%P+MG_GRIDS %   Ò  $  a   GRID_TYPE%B+MG_GRIDS %   ö  $  a   GRID_TYPE%R+MG_GRIDS &     $  a   GRID_TYPE%PX+MG_GRIDS &   >  $  a   GRID_TYPE%PY+MG_GRIDS &   b  $  a   GRID_TYPE%PZ+MG_GRIDS &       a   GRID_TYPE%DX+MG_GRIDS &       a   GRID_TYPE%DY+MG_GRIDS '       a   GRID_TYPE%DXU+MG_GRIDS '   ª    a   GRID_TYPE%DYV+MG_GRIDS &   ¶  $  a   GRID_TYPE%DZ+MG_GRIDS '   Ú  $  a   GRID_TYPE%DZW+MG_GRIDS '   þ  $  a   GRID_TYPE%ARX+MG_GRIDS '   "  $  a   GRID_TYPE%ARY+MG_GRIDS '   F    a   GRID_TYPE%ARZ+MG_GRIDS (   R    a   GRID_TYPE%BETA+MG_GRIDS (   ^    a   GRID_TYPE%GAMU+MG_GRIDS (   j    a   GRID_TYPE%GAMV+MG_GRIDS (   v  $  a   GRID_TYPE%ZXDY+MG_GRIDS (     $  a   GRID_TYPE%ZYDX+MG_GRIDS )   ¾   $  a   GRID_TYPE%ALPHA+MG_GRIDS -   â!  $  a   GRID_TYPE%DUMMY3DNZ+MG_GRIDS .   #  $  a   GRID_TYPE%DUMMY3DNZP+MG_GRIDS *   *$  $  a   GRID_TYPE%DUMMY3+MG_GRIDS 2   N%  <  a   GRID_TYPE%GATHERBUFFER2D+MG_GRIDS 0   &  T  a   GRID_TYPE%GATHERBUFFER+MG_GRIDS &   Þ'  $  a   GRID_TYPE%DU+MG_GRIDS &   )  $  a   GRID_TYPE%DV+MG_GRIDS &   &*  $  a   GRID_TYPE%DW+MG_GRIDS    J+  q       IP+MG_CST    »+  q       RP+MG_CST    ,,  @       NLEVS+MG_GRIDS    l,         GRID+MG_GRIDS    -  s       HLF+MG_CST &   z-  b       GROUPEQ+MPI_CONSTANTS *   Ü-  W   a   GROUPEQ%LHS+MPI_CONSTANTS *   3.  W   a   GROUPEQ%RHS+MPI_CONSTANTS %   .  b       INFOEQ+MPI_CONSTANTS )   ì.  V   a   INFOEQ%LHS+MPI_CONSTANTS )   B/  V   a   INFOEQ%RHS+MPI_CONSTANTS (   /  b       MESSAGEEQ+MPI_CONSTANTS ,   ú/  Y   a   MESSAGEEQ%LHS+MPI_CONSTANTS ,   S0  Y   a   MESSAGEEQ%RHS+MPI_CONSTANTS $   ¬0  b       WINEQ+MPI_CONSTANTS (   1  U   a   WINEQ%LHS+MPI_CONSTANTS (   c1  U   a   WINEQ%RHS+MPI_CONSTANTS )   ¸1  b       DATATYPEEQ+MPI_CONSTANTS -   2  Z   a   DATATYPEEQ%LHS+MPI_CONSTANTS -   t2  Z   a   DATATYPEEQ%RHS+MPI_CONSTANTS (   Î2  b       REQUESTEQ+MPI_CONSTANTS ,   03  Y   a   REQUESTEQ%LHS+MPI_CONSTANTS ,   3  Y   a   REQUESTEQ%RHS+MPI_CONSTANTS %   â3  b       FILEEQ+MPI_CONSTANTS )   D4  V   a   FILEEQ%LHS+MPI_CONSTANTS )   4  V   a   FILEEQ%RHS+MPI_CONSTANTS +   ð4  b       ERRHANDLEREQ+MPI_CONSTANTS /   R5  \   a   ERRHANDLEREQ%LHS+MPI_CONSTANTS /   ®5  \   a   ERRHANDLEREQ%RHS+MPI_CONSTANTS %   
6  b       COMMEQ+MPI_CONSTANTS )   l6  V   a   COMMEQ%LHS+MPI_CONSTANTS )   Â6  V   a   COMMEQ%RHS+MPI_CONSTANTS #   7  b       OPEQ+MPI_CONSTANTS '   z7  T   a   OPEQ%LHS+MPI_CONSTANTS '   Î7  T   a   OPEQ%RHS+MPI_CONSTANTS '   "8  b       GROUPNEQ+MPI_CONSTANTS +   8  W   a   GROUPNEQ%LHS+MPI_CONSTANTS +   Û8  W   a   GROUPNEQ%RHS+MPI_CONSTANTS &   29  b       INFONEQ+MPI_CONSTANTS *   9  V   a   INFONEQ%LHS+MPI_CONSTANTS *   ê9  V   a   INFONEQ%RHS+MPI_CONSTANTS )   @:  b       MESSAGENEQ+MPI_CONSTANTS -   ¢:  Y   a   MESSAGENEQ%LHS+MPI_CONSTANTS -   û:  Y   a   MESSAGENEQ%RHS+MPI_CONSTANTS %   T;  b       WINNEQ+MPI_CONSTANTS )   ¶;  U   a   WINNEQ%LHS+MPI_CONSTANTS )   <  U   a   WINNEQ%RHS+MPI_CONSTANTS *   `<  b       DATATYPENEQ+MPI_CONSTANTS .   Â<  Z   a   DATATYPENEQ%LHS+MPI_CONSTANTS .   =  Z   a   DATATYPENEQ%RHS+MPI_CONSTANTS )   v=  b       REQUESTNEQ+MPI_CONSTANTS -   Ø=  Y   a   REQUESTNEQ%LHS+MPI_CONSTANTS -   1>  Y   a   REQUESTNEQ%RHS+MPI_CONSTANTS &   >  b       FILENEQ+MPI_CONSTANTS *   ì>  V   a   FILENEQ%LHS+MPI_CONSTANTS *   B?  V   a   FILENEQ%RHS+MPI_CONSTANTS ,   ?  b       ERRHANDLERNEQ+MPI_CONSTANTS 0   ú?  \   a   ERRHANDLERNEQ%LHS+MPI_CONSTANTS 0   V@  \   a   ERRHANDLERNEQ%RHS+MPI_CONSTANTS &   ²@  b       COMMNEQ+MPI_CONSTANTS *   A  V   a   COMMNEQ%LHS+MPI_CONSTANTS *   jA  V   a   COMMNEQ%RHS+MPI_CONSTANTS $   ÀA  b       OPNEQ+MPI_CONSTANTS (   "B  T   a   OPNEQ%LHS+MPI_CONSTANTS (   vB  T   a   OPNEQ%RHS+MPI_CONSTANTS -   ÊB         FILL_HALO_2D+MG_MPI_EXCHANGE 5   cC  @      FILL_HALO_2D%PRESENT+MG_MPI_EXCHANGE 2   £C  =      FILL_HALO_2D%TRIM+MG_MPI_EXCHANGE 1   àC  @   a   FILL_HALO_2D%LEV+MG_MPI_EXCHANGE 1    D  ¤   a   FILL_HALO_2D%A2D+MG_MPI_EXCHANGE 6   ÄD  P   a   FILL_HALO_2D%LBC_NULL+MG_MPI_EXCHANGE -   E  ®       FILL_HALO_3D+MG_MPI_EXCHANGE 5   ÂE  @      FILL_HALO_3D%PRESENT+MG_MPI_EXCHANGE 2   F  =      FILL_HALO_3D%TRIM+MG_MPI_EXCHANGE 2   ?F  =      FILL_HALO_3D%SIZE+MG_MPI_EXCHANGE 1   |F  @   a   FILL_HALO_3D%LEV+MG_MPI_EXCHANGE /   ¼F  ¼   a   FILL_HALO_3D%P+MG_MPI_EXCHANGE 6   xG  P   a   FILL_HALO_3D%LBC_NULL+MG_MPI_EXCHANGE 3   ÈG         FILL_HALO_3D_RELAX+MG_MPI_EXCHANGE 8   UH  =      FILL_HALO_3D_RELAX%SIZE+MG_MPI_EXCHANGE 7   H  @   a   FILL_HALO_3D_RELAX%LEV+MG_MPI_EXCHANGE 5   ÒH  ñ  a   FILL_HALO_3D_RELAX%P+MG_MPI_EXCHANGE 6   ÃK  @   a   FILL_HALO_3D_RELAX%NX+MG_MPI_EXCHANGE 6   L  @   a   FILL_HALO_3D_RELAX%NY+MG_MPI_EXCHANGE 6   CL  @   a   FILL_HALO_3D_RELAX%NZ+MG_MPI_EXCHANGE -   L  p       FILL_HALO_4D+MG_MPI_EXCHANGE 2   óL  =      FILL_HALO_4D%SIZE+MG_MPI_EXCHANGE 1   0M  @   a   FILL_HALO_4D%LEV+MG_MPI_EXCHANGE 0   pM  Ô   a   FILL_HALO_4D%CA+MG_MPI_EXCHANGE $   DN  _       GATHER_2D+MG_GATHER (   £N  @   a   GATHER_2D%LEV+MG_GATHER &   ãN  ¤   a   GATHER_2D%X+MG_GATHER &   O  ¤   a   GATHER_2D%Y+MG_GATHER $   +P  _       GATHER_3D+MG_GATHER (   P  @   a   GATHER_3D%LEV+MG_GATHER &   ÊP  ¼   a   GATHER_3D%X+MG_GATHER &   Q  ¼   a   GATHER_3D%Y+MG_GATHER     BR  H       SET_HORIZ_GRIDS %   R  ]       MPI_OP+MPI_CONSTANTS -   çR  H   a   MPI_OP%MPI_VAL+MPI_CONSTANTS '   /S  ]       MPI_COMM+MPI_CONSTANTS /   S  H   a   MPI_COMM%MPI_VAL+MPI_CONSTANTS -   ÔS  ]       MPI_ERRHANDLER+MPI_CONSTANTS 5   1T  H   a   MPI_ERRHANDLER%MPI_VAL+MPI_CONSTANTS '   yT  ]       MPI_FILE+MPI_CONSTANTS /   ÖT  H   a   MPI_FILE%MPI_VAL+MPI_CONSTANTS *   U  ]       MPI_REQUEST+MPI_CONSTANTS 2   {U  H   a   MPI_REQUEST%MPI_VAL+MPI_CONSTANTS +   ÃU  ]       MPI_DATATYPE+MPI_CONSTANTS 3    V  H   a   MPI_DATATYPE%MPI_VAL+MPI_CONSTANTS &   hV  ]       MPI_WIN+MPI_CONSTANTS .   ÅV  H   a   MPI_WIN%MPI_VAL+MPI_CONSTANTS *   W  ]       MPI_MESSAGE+MPI_CONSTANTS 2   jW  H   a   MPI_MESSAGE%MPI_VAL+MPI_CONSTANTS '   ²W  ]       MPI_INFO+MPI_CONSTANTS /   X  H   a   MPI_INFO%MPI_VAL+MPI_CONSTANTS (   WX  ]       MPI_GROUP+MPI_CONSTANTS 0   ´X  H   a   MPI_GROUP%MPI_VAL+MPI_CONSTANTS 