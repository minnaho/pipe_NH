	  �!  b   k820309    ?          14.0        �qZ                                                                                                           
       mg_mpi.f90 MG_MPI                                                    
                                                          
                                                                                                         4#         @                                       	               #COMM    #RANK    #IERROR                                                                                                                                                                                                                                                   
          D            1140850688#         @                                  	     	               #COMM 
   #SIZE    #IERROR                                              
                                                                                                                        %         @                                                           #LHS    #RHS              
                                                     #MPI_GROUP              
                                                     #MPI_GROUP    %         @                                                           #LHS    #RHS              
                                                     #MPI_INFO              
                                                     #MPI_INFO    %         @                                                           #LHS    #RHS              
                                                     #MPI_MESSAGE              
                                                     #MPI_MESSAGE    %         @                                                           #LHS    #RHS              
                                                     #MPI_WIN              
                                                     #MPI_WIN    %         @                                                           #LHS    #RHS               
                                                     #MPI_DATATYPE              
                                                      #MPI_DATATYPE    %         @                                !                           #LHS "   #RHS $             
                                  "                   #MPI_REQUEST #             
                                  $                   #MPI_REQUEST #   %         @                                %                           #LHS &   #RHS (             
                                  &                   #MPI_FILE '             
                                  (                   #MPI_FILE '   %         @                                )                           #LHS *   #RHS ,             
                                  *                   #MPI_ERRHANDLER +             
                                  ,                   #MPI_ERRHANDLER +   %         @                                -                           #LHS .   #RHS 0             
                                  .                   #MPI_COMM /             
                                  0                   #MPI_COMM /   %         @                                1                           #LHS 2   #RHS 4             
                                  2                   #MPI_OP 3             
                                  4                   #MPI_OP 3   %         @                                5                           #LHS 6   #RHS 7             
                                  6                   #MPI_GROUP              
                                  7                   #MPI_GROUP    %         @                                8                           #LHS 9   #RHS :             
                                  9                   #MPI_INFO              
                                  :                   #MPI_INFO    %         @                                ;                           #LHS <   #RHS =             
                                  <                   #MPI_MESSAGE              
                                  =                   #MPI_MESSAGE    %         @                                >                           #LHS ?   #RHS @             
                                  ?                   #MPI_WIN              
                                  @                   #MPI_WIN    %         @                                A                           #LHS B   #RHS C             
                                  B                   #MPI_DATATYPE              
                                  C                   #MPI_DATATYPE    %         @                                D                           #LHS E   #RHS F             
                                  E                   #MPI_REQUEST #             
                                  F                   #MPI_REQUEST #   %         @                                G                           #LHS H   #RHS I             
                                  H                   #MPI_FILE '             
                                  I                   #MPI_FILE '   %         @                                J                           #LHS K   #RHS L             
                                  K                   #MPI_ERRHANDLER +             
                                  L                   #MPI_ERRHANDLER +   %         @                                M                           #LHS N   #RHS O             
                                  N                   #MPI_COMM /             
                                  O                   #MPI_COMM /   %         @                                P                           #LHS Q   #RHS R             
                                  Q                   #MPI_OP 3             
                                  R                   #MPI_OP 3              @ @                              S                       @ @                              T            #         @                                   U                     #         @                                  V                     #         @                                  W                                      @                           3     '                    #MPI_VAL X                �                               X                                    @                           /     '                    #MPI_VAL Y                �                               Y                                    @                           +     '                    #MPI_VAL Z                �                               Z                                    @                           '     '                    #MPI_VAL [                �                               [                                    @                           #     '                    #MPI_VAL \                �                               \                                    @                                '                    #MPI_VAL ]                �                               ]                                    @                                '                    #MPI_VAL ^                �                               ^                                    @                                '                    #MPI_VAL _                �                               _                                    @                                '                    #MPI_VAL `                �                               `                                    @                                '                    #MPI_VAL a                �                               a                      �         fn#fn    �   @   J   MG_CST    �   @   J   MPI    :  q       IP+MG_CST '   �  h       MPI_COMM_RANK+MPI_BASE ,     @   a   MPI_COMM_RANK%COMM+MPI_BASE ,   S  @   a   MPI_COMM_RANK%RANK+MPI_BASE .   �  @   a   MPI_COMM_RANK%IERROR+MPI_BASE -   �  z       MPI_COMM_WORLD+MPI_CONSTANTS '   M  h       MPI_COMM_SIZE+MPI_BASE ,   �  @   a   MPI_COMM_SIZE%COMM+MPI_BASE ,   �  @   a   MPI_COMM_SIZE%SIZE+MPI_BASE .   5  @   a   MPI_COMM_SIZE%IERROR+MPI_BASE &   u  b       GROUPEQ+MPI_CONSTANTS *   �  W   a   GROUPEQ%LHS+MPI_CONSTANTS *   .  W   a   GROUPEQ%RHS+MPI_CONSTANTS %   �  b       INFOEQ+MPI_CONSTANTS )   �  V   a   INFOEQ%LHS+MPI_CONSTANTS )   =  V   a   INFOEQ%RHS+MPI_CONSTANTS (   �  b       MESSAGEEQ+MPI_CONSTANTS ,   �  Y   a   MESSAGEEQ%LHS+MPI_CONSTANTS ,   N  Y   a   MESSAGEEQ%RHS+MPI_CONSTANTS $   �  b       WINEQ+MPI_CONSTANTS (   	  U   a   WINEQ%LHS+MPI_CONSTANTS (   ^  U   a   WINEQ%RHS+MPI_CONSTANTS )   �  b       DATATYPEEQ+MPI_CONSTANTS -   	  Z   a   DATATYPEEQ%LHS+MPI_CONSTANTS -   o	  Z   a   DATATYPEEQ%RHS+MPI_CONSTANTS (   �	  b       REQUESTEQ+MPI_CONSTANTS ,   +
  Y   a   REQUESTEQ%LHS+MPI_CONSTANTS ,   �
  Y   a   REQUESTEQ%RHS+MPI_CONSTANTS %   �
  b       FILEEQ+MPI_CONSTANTS )   ?  V   a   FILEEQ%LHS+MPI_CONSTANTS )   �  V   a   FILEEQ%RHS+MPI_CONSTANTS +   �  b       ERRHANDLEREQ+MPI_CONSTANTS /   M  \   a   ERRHANDLEREQ%LHS+MPI_CONSTANTS /   �  \   a   ERRHANDLEREQ%RHS+MPI_CONSTANTS %     b       COMMEQ+MPI_CONSTANTS )   g  V   a   COMMEQ%LHS+MPI_CONSTANTS )   �  V   a   COMMEQ%RHS+MPI_CONSTANTS #     b       OPEQ+MPI_CONSTANTS '   u  T   a   OPEQ%LHS+MPI_CONSTANTS '   �  T   a   OPEQ%RHS+MPI_CONSTANTS '     b       GROUPNEQ+MPI_CONSTANTS +     W   a   GROUPNEQ%LHS+MPI_CONSTANTS +   �  W   a   GROUPNEQ%RHS+MPI_CONSTANTS &   -  b       INFONEQ+MPI_CONSTANTS *   �  V   a   INFONEQ%LHS+MPI_CONSTANTS *   �  V   a   INFONEQ%RHS+MPI_CONSTANTS )   ;  b       MESSAGENEQ+MPI_CONSTANTS -   �  Y   a   MESSAGENEQ%LHS+MPI_CONSTANTS -   �  Y   a   MESSAGENEQ%RHS+MPI_CONSTANTS %   O  b       WINNEQ+MPI_CONSTANTS )   �  U   a   WINNEQ%LHS+MPI_CONSTANTS )     U   a   WINNEQ%RHS+MPI_CONSTANTS *   [  b       DATATYPENEQ+MPI_CONSTANTS .   �  Z   a   DATATYPENEQ%LHS+MPI_CONSTANTS .     Z   a   DATATYPENEQ%RHS+MPI_CONSTANTS )   q  b       REQUESTNEQ+MPI_CONSTANTS -   �  Y   a   REQUESTNEQ%LHS+MPI_CONSTANTS -   ,  Y   a   REQUESTNEQ%RHS+MPI_CONSTANTS &   �  b       FILENEQ+MPI_CONSTANTS *   �  V   a   FILENEQ%LHS+MPI_CONSTANTS *   =  V   a   FILENEQ%RHS+MPI_CONSTANTS ,   �  b       ERRHANDLERNEQ+MPI_CONSTANTS 0   �  \   a   ERRHANDLERNEQ%LHS+MPI_CONSTANTS 0   Q  \   a   ERRHANDLERNEQ%RHS+MPI_CONSTANTS &   �  b       COMMNEQ+MPI_CONSTANTS *     V   a   COMMNEQ%LHS+MPI_CONSTANTS *   e  V   a   COMMNEQ%RHS+MPI_CONSTANTS $   �  b       OPNEQ+MPI_CONSTANTS (     T   a   OPNEQ%LHS+MPI_CONSTANTS (   q  T   a   OPNEQ%RHS+MPI_CONSTANTS    �  @       MYRANK      @       NPROCS    E  H       MG_MPI_INIT    �  H       MPI_MYRANK    �  H       MPI_NPROCS %     ]       MPI_OP+MPI_CONSTANTS -   z  H   a   MPI_OP%MPI_VAL+MPI_CONSTANTS '   �  ]       MPI_COMM+MPI_CONSTANTS /     H   a   MPI_COMM%MPI_VAL+MPI_CONSTANTS -   g  ]       MPI_ERRHANDLER+MPI_CONSTANTS 5   �  H   a   MPI_ERRHANDLER%MPI_VAL+MPI_CONSTANTS '     ]       MPI_FILE+MPI_CONSTANTS /   i  H   a   MPI_FILE%MPI_VAL+MPI_CONSTANTS *   �  ]       MPI_REQUEST+MPI_CONSTANTS 2     H   a   MPI_REQUEST%MPI_VAL+MPI_CONSTANTS +   V  ]       MPI_DATATYPE+MPI_CONSTANTS 3   �  H   a   MPI_DATATYPE%MPI_VAL+MPI_CONSTANTS &   �  ]       MPI_WIN+MPI_CONSTANTS .   X  H   a   MPI_WIN%MPI_VAL+MPI_CONSTANTS *   �  ]       MPI_MESSAGE+MPI_CONSTANTS 2   �  H   a   MPI_MESSAGE%MPI_VAL+MPI_CONSTANTS '   E   ]       MPI_INFO+MPI_CONSTANTS /   �   H   a   MPI_INFO%MPI_VAL+MPI_CONSTANTS (   �   ]       MPI_GROUP+MPI_CONSTANTS 0   G!  H   a   MPI_GROUP%MPI_VAL+MPI_CONSTANTS 