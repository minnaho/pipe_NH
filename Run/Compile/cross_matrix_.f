      implicit none
      integer*4 max_name_size, max_names, max_string_size, input,iout
      parameter (max_name_size=62, max_names=956, max_string_size=92,
     &                                             input=11, iout=12)
      character*62 fname(max_names)
      character*62 modname
      character string*92, quote*1, double_quote*1, backslash*1
      integer*4 nsize(max_names), size,  lines, nmax,  n_disc,  n_excl,
     &                              last_arg, iargc, iocheck, i,j,k,n
      logical matrix(max_names,max_names), new_name, not_end_of_file
      character ctoupper,tab
      tab=char(9)
      last_arg=iargc()
      if (last_arg.eq.0) then
        write(*,*) 'ERROR IN cross_matrix: NO FILE NAMES ARE GIVEN'
        stop
      endif
      write(*,'(/1x,A,1x,A/)') 'This is CROSS_MATRIX: Creating',
     &                             'new version of Make.depend.'
      quote=char(39)
      double_quote=char(34)
      backslash=char(92)
      lines=0
      n_disc=0
      n_excl=0
      do j=1,max_names
        nsize(j)=0
        do i=1,max_name_size
          fname(j)(i:i)=' '
        enddo
        do i=1,max_names
          matrix(i,j)=.false.
        enddo
      enddo
      do n=1,last_arg
        call getarg(n,fname(n))
        i=max_name_size+1
  1     i=i-1
         if (fname(n)(i:i).eq.' ' .and. i.gt.1) goto 1
        nsize(n)=i
      enddo
      nmax=last_arg
      open(unit=iout, file='Make.depend', form='formatted')
      write(iout,'(3(A,1x,A/),A/A/A/A/A/A)')
     &         '# Make.depend: list of dependencies generated by',
     &         'cross_matrix.', '# WARNING: THIS IS A MACHINE',
     &         'GENERATED FILE: DO NOT EDIT !!!!!', '# To create',
     &         'or update this file use commands:', '#',
     &         '#        cross_matrix *.F',         '# or',
     &         '#        cross_matrix *.F *.h',     '# or',
     &         '#        make depend'
      n=0
  2    n=n+1
        not_end_of_file=.true.
      IF (n.GT.nmax) GOTO 29
      IF (fname(n)(nsize(n)-3:nsize(n)).EQ.'.mod') goto 2
        open(unit=input, file=fname(n), form='formatted',
     &                               status='old', err=3)
        goto 4
  3     continue
        n_excl=n_excl+1
        write(iout,'(A/A,2x,A/A)') '#',
     &         '# WARNING: File is not found:',fname(n)(1:nsize(n)),
     &         '# This file is excluded from the dependency list.'
        do i=1,nsize(n)
          fname(n)(i:i)=' '
        enddo
        nsize(n)=0
        goto 2
  4     string(1:1)=' '
         read(input,'(A)',iostat=iocheck,end=5) string
         lines=lines+1
        if (iocheck.eq.0 .and. string(1:1).ne.'#') goto 19
        goto 6
  5     not_end_of_file=.false.
  6     i=2
  7      if (string(i:i).ne.'i') then
           i=i+1
           if (i.lt.max_string_size) goto 7
         elseif (string(i:i+6) .eq. 'include') then
           i=i+7
  8        if (string(i:i).eq.double_quote) then
             j=i+2
  9          if (string(j:j).eq.double_quote) then
               size=j-i-1
               new_name=.true.
               do k=1,nmax
                 if (size .eq. nsize(k)) then
                   if (string(i+1:j-1) .eq. fname(k)(1:size)) then
                     new_name=.false.
                     matrix(k,n)=.true.
                   endif
                 endif
               enddo
               if (new_name) then
                 n_disc=n_disc+1
                 nmax=nmax+1
                 matrix(nmax,n)=.true.
                 nsize(nmax)=size
                 fname(nmax)(1:size)=string(i+1:j-1)
                 do i=size+1,max_name_size
                   fname(nmax)(i:i)=' '
                 enddo
               endif
             elseif (j.lt.max_string_size) then
               j=j+1
               goto 9
             endif
           elseif (i.lt.max_string_size) then
             i=i+1
             goto 8
           endif
         endif
 19      continue
       if (iocheck.eq.0 .and.
     &     ((string(1:1).ne.' ').AND.
     &     (string(1:1).ne.tab))) goto 4
         i=2
 17      if ((string(i:i).eq.' ').OR.(string(i:i).eq.tab)) then
        i=i+1
        if (i.lt.max_string_size) goto 17
      elseif ((string(i:i+3) .eq. 'USE ').or.
     &                  (string(i:i+3) .eq. 'use ')) then
        i=i+4
        j=i
 18        j=j+1
        if ((string(j:j).EQ.' ').or.(string(j:j).EQ.tab)) then
        j=j-1
        size=j-i+1
      do k=1,size
      modname(k:k)=ctoupper(string(i+k-1:i+k-1))
      end do
      modname(size+1:size+4)='.mod'
        size=size+4
              new_name=.true.
               do k=1,nmax
                 if (size .eq. nsize(k)) then
                   if (modname(1:size) .eq. fname(k)(1:size)) then
                     new_name=.false.
                     matrix(k,n)=.true.
                   endif
                 endif
               enddo
               if (new_name) then
                 n_disc=n_disc+1
                 nmax=nmax+1
                 matrix(nmax,n)=.true.
                 nsize(nmax)=size
                 fname(nmax)(1:size)=modname(1:size)
                 do i=size+1,max_name_size
                   fname(nmax)(i:i)=' '
                 enddo
               endif
          elseif (j.lt.max_string_size) then
               goto 18
          end if
      end if
        if (not_end_of_file) goto 4
       close (unit=input)
      if (n.lt.nmax) goto 2
 29      continue
      i=0
  10  do n=1,nmax
        do k=1,nmax
          if (matrix(k,n)) then
            do j=1,nmax
              if (matrix(j,k))
     &         matrix(j,n)=.true.
            enddo
          endif
        enddo
      enddo
      j=0
      do n=1,nmax
        do k=1,nmax
          if (matrix(k,n)) j=j+1
        enddo
      enddo
      if (i.ne.j) then
        i=j
        goto 10
      endif
      write(iout,'(A/A,5x,I4/A/A,18x,I4/A/A,1x,I4)') '#',
     &  '# Number of files given for dependency analysis:', last_arg,
     &  '#', '# Number of newly discovered files:', n_disc, '#',
     &  '# Number of files excluded from dependency analysis:',n_excl
      write(iout,'(A/A,3x,I4/A/A,9x,I6/A)') '#',
     &  '# Total number of files analyzed for dependencies:', nmax,
     &  '#', '# Total number of code lines in all files:', lines, '#'
      do n=1,last_arg
        if (nsize(n).gt.0) then
          write(iout,'(A1)') '#'
          k=0
  11      i=nsize(n)
           string(1:i)=fname(n)(1:i)
           if (string(i-1:i).eq.'.F') string(i:i)='o'
           i=i+1
           string(i:i)=':'
  12       i=i+1
            string(i:i)=' '
            if (k.eq.0) then
              if (fname(n)(nsize(n)-1:nsize(n)).eq.'.F') then
                string(i+1:i+nsize(n))=fname(n)(1:nsize(n))
                i=i+nsize(n)+1
                string(i:i)=' '
              endif
            endif
  13        k=k+1
            if (matrix(k,n)) then
              if (i+nsize(k).lt.max_string_size) then
                string(i+1:i+nsize(k))=fname(k)(1:nsize(k))
                i=i+nsize(k)
                goto 12
              else
                write(iout,'(A)') string(1:i)
                goto 11
              endif
            elseif (k.lt.nmax) then
              goto 13
            else
              write(iout,'(A)') string(1:i)
              if (i.gt.nsize(n)+2) then
                write(iout,'(A)') string(1:nsize(n)+1)
              endif
            endif
        endif
      enddo
      close (iout)
      stop
      end
      character function ctoupper(c)
      character c
      character temp
      temp=c
      IF ((ichar(c).GE.ichar('a')).AND.(ichar(c).LE.ichar('z')))
     &  temp=char(ichar(c)+ichar('A')-ichar('a'))
      ctoupper=temp
      return
      end