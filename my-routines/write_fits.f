!     subroutines to write the fits file
c-----------------------------------------------------------------------
      ! create fits file with one extension and  columns
      subroutine create_fits(filename)
      implicit none
      character* (*) filename
!     Internal variables
      integer status,unit,blocksize
      logical simple,extend
      integer bitpix,naxis,naxes(1)
!Don't touch!
c Initialize status
      status=0
c Delete the file if it already exists
      call deletefile(filename,status)
      if( status .gt. 0 )call printerror(status)
c Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)
c Create the new empty FITS file
      blocksize=1
      call ftinit(unit,filename,blocksize,status)
c Initialize parameters about the FITS primary image - need image for structure
      simple=.true.
      bitpix=16
      naxis=1
      naxes(1)=0
      extend=.true.
c Write the required header keywords to the file
      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
c The FITS file must always be closed before exiting the program.
c Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)
c Check for any error, and if so print out error messages.
      if (status .gt. 0) call printerror(status)
      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      ! create fits file with one extension and  columns
      subroutine write_param(temp_dim, en_dim, temp, en, filename)
      implicit none
      integer          :: temp_dim, en_dim
      double precision :: temp(temp_dim), en(en_dim)
      character* (*) filename
!     Internal variables
      integer status,unit,blocksize
      integer rownum,readwrite,colnum,tfields
      character (len=16) extname, name
      integer,parameter :: coldim = 2
      character (LEN=16) ttype(coldim),tform(coldim),tunit(coldim)
      integer nrows, varidat
c$$$! inverse of kbol (K * ev-1)
c$$$      double precision, parameter :: ikbol   = 1.16d4
c$$$! m_e c^2 (eV)
c$$$  double precision, parameter :: mec2  = 5.11d5

c Initialize status
      status=0
      blocksize = 1
c Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)
c Open the FITS file, with write access.
      readwrite=1
      call ftopen(unit,filename,readwrite,blocksize,status)
c Append/create a new empty HDU onto the end of the file and move to it.
      call ftcrhd(unit,status)
      ttype(1) = 'TEMP'
      ttype(2) = 'ENERGIES'
      write(name, '(I4,A1)') temp_dim, 'D'
      tform(1) = trim(name)

      write(name, '(I4,A1)') en_dim, 'D'
      tform(2) = trim(name)
      write(*,*) tform(1) , tform(2)
c$$$      tform(1) = '70D'
c$$$      tform(2) = '500D'
      tunit(1) = 'kT/mec2'
      tunit(2) = 'eV'

c Define parameters for the binary table (see the above data statements)
      tfields  = coldim
      nrows    = 1
      extname  = 'PARAMETERS'
      varidat  = 0
c FTPHBN writes all the required header keywords which define the
c structure of the binary table. NROWS and TFIELDS gives the number of
c rows and columns in the table, and the TTYPE, TFORM, and TUNIT arrays
c give the column name, format, and units, respectively of each column.
      call ftphbn(unit,nrows,tfields,ttype,tform,tunit,
     &     extname,varidat,status)
!     Filling the columns
            rownum = 1
            colnum = 1
            call ftpcld(unit, colnum, rownum, 1, temp_dim,
     &        temp, status)
            colnum = 2
            call ftpcld(unit, colnum, rownum, 1, en_dim,
     &       en, status)
c$$$c Write keywords to this extension
c$$$!ftpkyj to write integer
c$$$!ftpkye to write real (4)
c$$$!ftpkyd to write real (8)
c The FITS file must always be closed before exiting the program.
c Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)
c Check for any error, and if so print out error messages.
      if (status .gt. 0) call printerror(status)
      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
! Add one HDU in the 3rd position
!IMPORTANT: this routine leaves the fits file opened

      subroutine add_HDU(temp_dim, en_dim, filename, unit)
      implicit none
      integer :: unit, temp_dim, en_dim
      character* (*) filename

!     Internal variables
      integer status, blocksize, readwrite, tfields
      character (len=16) extname
      integer,parameter :: coldim = 4
      character (LEN=16) ttype(coldim),tform(coldim),tunit(coldim)
      integer i, nrows, varidat
      ttype(1) = 'KN_CROSS'
      ttype(2) = 'OUT_IND'
      ttype(3) = 'OUT_NUM'
      ttype(4) = 'SRF'

      tform(1) = '1D'
      tform(2) = '1I'
      tform(3) = '1I'
      tform(4) = '1PD'
      do i = 1, coldim
         tunit(i) = ''
      enddo
c Initialize status
      status=0
      blocksize = 1
c Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)
c Open the FITS file, with write access.
      readwrite=1
      call ftopen(unit,filename,readwrite,blocksize,status)
c Append/create a new empty HDU onto the end of the file and move to it.
      call ftcrhd(unit,status)

c Define parameters for the binary table (see the above data statements)
      tfields  = coldim
      nrows    = temp_dim * en_dim
      extname  = 'SRF'
      varidat  = 0
c FTPHBN writes all the required header keywords which define the
c structure of the binary table. NROWS and TFIELDS gives the number of
c rows and columns in the table, and the TTYPE, TFORM, and TUNIT arrays
c give the column name, format, and units, respectively of each column.
      call ftphbn(unit,nrows,tfields,ttype,tform,tunit,
     &     extname,varidat,status)
c$$$c Write keywords to this extension
c$$$!ftpkyj to write integer
c$$$!ftpkye to write real (4)
c$$$!ftpkyd to write real (8)
c$$$c The FITS file must always be closed before exiting the program.
c$$$c Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
c$$$      call ftclos(unit, status)
c$$$      call ftfiou(unit, status)
c Check for any error, and if so print out error messages.
      if (status .gt. 0) then
         call printerror(status)
      endif
      return
      end
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
      subroutine add_row_HDU(n, nmaxp, out_en_dim, out_en_ind,
     & kn_cross, srf, unit)
      implicit none

      integer          :: unit, n, nmaxp, out_en_dim, out_en_ind
      double precision :: srf(nmaxp)
      double precision :: kn_cross(1)
      integer status,readwrite,hdutype, blocksize
      integer colnum,rownum

c Initialize status
      status=0
c$$$      blocksize = 1
c$$$c Get an unused Logical Unit Number to use to open the FITS file.
c$$$      call ftgiou(unit,status)
c$$$c Open the FITS file, with write access.
c$$$      readwrite=1
c$$$      call ftopen(unit,filename,readwrite,blocksize,status)
C  Move to the last (3nd) HDU in the file (the paameter values table).
      call ftmahd(unit,3,hdutype,status)
!     Filling the columns
            rownum = n
            colnum = 1
            call ftpcld(unit, colnum, rownum, 1, 1,
     &        kn_cross(1), status)
            colnum = 2
            call ftpcli(unit, colnum, rownum, 1, 1,
     &       out_en_ind, status)
            colnum = 3
            call ftpcli(unit, colnum, rownum, 1, 1,
     &       out_en_dim, status)
            colnum = 4
!Fill the column with arrays of different size
            call ftpcld(unit, colnum, rownum, 1, out_en_dim,
     &           srf, status)
c$$$cThe FITS file must always be closed before exiting the program.
c$$$c Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
c$$$      call ftclos(unit, status)
c$$$      call ftfiou(unit, status)
c Check for any error, and if so print out error messages.
            if (status .gt. 0) then
               call printerror(status)
               write(*,*) 'adding row'
            endif
      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine deletefile(filename,status)
C  A simple little routine to delete a FITS file
      implicit none
      integer status,unit,blocksize
      character*(*) filename

C  Simply return if status is greater than zero
      if (status .gt. 0)return
C  Get an unused Logical Unit Number to use to open the FITS file
      call ftgiou(unit,status)
C  Try to open the file, to see if it exists
      call ftopen(unit,filename,1,blocksize,status)
      if (status .eq. 0)then
C         file was opened;  so now delete it
          call ftdelt(unit,status)
          !write(*,*)"Deleted a file"
      else if (status .eq. 103)then
C         file doesn't exist, so just reset status to zero and clear errors
          status=0
          !write(*,*)"Didn't delete a file"
          call ftcmsg
      else
C         there was some other error opening the file; delete the file anyway
          status=0
          call ftcmsg
          call ftdelt(unit,status)
      end if
C  Free the unit number for later reuse
      call ftfiou(unit, status)
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine printerror(status)
C  This subroutine prints out the descriptive text corresponding to the
C  error status value and prints out the contents of the internal
C  error message stack generated by FITSIO whenever an error occurs.
      integer status
      character errtext*30,errmessage*80
C  Check if status is OK (no error); if so, simply return
      if (status .le. 0)return
C  The FTGERR subroutine returns a descriptive 30-character text string that
C  corresponds to the integer error status number.  A complete list of all
C  the error numbers can be found in the back of the FITSIO User's Guide.
      call ftgerr(status,errtext)
      print *,'FITSIO Error Status =',status,': ',errtext
C  FITSIO usually generates an internal stack of error messages whenever
C  an error occurs.  These messages provide much more information on the
C  cause of the problem than can be provided by the single integer error
C  status value.  The FTGMSG subroutine retrieves the oldest message from
C  the stack and shifts any remaining messages on the stack down one
C  position.  FTGMSG is called repeatedly until a blank message is
C  returned, which indicates that the stack is empty.  Each error message
C  may be up to 80 characters in length.  Another subroutine, called
C  FTCMSG, is available to simply clear the whole error message stack in
C  cases where one is not interested in the contents.
      call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
          print *,errmessage
          call ftgmsg(errmessage)
      end do
      end
c-----------------------------------------------------------------------
