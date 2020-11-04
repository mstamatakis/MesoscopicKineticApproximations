module calculation_parser_module

  use constants_module, only: nnam0, maxwords2, lengthrecinp, csimfname, remchar, &
       iread, iwrite, cgenoutfname, inttermsmax
  use parser_module

  implicit none

  private
  type, public :: calculation_parser
     private

     real(8) temp, muIni, muFin, Dmu
     integer norigterms, ncorrterms, n_int, nshellsmax
     integer n_corr, nsites, nterms, n_eqns, n_terms, nsitesmax
     integer, dimension(inttermsmax, 2) :: int_term
     character(nnam0) approx

     character(nnam0) lattice_type, cluster_depth
     character(nnam0), allocatable :: corr_terms_names(:)
     character(nnam0), allocatable :: orig_terms_names(:)
     real(8), allocatable :: int_values(:)

   contains
     procedure, public :: read_setup => read_calculation_setup
     procedure, public :: cluster_setup => cluster_build

     ! Getters
     procedure, public :: get_temp
     procedure, public :: get_muIni
     procedure, public :: get_Dmu
     procedure, public :: get_muFin
     procedure, public :: get_approx
     procedure, public::  get_norigterms
     procedure, public::  get_ncorrterms
     procedure, public :: get_lattice_type
     procedure, public :: get_cluster_depth
     procedure, public :: get_orig_terms_names
     procedure, public :: get_corr_terms_names
     procedure, public :: get_int_values
     procedure, public :: get_nsites
     procedure, public :: get_nterms
     procedure, public :: get_int_term
     procedure, public :: get_n_int
     procedure, public :: get_n_corr
     procedure, public :: get_n_eqns
     procedure, public :: get_n_terms

   end type calculation_parser 

contains
  ! Getters
  function get_n_terms(this)
    implicit none
    class(calculation_parser), intent(in) :: this
    integer get_n_terms
    get_n_terms = this%n_terms
  end function get_n_terms
  
  function get_n_eqns(this)
    implicit none
    class(calculation_parser), intent(in) :: this
    integer get_n_eqns
    get_n_eqns = this%n_eqns
  end function get_n_eqns

  function get_n_corr(this)
    implicit none
    class(calculation_parser), intent(in) :: this
    integer get_n_corr
    get_n_corr = this%n_corr 
  end function get_n_corr

  function get_n_int(this)
    implicit none
    class(calculation_parser), intent(in) :: this
    integer get_n_int
    get_n_int = this%n_int
  end function get_n_int

  function get_int_term(this)
    implicit none
    class(calculation_parser), intent(in) :: this
    integer, dimension(inttermsmax,2) :: get_int_term
    integer i
    get_int_term = this%int_term
  end function get_int_term

  function get_nterms(this)
    implicit none
    class(calculation_parser), intent(in) :: this
    integer get_nterms
    get_nterms = this%nterms
  end function get_nterms

  function get_nsites(this)
    implicit none
    class(calculation_parser), intent(in) :: this
    integer get_nsites
    get_nsites = this%nsites
  end function get_nsites

  function get_temp(this)
    implicit none
    class(calculation_parser), intent(in) :: this
    real(8) get_temp
    get_temp = this%temp
  end function get_temp

  function get_muIni(this)
    implicit none
    class(calculation_parser), intent(in) :: this
    real(8) get_muIni
    get_muIni = this%muIni
  end function get_muIni

  function get_Dmu(this)
    implicit none
    class(calculation_parser), intent(in) :: this
    real(8) get_Dmu
    get_Dmu = this%dmu
  end function get_dmu

  function get_muFin(this)
    implicit none
    class(calculation_parser), intent(in) :: this
    real(8) get_muFin
    get_muFin = this%muFin
  end function get_muFin

  function get_approx(this)
    implicit none
    class(calculation_parser), intent(in) :: this
    character(nnam0) get_approx
    get_approx = this%approx
  end function get_approx

  function get_norigterms(this)
    implicit none
    class(calculation_parser), intent(in) :: this
    integer get_norigterms
    get_norigterms = this%norigterms
  end function get_norigterms

  function get_ncorrterms(this)
    implicit none
    class(calculation_parser), intent(in) :: this
    integer get_ncorrterms
    get_ncorrterms = this%ncorrterms
  end function get_ncorrterms

  function get_lattice_type(this)
    implicit none
    class(calculation_parser), intent(in) :: this
    character(nnam0) get_lattice_type
    get_lattice_type = this%lattice_type
  end function get_lattice_type

  function get_cluster_depth(this)
    implicit none
    class(calculation_parser), intent(in) :: this
    character(nnam0) get_cluster_depth
    get_cluster_depth = this%cluster_depth
  end function get_cluster_depth

  function get_corr_terms_names(this)
    implicit none
    class(calculation_parser), intent(in) :: this
    character(nnam0), allocatable :: get_corr_terms_names(:)
    allocate(get_corr_terms_names(0:this%ncorrterms))
    get_corr_terms_names = this%corr_terms_names
   end function get_corr_terms_names

  function get_orig_terms_names(this)
    implicit none
    class(calculation_parser), intent(in) :: this
    character(nnam0), allocatable :: get_orig_terms_names(:)
    allocate(get_orig_terms_names(0:this%norigterms))
     get_orig_terms_names = this%orig_terms_names
  end function get_orig_terms_names

  function get_int_values(this)
    implicit none
    class(calculation_parser), intent(in) :: this
    real(8), allocatable :: get_int_values(:)
    allocate(get_int_values(0:this%norigterms))
    get_int_values = this%int_values
  end function get_int_values

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine read_calculation_setup(this)
implicit none
class(calculation_parser), intent(inout) :: this

integer i, j, io, irec, nwords, nwords2, nwords3, p_nsites
integer iw(maxwords2), iw2(maxwords2), iw3(maxwords2), iw4(maxwords2)

character(lengthrecinp) recinput
character(lengthrecinp) tmpstr1, tmpstr2, tmpstr3, tmpstr4
character(1024) moreinfo

logical readtemp, readmu, readlatt, readclust, readnorig, readorig, readncorr, readcorr

! Parses the calculation setup file calculation_input.dat

! Level-1 keywords allowed:
!    lattice_type
!    cluster_depth
!    n_original_terms
!    original_terms
!    n_correction_terms
!    correction_terms
!    temperature
!    chempot_range

open(unit=iwrite,file=trim(cgenoutfname),status='unknown')

write(iwrite,'(/,a)') 'Calculation setup:'
write(iwrite,'(a)')   '~~~~~~~~~~~~~~~~~'

open(unit=iread,file=trim(csimfname),status='old')
readtemp = .false.
readlatt = .false.
readclust = .false.
readnorig = .false.
readorig = .false.
readncorr = .false.
readcorr = .false.
readmu = .false.
irec = 0
io = 0

call getrecord(iread,recinput,irec,io)

do while (io >= 0)

   call break_words(recinput,' ',remchar,nwords,iw)
   
   if (nwords > 0) then
      
       if (striccompare(recinput(iw(1):iw(2)),'lattice_type')) then
         
           if (readlatt) call error(-100,recinput(iw(1):iw(2))//';'//trim(int2str(irec))//';'//trim(csimfname))
          
           if (nwords /= 2) then
               moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
               call error(9,moreinfo)
           else               
               this%lattice_type = recinput(iw(3):iw(4))
           
               if (striccompare(trim(this%lattice_type),'hexagonal')) then                
                    write(iwrite,'(/,a)') '1. Lattice type= ' // trim(this%lattice_type)
               else                    
                    moreinfo = 'Unknown lattice directive ' // recinput(iw(3):iw(4)) // &
                     ' in line ' // trim(int2str(irec)) // '.'
                    call error(9,moreinfo)
               endif
               readlatt = .true.
            endif
            
       elseif (striccompare(recinput(iw(1):iw(2)),'cluster_depth')) then
           
           if (readclust) call error(-100,recinput(iw(1):iw(2))//';'//trim(int2str(irec))//';'//trim(csimfname))
          
           if (nwords /= 2) then
               moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
               call error(10,moreinfo)   
           else              
              this%cluster_depth = recinput(iw(3):iw(4))              
              if (this%cluster_depth.eq.'1NN') then
                  this%nsites = 7
                  this%approx = '1NN'
              elseif (this%cluster_depth.eq.'2NN') then
                  this%nsites = 13
                  this%approx = '2NN'
              elseif (this%cluster_depth.eq.'3NN') then
                  this%nsites = 19
                  this%approx = '3NN'
                    ! elseif  (this%cluster_depth.eq.'4NN') then
                    !     this%nsites = 31
                    ! elseif  (this%cluster_depth.eq.'5NN') then
                    !     this%nsites = 37
                    ! elseif  (this%cluster_depth.eq.'6NN') then
                    !     this%nsites = 43
                    ! elseif  (this%cluster_depth.eq.'7NN') then
                    !     this%nsites = 55
                    ! elseif  (this%cluster_depth.eq.'8NN') then
                    !     this%nsites = 61
                    ! elseif  (this%cluster_depth.eq.'9NN') then
                    !     this%nsites = 73
                    ! elseif  (this%cluster_depth.eq.'10NN') then
                    !     this%nsites = 85
                    ! elseif  (this%cluster_depth.eq.'11NN') then
                    !     this%nsites = 91
                    ! elseif  (this%cluster_depth.eq.'12NN') then
                    !    this%nsites = 97
               else
                  moreinfo = 'Unknown cluster directive ' // recinput(iw(3):iw(4)) // &
                  ' in line ' // trim(int2str(irec)) // '.'
                  call error(10,moreinfo)  
               endif       
                  write(iwrite,'(/,a)') '2. Cluster depth= ' // trim(this%approx)
                  readclust = .true.
                     this%nshellsmax = 3 !12 is the maximum number of shells in a cluster up to 12NN
                     this%nsitesmax = 19 !97 is the maxiumum number of sites in a cluster up to 12NN                  
           endif
                  
       elseif (striccompare(recinput(iw(1):iw(2)),'n_original_terms')) then
           
           if (readnorig) call error(-100,recinput(iw(1):iw(2))//';'//trim(int2str(irec))//';'//trim(csimfname))
           
           if (str_check_expression(recinput,nwords,iw,'AA/I4',.true.)) then
               this%norigterms = str2int(recinput(iw(3):iw(4)))            
               if ((this%norigterms < 0).or.(this%norigterms > 2)) then ! number of original terms < 0 or > 2 -> invalid!
                   moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                   call error(11,moreinfo)
               endif      
               readnorig = .true.
               write(iwrite,'(/,a)') '3. Number of original terms= ' // trim(int2str(this%norigterms))   
            endif

       elseif (striccompare(recinput(iw(1):iw(2)),'original_terms')) then
                
           if (readorig) call error(-100,recinput(iw(1):iw(2))//';'//trim(int2str(irec))//';'//trim(csimfname))

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Parsing:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  SB -1,172 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  1NN 0.3   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

            if (.not.readnorig) then
                moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                call error(12,moreinfo)
            else
                allocate(this%orig_terms_names(0:this%norigterms))
                allocate(this%int_values(0:this%norigterms))
            endif
           
            i=0
            do
              call getrecord(iread,recinput,irec,io)
              call break_words(recinput,' ',remchar,nwords,iw)
              if (nwords == 0) then
                  cycle
              elseif (str_check_expression(recinput,nwords,iw,'AA/R8',.true.)) then
                  i=i+1                     
                  if ((recinput(iw(1):iw(2)).ne.'SB').and.(recinput(iw(1):iw(2)).ne.'SBWI').and.&
                      (recinput(iw(1):iw(2)).ne.'1NN').and.(recinput(iw(1):iw(2)).ne.'2NN')) then
                      moreinfo = 'Unknown cluster directive ' // recinput(iw(1):iw(2)) // &
                      ' in line ' // trim(int2str(irec)) // '.'
                      call error(27,moreinfo)  
                  endif
                  this%orig_terms_names(i) = recinput(iw(1):iw(2))
                  this%int_values(i) = str2dbl(recinput(iw(3):iw(4)))                  
                  write(iwrite,'(/,a)') '('//recinput(iw(1):iw(2))//','//recinput(iw(3):iw(4))//')'
                  if (i == this%norigterms) then
                      exit
                  endif
              else
                  moreinfo = 'Unknown cluster directive ' // recinput(iw(3):iw(4)) // &
                  ' in line ' // trim(int2str(irec)) // '.'
                  call error(27,moreinfo)  
              endif
            enddo  
            readorig = .true.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
             
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Parsing:  (SB,-1.172) (1NN,0.3)!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
         !          if (.not.readnorig) then
         !              moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
         !              call error(12,moreinfo)
         !          endif
               
         !         if (nwords > this%norigterms + 1) then
         !             moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
         !             call error(13,moreinfo) 
         !         endif

         !         if (nwords < this%norigterms + 1) then
         !             moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
         !             call error(14,moreinfo) 
         !         endif           

         !         if (str_check_expression(recinput,nwords,iw,trim(int2str(this%norigterms+1)) // 'AA/',.true.)) then

         !            if (recinput(iw(3):iw(2*nwords)).ne.'(SB,-1.172) (1NN,0.3)') then
         !                 moreinfo = 'Unknown cluster directive ' // recinput(iw(3):iw(2*nwords)) // &
         !                 ' in line ' // trim(int2str(irec)) // '.'
         !                 call error(26,moreinfo)  
         !            endif 
                     
         !             tmpstr2 = recinput(iw(3):iw(2*nwords))
         !             call break_words(trim(tmpstr2),'(',remchar,nwords,iw2)
                   
         !             allocate(this%orig_terms_names(0:this%norigterms))
         !             allocate(this%int_values(0:this%norigterms))
         !             i=0
         !             do j = 1,2*this%norigterms,2
         !                 i=i+1
                      
         !                 tmpstr3 = tmpstr2(iw2(j):iw2(j+1))
         !                 call break_words(trim(tmpstr3),')',remchar,nwords,iw3)
                       
         !                 tmpstr4 = tmpstr3(iw3(1):iw3(2))
         !                 call break_words(trim(tmpstr4),',',remchar,nwords,iw4)
                       
         !                 this%orig_terms_names(i) = tmpstr4(iw4(1):iw4(2))
         !                 this%int_values(i) = str2dbl(tmpstr4(iw4(3):iw4(4)))

         !                 write(iwrite,'(/,a)') '('//tmpstr3(iw3(1):iw3(2))//')'
                      
         !             enddo
         !             readorig = .true.
                    
         !         endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             
       elseif (striccompare(recinput(iw(1):iw(2)),'n_corrections')) then
                 
           if (readncorr) call error(-100,recinput(iw(1):iw(2))//';'//trim(int2str(irec))//';'//trim(csimfname))

           if (str_check_expression(recinput,nwords,iw,'AA/I4',.true.)) then
               this%ncorrterms = str2int(recinput(iw(3):iw(4)))
                    
               if ((this%ncorrterms < 0).or.(this%ncorrterms > 3)) then ! number of correction terms < 0 or > 3 -> invalid!
                   moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                   call error(15,moreinfo)
               endif

               if ((this%approx.eq.'1NN').and.(this%ncorrterms > 2)) then ! number of correction terms  > 2 -> invalid!
                   moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                   call error(16,moreinfo)
               endif         
                    readncorr = .true.
                     write(iwrite,'(/,a)') '4. Number of correction terms= ' // trim(int2str(this%ncorrterms))
           endif

       elseif (striccompare(recinput(iw(1):iw(2)),'corrections')) then
                 
           if (readcorr) call error(-100,recinput(iw(1):iw(2))//';'//trim(int2str(irec))//';'//trim(csimfname))
                  
           if (.not.readncorr) then
               moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
               call error(17,moreinfo)
           endif
               
           if (nwords > this%ncorrterms + 1) then
               moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
               call error(18,moreinfo) 
           endif

           if (nwords < this%ncorrterms + 1) then
               moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
               call error(19,moreinfo) 
           endif 
                 
           if (str_check_expression(recinput,nwords,iw, trim(int2str(this%ncorrterms+1)) // 'AA/',.true.)) then
                      
               if (this%ncorrterms.eq.1) then
                   if ((recinput(iw(3):iw(4)).eq.'SB').or.(recinput(iw(3):iw(4)).eq.'SBWI')) then
                       allocate(this%corr_terms_names(0:this%ncorrterms))
                       this%corr_terms_names(1) =  recinput(iw(3):iw(4))
                       write(iwrite,'(/,a)') trim(this%corr_terms_names(1)) 
                       readcorr = .true.
                   else
                       moreinfo = 'Unknown cluster directive ' // recinput(iw(3):iw(4)) // &
                       ' in line ' // trim(int2str(irec)) // '.'
                       call error(23,moreinfo)  
                   endif
                endif
                        
                if (this%ncorrterms.eq.2) then
                    if ((recinput(iw(3):iw(4)).eq.'SB').and.(recinput(iw(5):iw(6)).eq.'1NN')) then
                        allocate(this%corr_terms_names(0:this%ncorrterms))
                        this%corr_terms_names(1) =  recinput(iw(3):iw(4))
                        write(iwrite,'(/,a)') trim(this%corr_terms_names(1)) 
                        this%corr_terms_names(2) = recinput(iw(5):iw(6))
                        write(iwrite,'(/,a)') trim(this%corr_terms_names(2)) 
                        readcorr = .true.
                     else
                        moreinfo = 'Unknown cluster directive ' // recinput(iw(3):iw(4))//' '// recinput(iw(5):iw(6))// &
                        ' in line ' // trim(int2str(irec)) // '.'
                        call error(24,moreinfo)
                     endif
                 endif                    
                    
                 if (this%ncorrterms.eq.3) then
                     if ((recinput(iw(3):iw(4)).eq.'SB').and.(recinput(iw(5):iw(6)).eq.'1NN').and.&
                         (recinput(iw(7):iw(8)).eq.'2NN')) then
                         allocate(this%corr_terms_names(0:this%ncorrterms))
                         this%corr_terms_names(1) =  recinput(iw(3):iw(4))
                         write(iwrite,'(/,a)') trim(this%corr_terms_names(1)) 
                         this%corr_terms_names(2) = recinput(iw(5):iw(6))
                         write(iwrite,'(/,a)') trim(this%corr_terms_names(2)) 
                         this%corr_terms_names(3) = recinput(iw(7):iw(8))
                         write(iwrite,'(/,a)') trim(this%corr_terms_names(3)) 
                         readcorr = .true.
                      else
                         moreinfo = 'Unknown cluster directive ' // recinput(iw(3):iw(4))//' ' // &
                         recinput(iw(5):iw(6))//' '// recinput(iw(7):iw(8))//& 
                         ' in line ' // trim(int2str(irec)) // '.'
                         call error(25,moreinfo)
                      endif                           
                   endif
           endif                

       elseif (striccompare(recinput(iw(1):iw(2)),'temperature')) then

           if (readtemp) call error(-100,recinput(iw(1):iw(2))//';'//trim(int2str(irec))//';'//trim(csimfname))
               
           if (str_check_expression(recinput,nwords,iw,'AA/R8',.true.)) then
               this%temp = str2dbl(recinput(iw(3):iw(4)))                      
               if (this%temp <= 0) then  ! temperature <= 0 -> invalid!
                   moreinfo = 'Occurred while parsing line ' // trim(int2str(irec)) // '.'
                   call error(20,moreinfo)                        
               endif
               readtemp = .true.        
               write(iwrite,'(/,a)') '5. Temperature= ' // trim(recinput(iw(3):iw(4)))    
            else
               moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
               call error(21,moreinfo) 
            endif
    
       elseif (striccompare(recinput(iw(1):iw(2)),'chempot_range')) then
                  
           if (readmu) call error(-100,recinput(iw(1):iw(2))//';'//trim(int2str(irec))//';'//trim(csimfname))
                  
           tmpstr1 = recinput(iw(3):iw(2*nwords))
           call break_words(trim(tmpstr1),':',remchar,nwords2,iw2)
           if (str_check_expression(tmpstr1,nwords2,iw2,'3R8',.true.)) then
                 
               this%muIni = str2dbl(tmpstr1(iw2(1):iw2(2)))
               this%Dmu = str2dbl(tmpstr1(iw2(3):iw2(4)))
               this%muFin = str2dbl(tmpstr1(iw2(5):iw2(6)))
               readmu = .true.
                     
               write(iwrite,'(/,a)') '6. Chemical potential range=' // trim(tmpstr1(iw2(1):iw2(2)))// ':' &
               // trim(tmpstr1(iw2(3):iw2(4))) // ':' // trim(tmpstr1(iw2(5):iw2(6)))
                     
           else         
               moreinfo = 'Invalid expression "' // recinput(iw(1):iw(2*nwords)) // '" in line ' // trim(int2str(irec)) // '.'
               call error(22,moreinfo)                    
            endif
            
       endif
              
   endif

   call getrecord(iread,recinput,irec,io)

enddo

close(iread)

    ! Check for validity of input
    if (.not.readlatt) then ! lattice type not specified
        call error(1)
    endif 

    if (.not.readclust) then ! cluster depth  not specified
        call error(2)
     endif
     
    if (.not.readnorig) then ! number of original terms not specified
        call error(3)
    endif

    if (.not.readorig) then ! original term names not specified
        call error(4)
    endif

    if (.not.readncorr) then ! number of corrections not specified
        call error(5)
    endif

    if (.not.readcorr) then ! correction names not specified
        call error(6)
    endif

    if (.not.readtemp) then ! temperature not specified
        call error(7)
    endif

    if (.not.readmu) then ! raneg of the chemical potential not specified
        call error(8) 
     endif
     
return

end subroutine read_calculation_setup

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine cluster_build(this)
implicit none
class(calculation_parser), intent(inout) :: this
integer n1, n2, n3
integer, dimension(this%nsitesmax,this%nshellsmax) :: coord
character(nnam0) cluster_depth
real(8), dimension(this%nsitesmax) :: x
real(8), dimension(this%nsitesmax) :: y

call  cluster(this,n1,n2,n3,x,y,coord)

if (this%cluster_depth.eq.'1NN') then
   
    if (this%ncorrterms.eq.1) then
      
        if (this%corr_terms_names(1).eq.'SB') then      
            this%n_corr = 1
            this%n_eqns = 1
            this%n_terms = 2
            call int_01(this,n1,coord)
         else if (this%corr_terms_names(1).eq.'SBWI') then       
            this%n_corr = 1
            this%n_eqns = 1
            this%n_terms = 2
            call int_01(this,n1,coord)
            call int_11_1nn(this,n1,x,y,coord)    
         endif
         
    elseif (this%ncorrterms.eq.2) then

        if ((this%corr_terms_names(1).eq.'SB').and.(this%corr_terms_names(2).eq.'1NN')) then              
            this%n_corr = 2
            this%n_eqns = 2
            this%n_terms = 4
            call int_01(this,n1,coord)
            call int_11_1nn(this,n1,x,y,coord)   
        endif
            
    endif

endif

if (this%cluster_depth.eq.'2NN') then

    if (this%ncorrterms.eq.1) then

        if (this%corr_terms_names(1).eq.'SB') then      
            this%n_corr = 2
            this%n_eqns = 2
            this%n_terms = 3
            call int_01(this,n1,coord)                
        elseif (this%corr_terms_names(1).eq.'SBWI') then                
            this%n_corr = 2
            this%n_eqns = 2
            this%n_terms = 3              
            call int_01(this,n1,coord)
            call int_11_1nn(this,n1,x,y,coord)
            call int_12_1nn(this,n1,n2,x,y,coord)
            call int_22_1nn(this,n2,x,y,coord) !Not present    
        endif

     elseif (this%ncorrterms.eq.2) then

         if ((this%corr_terms_names(1).eq.'SB').and.(this%corr_terms_names(2).eq.'1NN')) then
             this%n_corr = 4
             this%n_eqns = 4
             this%n_terms = 6
             call int_01(this,n1,coord)
             call int_11_1nn(this,n1,x,y,coord)                       
             call int_12_1nn(this,n1,n2,x,y,coord)
             call int_22_1nn(this,n2,x,y,coord) !Not present    
         endif

     elseif (this%ncorrterms.eq.3) then
            
         if ((this%corr_terms_names(1).eq.'SB').and.(this%corr_terms_names(2).eq.'1NN').and. &
             (this%corr_terms_names(3).eq.'2NN')) then               
             this%n_corr = 6
             this%n_eqns = 6
             this%n_terms = 9
             call int_01(this,n1,coord)                  
             call int_11_1nn(this,n1,x,y,coord)
             call int_12_1nn(this,n1,n2,x,y,coord)
             call int_22_1nn(this,n2,x,y,coord) !Not present
             ! Encode all the 2NN interaction terms  
             call int_11_2nn(this,n1,x,y,coord)
             call int_12_2nn(this,n1,n2,x,y,coord) !Not present
             call int_22_2nn(this,n2,x,y,coord)
          endif
          
     endif

endif


if (this%cluster_depth.eq.'3NN') then

    if (this%ncorrterms.eq.1) then

        if (this%corr_terms_names(1).eq.'SB') then
            this%n_corr = 3
            this%n_eqns = 3
            this%n_terms = 4
            call int_01(this,n1,coord)            
        elseif (this%corr_terms_names(1).eq.'SBWI') then
            this%n_corr = 3
            this%n_eqns = 3
            this%n_terms = 4         
            call int_01(this,n1,coord)
            call int_11_1nn(this,n1,x,y,coord)
            call int_12_1nn(this,n1,n2,x,y,coord)
            call int_22_1nn(this,n2,x,y,coord) !Not needed
            call int_13_1nn(this,n2,n3,x,y,coord) 
            call int_23_1nn(this,n2,n3,x,y,coord)
            call int_33_1nn(this,n3,x,y,coord) !Not present
         endif
               
     elseif (this%ncorrterms.eq.2) then
             
         if ((this%corr_terms_names(1).eq.'SB').and.(this%corr_terms_names(2).eq.'1NN')) then                    
             this%n_corr = 7
             this%n_eqns = 7
             this%n_terms = 9
             call int_01(this,n1,coord)
             call int_11_1nn(this,n1,x,y,coord)
             call int_12_1nn(this,n1,n2,x,y,coord)
             call int_22_1nn(this,n2,x,y,coord) !Not needed
             call int_13_1nn(this,n2,n3,x,y,coord) 
             call int_23_1nn(this,n2,n3,x,y,coord)
             call int_33_1nn(this,n3,x,y,coord) !Not present
         endif
              
      elseif (this%ncorrterms.eq.3) then
         
          if ((this%corr_terms_names(1).eq.'SB').and.(this%corr_terms_names(2).eq.'1NN').and. &
              (this%corr_terms_names(3).eq.'2NN')) then     
              this%n_corr = 10
              this%n_eqns = 10
              this%n_terms = 13
              call int_01(this,n1,coord)
              call int_11_1nn(this,n1,x,y,coord)
              call int_12_1nn(this,n1,n2,x,y,coord)
              call int_22_1nn(this,n2,x,y,coord) !Not needed                  
              call int_13_1nn(this,n2,n3,x,y,coord) 
              call int_23_1nn(this,n2,n3,x,y,coord)
              call int_33_1nn(this,n3,x,y,coord) !Not present             
              ! Encode all the 2NN interaction terms                  
              call int_11_2nn(this,n1,x,y,coord)
              call int_12_2nn(this,n1,n2,x,y,coord) !Not present
              call int_22_2nn(this,n2,x,y,coord)                  
              call int_13_2nn(this,n1,n3,x,y,coord)
              call int_23_2nn(this,n2,n3,x,y,coord) !Not present
              call int_33_2nn(this,n3,x,y,coord) !Not present     
           endif
               
      endif
        
endif

return
end subroutine cluster_build

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
subroutine int_01(this,n1,coord)
    implicit none
    class(calculation_parser), intent(inout) :: this
    integer n1, j
    integer, dimension(this%nsitesmax,this%nshellsmax) :: coord
   
    do j = 1,n1
        this%int_term(j,1:2) = (/1,coord(j,1)/)
    enddo
    this%n_int = n1              
    this%nterms = this%nsites + this%n_int  

    return          
 end subroutine int_01

 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
subroutine int_11_1nn(this,n1,x,y,coord)
    implicit none
    class(calculation_parser), intent(inout) :: this
    integer n1, i, j
    integer, dimension(this%nsitesmax,this%nshellsmax) :: coord
    real(8), dimension(this%nsitesmax) :: x
    real(8), dimension(this%nsitesmax) :: y
    real(8) dv
   
    do i = 1,n1-1                
        do j = i+1,n1      
          dv = sqrt(((x(coord(i,1))-x(coord(j,1)))**2)+ &
               ((y(coord(i,1))-y(coord(j,1)))**2))                 
          if (dv.le.1.0d0) then          
          this%n_int = this%n_int+1
          this%int_term(this%n_int,1:2) = (/coord(i,1),coord(j,1)/)  
          endif                    
       enddo                
     enddo
     this%nterms = this%nsites+this%n_int
     
     return          
end subroutine int_11_1nn

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
 
subroutine int_12_1nn(this,n1,n2,x,y,coord)
    implicit none
    class(calculation_parser), intent(inout) :: this
    integer n1, n2, i, j
    integer, dimension(this%nsitesmax,this%nshellsmax) :: coord
    real(8), dimension(this%nsitesmax) :: x
    real(8), dimension(this%nsitesmax) :: y
    real(8) dv
    
     do i = 1,n1                 
         do j = 1,n2                    
             dv = sqrt(((x(coord(i,1))-x(coord(j,2)))**2)+ &
                  ((y(coord(i,1))-y(coord(j,2)))**2))                  
             if (dv.le.1.0d0) then          
                 this%n_int = this%n_int+1
                 this%int_term(this%n_int,1:2) = (/coord(i,1),coord(j,2)/)                      
             endif                   
         enddo          
      enddo
      this%nterms = this%nsites+this%n_int

    return          
end subroutine int_12_1nn

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
subroutine int_13_1nn(this,n1,n3,x,y,coord)
    implicit none
    class(calculation_parser), intent(inout) :: this
    integer n1, n3, i, j
    integer, dimension(this%nsitesmax,this%nshellsmax) :: coord
    real(8), dimension(this%nsitesmax) :: x
    real(8), dimension(this%nsitesmax) :: y
    real(8) dv
    
    do i = 1,n1                  
        do j = 1,n3                    
            dv=sqrt(((x(coord(i,1))-x(coord(j,3)))**2)+ &
               ((y(coord(i,1))-y(coord(j,3)))**2))                  
            if (dv.le.1.0d0) then          
                this%n_int = this%n_int+1
                this%int_term(this%n_int,1:2) = (/coord(i,1),coord(j,3)/)                     
         endif                  
       enddo                 
    enddo              
    this%nterms = this%nsites+this%n_int

    return          
end subroutine int_13_1nn

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
subroutine int_22_1nn(this,n2,x,y,coord)
    implicit none
    class(calculation_parser), intent(inout) :: this
    integer n2, i, j
    integer, dimension(this%nsitesmax,this%nshellsmax) :: coord
    real(8), dimension(this%nsitesmax) :: x
    real(8), dimension(this%nsitesmax) :: y
    real(8) dv
    
    do i = 1,n2-1                 
        do j = i+1,n2                    
            dv=sqrt(((x(coord(i,2))-x(coord(j,2)))**2)+ &
               ((y(coord(i,2))-y(coord(j,2)))**2))                   
            if (dv.le.1.0d0) then         
               this%n_int = this%n_int+1
               this%int_term(this%n_int,1:2) = (/coord(i,2),coord(j,2)/)                      
            endif                    
       enddo                
    enddo
    this%nterms = this%nsites+this%n_int

    return          
end subroutine int_22_1nn

 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine int_23_1nn(this,n2,n3,x,y,coord)
    implicit none
    class(calculation_parser), intent(inout) :: this
    integer n2, n3, i, j
    integer, dimension(this%nsitesmax,this%nshellsmax) :: coord
    real(8), dimension(this%nsitesmax) :: x
    real(8), dimension(this%nsitesmax) :: y
    real(8) dv
    
    do i = 1,n2                 
        do j = 1,n3                     
            dv=sqrt(((x(coord(i,2))-x(coord(j,3)))**2)+ &
               ((y(coord(i,2))-y(coord(j,3)))**2))   
            if (dv.le.1.0d0) then         
                this%n_int = this%n_int+1
                this%int_term(this%n_int,1:2) = (/coord(i,2),coord(j,3)/)                      
            endif                    
        enddo                 
    enddo     
    this%nterms = this%nsites+this%n_int
   
    return
end subroutine int_23_1nn

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
subroutine int_33_1nn(this,n3,x,y,coord)
    implicit none
    class(calculation_parser), intent(inout) :: this
    integer n3, i, j
    integer, dimension(this%nsitesmax,this%nshellsmax) :: coord
    real(8), dimension(this%nsitesmax) :: x
    real(8), dimension(this%nsitesmax) :: y
    real(8) dv
    
    do i = 1,n3-1                 
        do j = i+1,n3                     
          dv=sqrt(((x(coord(i,3))-x(coord(j,3)))**2)+ &
             ((y(coord(i,3))-y(coord(j,3)))**2))                  
          if (dv.le.1.0d0) then          
              this%n_int = this%n_int+1
              this%int_term(this%n_int,1:2) = (/coord(i,3),coord(j,3)/)                     
          endif                    
       enddo              
    enddo
    this%nterms = this%nsites+this%n_int

    return          
end subroutine int_33_1nn

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine int_11_2nn(this,n1,x,y,coord)
    implicit none
    class(calculation_parser), intent(inout) :: this
    integer n1, i, j
    integer, dimension(this%nsitesmax,this%nshellsmax) :: coord
    real(8), dimension(this%nsitesmax) :: x
    real(8), dimension(this%nsitesmax) :: y
    real(8) dv
    
    do i = 1,n1-1                 
        do j = i+1,n1                     
          dv=sqrt(((x(coord(i,1))-x(coord(j,1)))**2)+ &
             ((y(coord(i,1))-y(coord(j,1)))**2))                
          if ((dv.gt.1.0d0).and.(dv.le.sqrt(3.0d0))) then
              this%n_int = this%n_int+1
              this%int_term(this%n_int,1:2) = (/coord(i,1),coord(j,1)/)                   
           endif                   
        enddo                 
    enddo

    this%nterms = this%nsites+this%n_int

    return          
end subroutine int_11_2nn

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
subroutine int_12_2nn(this,n1,n2,x,y,coord)
    implicit none
    class(calculation_parser), intent(inout) :: this
    integer n1, n2, i, j
    integer, dimension(this%nsitesmax,this%nshellsmax) :: coord
    real(8), dimension(this%nsitesmax) :: x
    real(8), dimension(this%nsitesmax) :: y
    real(8) dv
    
     do i = 1,n1                  
         do j = 1,n2                    
             dv=sqrt(((x(coord(i,1))-x(coord(j,2)))**2)+ &
               ((y(coord(i,1))-y(coord(j,2)))**2))                  
             if ((dv.gt.1.0d0).and.(dv.le.sqrt(3.0d0))) then         
                this%n_int = this%n_int+1
                this%int_term(this%n_int,1:2) = (/coord(i,1),coord(j,2)/)                   
           endif                    
         enddo                 
     enddo
     this%nterms = this%nsites+this%n_int

    return          
end subroutine int_12_2nn

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine int_22_2nn(this,n2,x,y,coord)
    implicit none
    class(calculation_parser), intent(inout) :: this
    integer n2, i, j
    integer, dimension(this%nsitesmax,this%nshellsmax) :: coord
    real(8), dimension(this%nsitesmax) :: x
    real(8), dimension(this%nsitesmax) :: y
    real(8) dv
    
    do i = 1,n2-1                 
        do j = i+1,n2                     
            dv=sqrt(((x(coord(i,2))-x(coord(j,2)))**2)+ &
               ((y(coord(i,2))-y(coord(j,2)))**2))
            if ((dv.gt.1.0d0).and.(dv.le.sqrt(3.0d0))) then                          
                this%n_int = this%n_int+1
                this%int_term(this%n_int,1:2) = (/coord(i,2),coord(j,2)/)                      
            endif                         
         enddo                 
    enddo
    this%nterms = this%nsites+this%n_int

    return          
end subroutine int_22_2nn

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
subroutine int_13_2nn(this,n1,n3,x,y,coord)
    implicit none
    class(calculation_parser), intent(inout) :: this
    integer n1, n3, i, j
    integer, dimension(this%nsitesmax,this%nshellsmax) :: coord
    real(8), dimension(this%nsitesmax) :: x
    real(8), dimension(this%nsitesmax) :: y
    real(8) dv
  
    do i = 1,n1                
        do j = 1,n3                    
            dv=sqrt(((x(coord(i,1))-x(coord(j,3)))**2)+ &
               ((y(coord(i,1))-y(coord(j,3)))**2))                   
            if ((dv.gt.1.0d0).and.(dv.le.sqrt(3.0d0))) then                
                this%n_int = this%n_int+1
                this%int_term(this%n_int,1:2) = (/coord(i,1),coord(j,3)/)                                               
            endif                    
         enddo                
    enddo
    this%nterms = this%nsites+this%n_int

   return          
end subroutine int_13_2nn

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine int_23_2nn(this,n2,n3,x,y,coord)
    implicit none
    class(calculation_parser), intent(inout) :: this
    integer n2, n3, i, j
    integer, dimension(this%nsitesmax,this%nshellsmax) :: coord
    real(8), dimension(this%nsitesmax) :: x
    real(8), dimension(this%nsitesmax) :: y
    real(8) dv
  
    do i = 1,n2                  
        do j = 1,n3        
          dv=sqrt(((x(coord(i,2))-x(coord(j,3)))**2)+ &
             ((y(coord(i,2))-y(coord(j,3)))**2))                  
          if ((dv.gt.1.0d0).and.(dv.le.sqrt(3.0d0))) then                
              this%n_int = this%n_int+1
              this%int_term(this%n_int,1:2) = (/coord(i,2),coord(j,3)/)                      
          endif                   
        enddo                 
    enddo
   this%nterms = this%nsites+this%n_int

   return          
end subroutine int_23_2nn 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine int_33_2nn(this,n3,x,y,coord)
    implicit none
    class(calculation_parser), intent(inout) :: this
    integer n3, i, j
    integer, dimension(this%nsitesmax,this%nshellsmax) :: coord
    real(8), dimension(this%nsitesmax) :: x
    real(8), dimension(this%nsitesmax) :: y
    real(8) dv
    
     do i = 1,n3-1                 
         do j = i+1,n3                     
             dv=sqrt(((x(coord(i,3))-x(coord(j,3)))**2)+ &
               ((y(coord(i,3))-y(coord(j,3)))**2))                  
             if ((dv.gt.1.0d0).and.(dv.le.sqrt(3.0d0))) then                          
                 this%n_int = this%n_int+1
                 this%int_term(this%n_int,1:2) = (/coord(i,3),coord(j,3)/)                         
             endif                    
         enddo                
     enddo
     this%nterms = this%nsites+this%n_int

    return          
end subroutine int_33_2nn

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
subroutine cluster(this,n1,n2,n3,x,y,coord)
    implicit none
    class(calculation_parser), intent(inout) :: this
    integer, intent(out) :: n1,n2,n3
    integer, dimension(this%nsitesmax,this%nshellsmax) :: coord
    real(8), dimension(this%nsitesmax) :: x
    real(8), dimension(this%nsitesmax) :: y
    integer i
    real(8) dv

    !!!!!!!!!! Center                
    x(1) = 0.0                   !******* 0NN site (1)***************
    y(1) = 0.0
    !!!!!!!!!!! first shell        
    x(2) = 1.0                   !******* 1NN site (2)****************
    y(2) = 0.0
    
    x(3) = 1.0/2.0               !******* 1NN site (3)*****************
    y(3) = -sqrt(3.0)/2.0

    x(4) = -1.0/2.0              !******* 1NN site (4)*****************
    y(4) = -sqrt(3.0)/2.0

    x(5) = -1.0                  !******* 1NN site (5)*****************
    y(5) = 0.0

    x(6) = -1.0/2.0              !******* 1NN site (6)*****************
    y(6) = sqrt(3.0)/2.0 
   
    x(7) = 1.0/2.0               !******* 1NN site (7)*****************
    y(7) = sqrt(3.0)/2.0
    !!!!!!!!!!!! second shell       
    x(8) = 3.0/2.0               !******* 2NN site (8)****************
    y(8) = -sqrt(3.0)/2.0

    x(9) = 0.0                   !******* 2NN site (9)*****************
    y(9) = -sqrt(3.0)

    x(10) = -3.0/2.0             !******* 2NN site (10)*****************
    y(10) = -sqrt(3.0)/2.0

    x(11) = -3.0/2.0             !******* 2NN site (11)*****************
    y(11) = sqrt(3.0)/2.0

    x(12) = 0.0                  !******* 2NN site (12)*****************
    y(12) = sqrt(3.0) 
   
    x(13) = 3.0/2.0              !******* 2NN site (13)*****************
    y(13) = sqrt(3.0)/2.0
    !!!!!!!!!!!! third shell       
    x(14) = 2.0                   !******* 3NN site (14)****************
    y(14) = 0.0

    x(15) = 1.0                   !******* 3NN site (15)*****************
    y(15) = -sqrt(3.0)

    x(16) = -1.0                 !******* 3NN site (16)*****************
    y(16) = -sqrt(3.0)

    x(17) = -2.0                 !******* 3NN site (17)*****************
    y(17) = 0.0

    x(18) = -1.0                 !******* 3NN site (18)*****************
    y(18) = sqrt(3.0) 
   
    x(19) = 1.0                  !******* 3NN site (19)*****************
    y(19) = sqrt(3.0)

    n1 = 0
    n2 = 0
    n3 = 0
    do i = 2,this%nsitesmax      
        dv = sqrt(((x(1)-x(i))**2)+((y(1)-y(i))**2))       
        if (dv.le.1.0d0) then
           n1 = n1+1
           coord(n1,1) = i
        endif   
        if (dv.gt.1.0d0) then
            if (dv.le.sqrt(3.0d0)) then
                n2 = n2+1
                coord(n2,2) = i
            endif
        endif          
        if (dv.gt.sqrt(3.0d0)) then
            if (dv.le.2.0d0) then
               n3 = n3+1
               coord(n3,3) = i
            endif
        endif     
    enddo
    
    return
end subroutine cluster

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine error(ierror,moreinfo,errorout)

use err_warn_module
use constants_module, only : iwrite

implicit none

integer, parameter           :: errCallerID = 2000
character(len=callerNameLen) :: errCallerName = 'calculation_parser_module'

integer :: ierror, ierror_local
character(*), optional, intent(in) :: moreinfo
integer, optional, intent(in) :: errorout

character(512) :: errormessage

ierror_local = ierror ! we need a local variable, ierror_local, since we may need to override it later

select case (ierror_local)

     case (-100)
        errormessage = '' ! error for invalid duplicate specification
        
     case (1)
        errormessage = 'lattice type has not been specified in ' // trim(csimfname) // '.'

     case (2)
        errormessage = 'cluster depth has not been specified in ' // trim(csimfname) // '.'

     case (3)
        errormessage = 'n_original_terms have not been specified in ' // trim(csimfname) // '.'

     case (4)
        errormessage = 'original_terms have not been specified in ' // trim(csimfname) // '.'

     case (5)
        errormessage = 'n_corrections have not been specified in ' // trim(csimfname) // '.'

      case (6)
         errormessage = 'corrections have not been specified in ' // trim(csimfname) // '.'

      case (7)
         errormessage = 'temperature has not been specified in ' // trim(csimfname) // '.'
        
      case (8)
         errormessage = 'chempot_range has not been specified in ' // trim(csimfname) // '.'
         
      case (9)
         errormessage = ' must start with the keyword lattice_type followed by ' // &
                       'the following directive: hexagonal.'       
      case (10)
         errormessage =' must start with the keyword cluster depth followed by ' // &
                       'one of the following directives: 1NN, 2NN, or 3NN.'
        
      case (11)
         errormessage = 'invalid number of original terms has been specified in ' // trim(csimfname) // '.' 

       case(12)
          errormessage = 'original_terms keyword must be preceded by n_original_terms keyword in ' // trim(csimfname) // '.'

       case(13)   
          errormessage = 'more than n_original_term names appear in ' // trim(csimfname) // '.'
          
       case(14)   
          errormessage = 'less than n_original_term names appear in ' // trim(csimfname) // '.'

       case (15)
          errormessage = 'invalid number of correction terms has been specified in ' // trim(csimfname) // '.'
          
       case (16)
          errormessage =' cluster depth is up to 1NN and the n_corrections must not exceed two in '// trim(csimfname) // '.'
          
       case(17)
          errormessage = 'corrections keyword must be preceded by n_corrections keyword in ' // trim(csimfname) // '.'

       case(18)   
          errormessage = 'more than n_correction names appear in ' // trim(csimfname) // '.'

       case(19)   
          errormessage = 'less than n_correction names appear in ' // trim(csimfname) // '.'

       case (20)
          errormessage = 'invalid non-positive temperature has been specified in ' // trim(csimfname) // '.'
   
       case (21)
          errormessage =  'temperature keyword must be followed by a single real in ' // trim(csimfname) // '.'

       case (22)
          errormessage = 'chemical potential keyword must be followed by three reals in ' // trim(csimfname) // '.'

       case (23)
          errormessage = ' must start with the keyword corrections followed by ' // &
                        'the following directive: SB or SBWI.'
       case (24)
          errormessage = ' must start with the keyword corrections followed by ' // &
                        'the following directive: SB 1NN.'
       case (25)
          errormessage = ' must start with the keyword corrections followed by ' // &
                        'the following directive: SB 1NN 2NN.'
       !case (26)
       !   errormessage = ' must start with the keyword original_terms  followed by ' // &
       !                 'the following directive: (SB,-1.172) (1NN,0.3) .'

       case (27)
          errormessage =' must start with the keyword original terms followed by ' // &
                        'allowed directives.'   

end select

call errwarn_globj%error(errCallerID, errCallerName, ierror_local, errorout, errormessage, moreinfo)

return

end subroutine error

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


end module calculation_parser_module
