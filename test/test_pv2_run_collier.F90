program run_loops
  use COLLIER

  Implicit none

  double precision :: p2, m12, m22, mu2
  double complex :: A01,A02,B0,B1,B00,C0
  double complex, allocatable :: Bcoeff(:,:),Bcoeffuv(:,:)
  double complex, allocatable :: Ccoeff(:,:,:),Ccoeffuv(:,:,:)
  double precision, allocatable :: Berr(:), Cerr(:)
  integer :: Brank, Crank, st
  character(len=255) :: read_line, file_name
  logical :: file_exists=.False.

  ! for B coefficients
  Brank = 2
  allocate(Bcoeff(0:Brank/2,0:Brank))
  allocate(Bcoeffuv(0:Brank/2,0:Brank))
  allocate(Berr(0:Brank))

  ! for C coefficients
  Crank = 3
  allocate(Ccoeff(0:Crank/2,0:Crank,0:Crank))
  allocate(Ccoeffuv(0:Crank/2,0:Crank,0:Crank))
  allocate(Cerr(0:Crank))

  ! open input file
  call get_command_argument(1, file_name)
  inquire(file=trim(file_name), exist=file_exists)
  if (.not.file_exists) then
     write(*,*) "Error: file ",trim(file_name)," does not exist!"
     stop 1
  end if
  open(100, file=trim(file_name), iostat=st)
  if (st.ne.0) then
     write(*,*) "Error: could not open file."
     stop 1
  end if

  call Init_cll(6,6,"")

  write (*,*) "# A0(m0^2) A0(m1^2) B0(p^2,m0^2,m1^2) B1(p^2,m0^2,m1^2) B00(p^2,m0^2,m1^2)"

  do
     read(100, "(a255)", iostat=st) read_line
     if (st.ne.0) then
        exit
     end if

     if (read_line(1:1).eq."#") cycle
     if (read_line.eq."") cycle

     backspace(100) ! jump to beginning of the line
     read(100,*) p2, m12, m22, mu2

     call SetMuUV2_cll(mu2)
     call SetMuIR2_cll(mu2)

     call A0_cll(A01,dcmplx(m12))
     call A0_cll(A02,dcmplx(m22))
     call B_cll(Bcoeff,Bcoeffuv,dcmplx(p2),dcmplx(m12),dcmplx(m22),Brank,Berr)
     call C_cll(Ccoeff,Ccoeffuv,dcmplx(0.0),dcmplx(0.0),dcmplx(0.0), &
                dcmplx(p2),dcmplx(m12),dcmplx(m22),Crank,Cerr)

     B0 = Bcoeff(0,0)
     B1 = Bcoeff(0,1)
     B00 = Bcoeff(1,0)
     C0 = Ccoeff(0,0,0)

     if (p2.eq.m12.and.p2.eq.0.0) C0 = 0
     if (p2.eq.m22.and.p2.eq.0.0) C0 = 0
     if (m12.eq.m22.and.m12.eq.0.0) C0 = 0

     write(*,*) real(A01), real(A02), real(B0), real(B1), real(B00), real(C0)
  end do

  close(100)

end program run_loops
