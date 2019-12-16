module polymer_reader

implicit none

contains

subroutine read_orbitals(c) bind (C, name='read_orbitals_')

    use iso_c_binding

    implicit none

    double precision, intent(out) :: c(*)
    integer :: numorb, maxngrid, natom
    integer :: i, j, io
    logical :: e

    open(7, file='ORBITALS.DAT.BIN', status='old', form='unformatted')
    inquire(7, exist=e)

    if (.not. e) then
        write(*,*) 'Orbital coefficients file was not found! Exiting..'
        stop
    end if 

    read(7) numorb, maxngrid, natom
    write(*, '(a, 1x, 3(i6, 1x))') 'Metadata from orbital file', numorb, maxngrid, natom

    do i = 1, numorb * maxngrid
      read(7) c(i)
    end do

    close(7)

end subroutine read_orbitals

subroutine read_integrals(eri_index, eri) bind (C, name='read_integrals_')

    use iso_c_binding

    implicit none

    integer, intent(out) :: eri_index(*)
    double precision, intent(out) :: eri(*)
    double precision :: e_nuc
    integer :: numorb, numpair, numeri
    integer :: i, idx(4)
    logical :: e

    open(7, file='FCIDUMP.BIN', status='old', form='unformatted')
    inquire(7, exist=e)

    if (.not. e) then
        write(*,*) 'FCIDUMP was not found! Exiting..'
        stop
    end if 

    read(7) numorb
    write(*, '(a, 1x, i6)') 'Number of orbitals as extracted from FCIDUMP', numorb

    ! Number of pairs 
    numpair = numorb * (numorb + 1) / 2
    write(*,'(a, i6)') 'Number of Hcore matrix elements', numpair

    ! Number of ERI-s
    numeri = numpair * (numpair + 1) / 2
    write(*,'(a, i6)') 'Number of ERI-s', numeri

    read(7) e_nuc, idx(4)
    write(*, '(a, f20.10, 4(i6,1x))') 'Nuclear repulsion energy = ', e_nuc

    ! Read core Hamiltonian
    do i = 1, numpair
      read(7) eri(i), idx
      eri_index(4*(i - 1) + 1 : 4*(i - 1) + 4) = idx(1:4)
    end do
    write(*,*) 'Finished reading core Hamiltonian'

    ! Read electron repulsion integrals
    do i = 1, numeri
      read(7) eri(numpair + i), idx(4)
      eri_index(4*(i - 1) + 4*numpair + 1 : 4*(i - 1) + 4*numpair + 4) = idx(1:4)
    end do 

    write(*,*) 'Finished reading electron replsion integrals'

end subroutine read_integrals

end module polymer_reader
