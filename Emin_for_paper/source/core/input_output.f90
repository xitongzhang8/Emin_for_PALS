module input_output
  ! Contains all of the subroutines used to write and read to files, 
  ! both text and binary.
  !
  ! Input Values: 
  ! 10 : input coordinates file
  ! 20 : scratch for minima storage 
  ! 30 : sorted output 
  ! 50 : convergence data 
  use global_variables
  implicit none
  integer unique_minima
  character(len=40) :: file_name
  integer, parameter :: coords_in_io = 10
  integer, parameter :: coords_scratch_io = 20
  integer, parameter :: coords_out_io = 30
  integer, parameter :: convergence_io = 50

  integer :: ioerr
  contains
    subroutine open_io
      ! Initialise IO and read coordinates  
      !
      ! We have the coordinate read singlet code, 
      ! however, the default read is actually pretty intelligent.
      ! 
      open(unit=coords_in_io, file="coords", form="formatted", status="old", action="read")
      read(coords_in_io, *) X
      close(unit=coords_in_io)

      ! For saving working arrays
      ! Scratch: deleted on close.
      open(unit=coords_scratch_io, form="unformatted", status="scratch", action="readwrite", &
        iostat=ioerr, access='direct', recl=8*N )
      unique_minima = 0

      ! For printing convergence and other debug messages 
      open(unit=convergence_io, file="GMIN_OUT", form="formatted", status="replace", action="write")
    end subroutine open_io

    subroutine close_io()
      ! Close all IO and commit scratch file.
      close(coords_out_io)
      close(convergence_io)
      write(*,*) "Writing Database"

      call write_database()
    end subroutine close_io

    subroutine store_minima()
      ! If the minima is unique, write it to disk so that X may be reused.
      ! First, a check is performed for uniqueness
       if (is_unique(E)) then 
         unique_minima = unique_minima + 1
         e_store(unique_minima) = E
         write(coords_scratch_io, rec=unique_minima) X
       end if
    end subroutine store_minima

    subroutine write_database()
      ! Write the scratch array to file and sort the outputs.
      !
      integer, dimension(max_runs) :: sort_indicies
      integer :: i, indx
      ! This seems like a stupid way of doing things as they are nearly 
      ! identical, however, the formatted value needs a format, this breaks the 
      ! binary version...
      if (binary_io) then
        open(unit=coords_out_io, file='lowest', form='unformatted', &
          status='replace', action='write')
      else
        open(unit=coords_out_io, file='lowest', form='formatted', &
          status='replace', action='write')
      end if

      ! We pull the arrays out in the sorted order as given 
      ! by arg_sort
      sort_indicies = arg_sort(unique_minima)

      ! Read arrays and write to output file
      do i=1,unique_minima
        indx = sort_indicies(i)
        read(unit=coords_scratch_io, iostat=ioerr, rec=indx) X
        E = e_store(indx)
        if (binary_io) then
          ! Binary is far simpler
          write(unit=coords_out_io) X, E
        else
          ! Requires formatting for formatted output.
          call coordinate_write_formatted(i)
        end if
        if (i .eq. minima_save_limit) then
          ! Use this to only print off the first n minima
          exit
        end if
      end do

      close(coords_out_io)
      close(coords_scratch_io)
    end subroutine write_database

    subroutine coordinate_write_formatted(indx)
      ! Pretty output for the formatted coordinates
      !
      ! TODO: Maybe pass one of the indices?
      integer, intent(in) :: indx
      write(unit=coords_out_io,fmt='(A)') 
      write(unit=coords_out_io,fmt='(A, I3)') "Quench Number", indx
      write(unit=coords_out_io,fmt='(A, ES17.10, A, ES17.10, A, ES17.10)') "Energy is ", E
      write(unit=coords_out_io,fmt='(A, ES17.10)') "RMS is ", RMS
      write(unit=coords_out_io,fmt='(G16.10)') X 
      write(unit=coords_out_io,fmt='(A)') 
    end subroutine coordinate_write_formatted

    subroutine coordinate_read_singlet()
      ! Read the coordinates from a file.
      ! In this version, all coordinates are to be given on a separate line. 
      !
      ! Don't call this directly! It should be called by detect_file_format
      ! 
      integer :: i

      ! Read lines from coordinates file 
      open(unit=coords_in_io, file="coords", form="formatted", status="old", action="read")
      do i=1,N
        read(unit=coords_in_io, fmt="(F14.10)") X(i)
      end do
      close(unit=coords_in_io)
    end subroutine coordinate_read_singlet

    subroutine coordinate_read_triplet() 
      ! Read coordinates from the old 3-coordinate style.
      ! This is for backwards compatibility.
      !
    end subroutine coordinate_read_triplet 

    function arg_sort(unique_minima_l) result (indicies)
      ! Return the order of the indices that will sort the array.
      !
      ! Primarily designed to allow for the sorting of coordinate arrays in 
      ! the right order
      !
      ! Modified from
      ! https://github.com/certik/fortran-utils/blob/master/src/sorting.f90
      integer, intent(in) :: unique_minima_l
      integer :: indicies(unique_minima_l)
      double precision, dimension(unique_minima_l) :: unique_energies 

      integer :: i, imin, temp1
      double precision :: temp2, es2(unique_minima_l)

      unique_energies = e_store(:unique_minima_l)
      es2 = unique_energies

      ! Can do this implicitly
      do i=1,unique_minima
        indicies(i) = i
      end do

      do i=1, unique_minima-1
        ! Find the ith smallest energy
        imin = minloc(es2(i:),1) + i - 1
        ! Swap to position i in the temp energy and index array
        if (imin .ne. i) then 
          temp2 = es2(i)
          es2(i) = es2(imin)
          es2(imin) = temp2
          temp1     = indicies(i)
          indicies(i) = indicies(imin)
          indicies(imin) = temp1
        end if
      end do
    end function arg_sort

    logical function is_unique(e_test) 
      ! Test that the minima is not already on file.
      ! Currently uses a relative energy check, given by e_sim_max 
      !
      double precision, intent(in) :: e_test
      double precision :: relative_energy
      integer :: i

      ! If this is the first run then it is trivially unique.
      if (unique_minima .eq. 0) then
        is_unique = .true.
        if (debug) then
        write(*,*) "First mimima added"
        end if
        return
      end if
      ! Does a zero length loop do anything? 
      ! If not, we may be able to remove the if statement above.

      ! Else we check against all other energies
      do i=1,unique_minima
        relative_energy = abs((e_store(i) - e_test)/e_store(i))
        if (relative_energy < max_similarity) then
          is_unique = .false.
          ! Minima is already present
          if (debug) then
            write(*,'(A,G14.8,A,G14.8,A,G10.4)') "Minima already present  ", e_test, "   failed against ", &
              e_store(i), " Similarity ", relative_energy
          end if
          return

        end if
      end do

      if (debug) then
        is_unique = .true.
      end if
      write(*,*) "Minima is unique", e_test
    end function is_unique

    logical function is_new_coordinate_file()
      ! Detect if the file uses the new single coordinate format.
      ! Old uses triplets
      ! 
      ! TODO: This is a dummy at the moment
      !
      is_new_coordinate_file = .true.
    end function is_new_coordinate_file

    subroutine TEST_add_unique()
      ! Test the unique check
      E = 10d0 
      ! First allowed
      E = 40.5d0
      call store_minima()
      ! Rejected
      call store_minima()
      E = 41d0
      ! Allowed
      call store_minima()
      E = 10d0
      ! Rejected
      call store_minima()
      E = 40d0
      ! Allowed
      call store_minima()
      ! Rejected
      E = 40d0 + 40D-5
      call store_minima()
      
      write(*,*) "energy_store"
      write(*,'(*(G10.4))') e_store
      
      write(*,*) "unique energy"
      write(*,'(*(G10.4))') e_store(:unique_minima)

      ! Bonus, arg sort 
      write(*,*) arg_sort(unique_minima)
    end subroutine TEST_add_unique
end module


