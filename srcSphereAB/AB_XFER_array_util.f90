!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!-*- F90 -*- so emacs thinks this is an f90 file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Name: AB_XFER_array_util.F90
!
! Purpose: Holds various handy AB utilities provided for the users
!          convienience. These utilities deal with arrays.
!
! History:
! 2/22/01 Robert Oehmke - created
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module AB_XFER_array_util
  use AB_module
  implicit none

contains
  subroutine AB_array3_gc_add_size(long, lat, alt, gcn, size)
    integer, intent(inout) :: size(ab_num_nbrs)
    integer, intent(in) :: long, lat, alt, gcn

    size(ab_north) = size(ab_north) + long*gcn*alt
    size(ab_south) = size(ab_south) + long*gcn*alt
    size(ab_east) = size(ab_east) + lat*gcn*alt
    size(ab_west) = size(ab_west) + lat*gcn*alt
    size(ab_northeast) = size(ab_northeast) + gcn*gcn*alt
    size(ab_northwest) = size(ab_northwest) + gcn*gcn*alt
    size(ab_southeast) = size(ab_southeast) + gcn*gcn*alt
    size(ab_southwest) = size(ab_southwest) + gcn*gcn*alt
  end subroutine AB_array3_gc_add_size

  subroutine AB_array3_gc_pack(long, lat, alt, gcn, v_in, dir, pole, p, out_array)
    use AB_module
    integer, intent(in) :: long, lat, alt, gcn
    integer, intent(in) :: dir
    logical, intent(in) :: pole
    integer, intent(inout) :: p
    real, dimension(1 - gcn:, 1 - gcn:, 1:), intent(in) :: v_in
    real, dimension(:), intent(inout) :: out_array
    integer :: i, j, a

    select case (dir)

    case (ab_north)
      if (pole) then
        do a = 1, alt
          do i = lat, lat - gcn + 1, -1
            do j = 1, long
              out_array(p) = v_in(j, i, a)
              p = p + 1
            end do
          end do
        end do
      else
        do a = 1, alt
          do i = lat - gcn + 1, lat
            do j = 1, long
              out_array(p) = v_in(j, i, a)
              p = p + 1
            end do
          end do
        end do
      end if

    case (ab_south)
      if (pole) then
        do a = 1, alt
          do i = gcn, 1, -1
            do j = 1, long
              out_array(p) = v_in(j, i, a)
              p = p + 1
            end do
          end do
        end do
      else
        do a = 1, alt
          do i = 1, gcn
            do j = 1, long
              out_array(p) = v_in(j, i, a)
              p = p + 1
            end do
          end do
        end do
      end if

    case (ab_east)
      do a = 1, alt
        do i = 1, lat
          do j = long - gcn + 1, long
            out_array(p) = v_in(j, i, a)
            p = p + 1
          end do
        end do
      end do

    case (ab_west)
      do a = 1, alt
        do i = 1, lat
          do j = 1, gcn
            out_array(p) = v_in(j, i, a)
            p = p + 1
          end do
        end do
      end do

    case (ab_northeast)
      if (pole) then
        do a = 1, alt
          do i = lat, lat - gcn + 1, -1
            do j = long - gcn + 1, long
              out_array(p) = v_in(j, i, a)
              p = p + 1
            end do
          end do
        end do
      else
        do a = 1, alt
          do i = lat - gcn + 1, lat
            do j = long - gcn + 1, long
              out_array(p) = v_in(j, i, a)
              p = p + 1
            end do
          end do
        end do
      end if

    case (ab_northwest)
      if (pole) then
        do a = 1, alt
          do i = lat, lat - gcn + 1, -1
            do j = 1, gcn
              out_array(p) = v_in(j, i, a)
              p = p + 1
            end do
          end do
        end do
      else
        do a = 1, alt
          do i = lat - gcn + 1, lat
            do j = 1, gcn
              out_array(p) = v_in(j, i, a)
              p = p + 1
            end do
          end do
        end do
      end if

    case (ab_southeast)
      if (pole) then
        do a = 1, alt
          do i = gcn, 1, -1
            do j = long - gcn + 1, long
              out_array(p) = v_in(j, i, a)
              p = p + 1
            end do
          end do
        end do
      else
        do a = 1, alt
          do i = 1, gcn
            do j = long - gcn + 1, long
              out_array(p) = v_in(j, i, a)
              p = p + 1
            end do
          end do
        end do
      end if

    case (ab_southwest)
      if (pole) then
        do a = 1, alt
          do i = gcn, 1, -1
            do j = 1, gcn
              out_array(p) = v_in(j, i, a)
              p = p + 1
            end do
          end do
        end do
      else
        do a = 1, alt
          do i = 1, gcn
            do j = 1, gcn
              out_array(p) = v_in(j, i, a)
              p = p + 1
            end do
          end do
        end do
      end if

    end select

  end subroutine AB_array3_gc_pack

  subroutine AB_array3_gc_unpack(long, lat, alt, gcn, v_out, dir, p, in_array)
    use AB_module
    integer, intent(in) :: long, lat, alt, gcn
    integer, intent(in) :: dir
    integer, intent(inout) :: p
    real, dimension(1 - gcn:, 1 - gcn:, 1:), intent(inout) :: v_out
    real, dimension(:), intent(in) :: in_array
    integer :: i, j, a

    select case (dir)

    case (ab_north)
      do a = 1, alt
        do i = lat + 1, lat + gcn
          do j = 1, long
            v_out(j, i, a) = in_array(p)
            p = p + 1
          end do
        end do
      end do

    case (ab_south)
      do a = 1, alt
        do i = 1 - gcn, 0
          do j = 1, long
            v_out(j, i, a) = in_array(p)
            p = p + 1
          end do
        end do
      end do

    case (ab_east)
      do a = 1, alt
        do i = 1, lat
          do j = long + 1, long + gcn
            v_out(j, i, a) = in_array(p)
            p = p + 1
          end do
        end do
      end do

    case (ab_west)
      do a = 1, alt
        do i = 1, lat
          do j = 1 - gcn, 0
            v_out(j, i, a) = in_array(p)
            p = p + 1
          end do
        end do
      end do

    case (ab_northeast)
      do a = 1, alt
        do i = lat + 1, lat + gcn
          do j = long + 1, long + gcn
            v_out(j, i, a) = in_array(p)
            p = p + 1
          end do
        end do
      end do

    case (ab_northwest)
      do a = 1, alt
        do i = lat + 1, lat + gcn
          do j = 1 - gcn, 0
            v_out(j, i, a) = in_array(p)
            p = p + 1
          end do
        end do
      end do

    case (ab_southeast)
      do a = 1, alt
        do i = 1 - gcn, 0
          do j = long + 1, long + gcn
            v_out(j, i, a) = in_array(p)
            p = p + 1
          end do
        end do
      end do

    case (ab_southwest)
      do a = 1, alt
        do i = 1 - gcn, 0
          do j = 1 - gcn, 0
            v_out(j, i, a) = in_array(p)
            p = p + 1
          end do
        end do
      end do

    end select

  end subroutine AB_array3_gc_unpack

  subroutine AB_array4_gc_add_size(long, lat, alt, eta, gcn, size)
    use AB_module
    integer, intent(inout) :: size(ab_num_nbrs)
    integer, intent(in) :: long, lat, alt, eta, gcn

    size(ab_north) = size(ab_north) + long*gcn*alt*eta
    size(ab_south) = size(ab_south) + long*gcn*alt*eta
    size(ab_east) = size(ab_east) + lat*gcn*alt*eta
    size(ab_west) = size(ab_west) + lat*gcn*alt*eta
    size(ab_northeast) = size(ab_northeast) + gcn*gcn*alt*eta
    size(ab_northwest) = size(ab_northwest) + gcn*gcn*alt*eta
    size(ab_southeast) = size(ab_southeast) + gcn*gcn*alt*eta
    size(ab_southwest) = size(ab_southwest) + gcn*gcn*alt*eta
  end subroutine AB_array4_gc_add_size

  subroutine AB_array4_gc_pack(long, lat, alt, eta, gcn, v_in, dir, pole, p, out_array)
    use AB_module
    integer, intent(in) :: long, lat, alt, eta, gcn
    integer, intent(in) :: dir
    logical, intent(in) :: pole
    integer, intent(inout) :: p
    real, dimension(1 - gcn:, 1 - gcn:, 1:, 1:), intent(in) :: v_in
    real, dimension(:), intent(inout) :: out_array
    integer :: i, j, a, e

    select case (dir)
    case (ab_north)
      if (pole) then
        do e = 1, eta
          do a = 1, alt
            do i = lat, lat - gcn + 1, -1
              do j = 1, long
                out_array(p) = v_in(j, i, a, e)
                p = p + 1
              end do
            end do
          end do
        end do
      else
        do e = 1, eta
          do a = 1, alt
            do i = lat - gcn + 1, lat
              do j = 1, long
                out_array(p) = v_in(j, i, a, e)
                p = p + 1
              end do
            end do
          end do
        end do
      end if

    case (ab_south)
      if (pole) then
        do e = 1, eta
        do a = 1, alt
          do i = gcn, 1, -1
            do j = 1, long
              out_array(p) = v_in(j, i, a, e)
              p = p + 1
            end do
          end do
        end do
        end do
      else
        do e = 1, eta
        do a = 1, alt
          do i = 1, gcn
            do j = 1, long
              out_array(p) = v_in(j, i, a, e)
              p = p + 1
            end do
          end do
        end do
        end do
      end if

    case (ab_east)
      do e = 1, eta
      do a = 1, alt
        do i = 1, lat
          do j = long - gcn + 1, long
            out_array(p) = v_in(j, i, a, e)
            p = p + 1
          end do
        end do
      end do
      end do

    case (ab_west)
      do e = 1, eta
      do a = 1, alt
        do i = 1, lat
          do j = 1, gcn
            out_array(p) = v_in(j, i, a, e)
            p = p + 1
          end do
        end do
      end do
      end do

    case (ab_northeast)
      if (pole) then
        do e = 1, eta
        do a = 1, alt
          do i = lat, lat - gcn + 1, -1
            do j = long - gcn + 1, long
              out_array(p) = v_in(j, i, a, e)
              p = p + 1
            end do
          end do
        end do
        end do
      else
        do e = 1, eta
        do a = 1, alt
          do i = lat - gcn + 1, lat
            do j = long - gcn + 1, long
              out_array(p) = v_in(j, i, a, e)
              p = p + 1
            end do
          end do
        end do
        end do
      end if

    case (ab_northwest)
      if (pole) then
        do e = 1, eta
        do a = 1, alt
          do i = lat, lat - gcn + 1, -1
            do j = 1, gcn
              out_array(p) = v_in(j, i, a, e)
              p = p + 1
            end do
          end do
        end do
        end do
      else
        do e = 1, eta
        do a = 1, alt
          do i = lat - gcn + 1, lat
            do j = 1, gcn
              out_array(p) = v_in(j, i, a, e)
              p = p + 1
            end do
          end do
        end do
        end do
      end if

    case (ab_southeast)
      if (pole) then
        do e = 1, eta
        do a = 1, alt
          do i = gcn, 1, -1
            do j = long - gcn + 1, long
              out_array(p) = v_in(j, i, a, e)
              p = p + 1
            end do
          end do
        end do
        end do
      else
        do e = 1, eta
        do a = 1, alt
          do i = 1, gcn
            do j = long - gcn + 1, long
              out_array(p) = v_in(j, i, a, e)
              p = p + 1
            end do
          end do
        end do
        end do
      end if

    case (ab_southwest)
      if (pole) then
        do e = 1, eta
        do a = 1, alt
          do i = gcn, 1, -1
            do j = 1, gcn
              out_array(p) = v_in(j, i, a, e)
              p = p + 1
            end do
          end do
        end do
        end do
      else
        do e = 1, eta
        do a = 1, alt
          do i = 1, gcn
            do j = 1, gcn
              out_array(p) = v_in(j, i, a, e)
              p = p + 1
            end do
          end do
        end do
        end do
      end if

    end select

  end subroutine AB_array4_gc_pack

  subroutine AB_array4_gc_unpack(long, lat, alt, eta, gcn, v_out, dir, p, in_array)
    use AB_module
    integer, intent(in) :: long, lat, alt, eta, gcn
    integer, intent(in) :: dir
    integer, intent(inout) :: p
    real, dimension(1 - gcn:, 1 - gcn:, 1:, 1:), intent(inout) :: v_out
    real, dimension(:), intent(in) :: in_array
    integer :: i, j, a, e

    select case (dir)

    case (ab_north)
      do e = 1, eta
      do a = 1, alt
        do i = lat + 1, lat + gcn
          do j = 1, long
            v_out(j, i, a, e) = in_array(p)
            p = p + 1
          end do
        end do
      end do
      end do

    case (ab_south)
      do e = 1, eta
      do a = 1, alt
        do i = 1 - gcn, 0
          do j = 1, long
            v_out(j, i, a, e) = in_array(p)
            p = p + 1
          end do
        end do
      end do
      end do

    case (ab_east)
      do e = 1, eta
      do a = 1, alt
        do i = 1, lat
          do j = long + 1, long + gcn
            v_out(j, i, a, e) = in_array(p)
            p = p + 1
          end do
        end do
      end do
      end do

    case (ab_west)
      do e = 1, eta
      do a = 1, alt
        do i = 1, lat
          do j = 1 - gcn, 0
            v_out(j, i, a, e) = in_array(p)
            p = p + 1
          end do
        end do
      end do
      end do

    case (ab_northeast)
      do e = 1, eta
      do a = 1, alt
        do i = lat + 1, lat + gcn
          do j = long + 1, long + gcn
            v_out(j, i, a, e) = in_array(p)
            p = p + 1
          end do
        end do
      end do
      end do

    case (ab_northwest)
      do e = 1, eta
      do a = 1, alt
        do i = lat + 1, lat + gcn
          do j = 1 - gcn, 0
            v_out(j, i, a, e) = in_array(p)
            p = p + 1
          end do
        end do
      end do
      end do

    case (ab_southeast)
      do e = 1, eta
      do a = 1, alt
        do i = 1 - gcn, 0
          do j = long + 1, long + gcn
            v_out(j, i, a, e) = in_array(p)
            p = p + 1
          end do
        end do
      end do
      end do

    case (ab_southwest)
      do e = 1, eta
      do a = 1, alt
        do i = 1 - gcn, 0
          do j = 1 - gcn, 0
            v_out(j, i, a, e) = in_array(p)
            p = p + 1
          end do
        end do
      end do
      end do

    end select

  end subroutine AB_array4_gc_unpack

end module AB_XFER_array_util

