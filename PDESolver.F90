!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! *
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! *
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
SUBROUTINE PDESolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
     implicit none
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
   integer, parameter :: maxrows = 1000000
   integer :: a, b, nsteps
   integer :: ti_noel, ti_npt, ti_ntens
   real(8) :: ti_time, ti_v1, ti_v2, ti_v3, ti_v4
   logical :: ti_keep
   integer :: i, j, nrows, ios,k,nume,numip,numts, row
   integer :: noel(maxrows), npt(maxrows), ntens(maxrows)
   real(8) :: tval(maxrows), beta, density, Cp
   real(8) :: v1(maxrows), v2(maxrows), v3(maxrows), v4(maxrows)
   real(8):: dPlastic,pwi,cumwork,heatinc,cumheat,Trise,cumT
   real(8) :: v5(maxrows), v6(maxrows), v7(maxrows), v8(maxrows)
   real(8) :: v9(maxrows), v10(maxrows), v11(maxrows)

   logical :: keep(maxrows)

   character(len=200) :: line, header

   nrows = 0

   open(unit=10,file='PDEOutput.csv',status='old',action='read')
   open(unit=20,file='converged.csv',status='replace')
      read(10,'(A)',iostat=ios) header
   !-----------------------------------------------
   ! Read CSV file line by line
   !-----------------------------------------------
   do
      read(10,'(A)',iostat=ios) line
      if (ios /= 0) exit

      ! skip blank lines
      if (len_trim(line) == 0) cycle

      nrows = nrows + 1

      read(line,*) tval(nrows), noel(nrows), npt(nrows), ntens(nrows), &
                   v1(nrows), v2(nrows), v3(nrows), v4(nrows)

      keep(nrows) = .true.
   end do

   close(10)

   !-----------------------------------------------
   ! Mark earlier duplicates as false
   ! Keep only LAST occurrence
   !-----------------------------------------------
   do i = 1, nrows-1
      do j = i+1, nrows
         if ( tval(i)  == tval(j)  .and. &
              noel(i)  == noel(j)  .and. &
              npt(i)   == npt(j)   .and. &
              ntens(i) == ntens(j) ) then

            keep(i) = .false.
            exit
         end if
      end do
   end do

!-------------------------------------------------
! Sort retained rows before writing:
! order = noel, npt, ntens
!-------------------------------------------------

do a = 1, nrows-1
   do b = a+1, nrows

      if (keep(a) .and. keep(b)) then

if ( noel(b) < noel(a) .or. &
   (noel(b) == noel(a) .and. npt(b) < npt(a)) .or. &
   (noel(b) == noel(a) .and. npt(b) == npt(a) .and. ntens(b) < ntens(a)) .or. &
   (noel(b) == noel(a) .and. npt(b) == npt(a) .and. ntens(b) == ntens(a) &
        .and. tval(b) < tval(a)) ) then

            !--- swap integers
            ti_noel  = noel(a);  noel(a)  = noel(b);  noel(b)  = ti_noel
            ti_npt   = npt(a);   npt(a)   = npt(b);   npt(b)   = ti_npt
            ti_ntens = ntens(a); ntens(a) = ntens(b); ntens(b) = ti_ntens

            !--- swap reals
            ti_time = tval(a); Tval(a)=tval(b); tval(b)=ti_time
            ti_v1   = v1(a);   v1(a)=v1(b);     v1(b)=ti_v1
            ti_v2   = v2(a);   v2(a)=v2(b);     v2(b)=ti_v2
            ti_v3   = v3(a);   v3(a)=v3(b);     v3(b)=ti_v3
            ti_v4   = v4(a);   v4(a)=v4(b);     v4(b)=ti_v4

            !--- swap logical
            ti_keep = keep(a); keep(a)=keep(b); keep(b)=ti_keep

         end if
      end if

   end do
end do

  ! Calculate PDE for each noel, npt, ntens
density = 7850.0
beta = 0.90
Cp = 470.0
nume = noel(nrows)
numip = npt(nrows)
numts = ntens(nrows)
write(*,*) nrows
write(*,*) 'Elements:',nume,' max npt:',numip,' max ntens:',numts
write(*,*) 'Calculating PDE'

   !-----------------------------------------------
   ! Write final converged rows 
   !-----------------------------------------------
j = 0

write(20,*) trim(header) // ',dPlastic,pwi,cumwork,heatinc,cumH,Trise,cumT'
   do i = 1, nrows
      if (keep(i)) then
         write(20,'(F12.5,",",I8,",",I8,",",I8,",",11(1PE16.8,","))') &
             tval(i), noel(i), npt(i), ntens(i), &
             v1(i), v2(i), v3(i), v4(i),&
v5(i),v6(i),v7(i),v8(i),v9(i),v10(i),v11(i)
     j = j + 1
      end if
   end do
   close(20)



  ! Calculate PDE for each noel, npt, ntens
density = 7850.0
beta = 0.90
Cp = 470.0
write(*,*) 'Converged Rows:',j
nrows = j
Row = 1
beta = 0.9
nsteps = j/(nume*numip*numts)
! Reload Data 
  open(unit=20,file='converged.csv',status='old')
read(20,fmt='(A200)') header
do i = 1, nrows
 read(20,*) tval(i), noel(i), npt(i), ntens(i), &
             v1(i), v2(i), v3(i), v4(i),&
v5(i),v6(i),v7(i),v8(i),v9(i),v10(i),v11(i)
end do
do a = 1,nsteps
do i = 1,nume
 Do j= 1,numip
  Do  k = 1,numts
    if (tval(row)  .eq.  0.0) then
        dPlastic = 0.0
        pwi = 0.0
        cumwork = 0.0
        heatinc = 0.0
        cumheat = 0.0
        Trise = 0.0
        CumT = 0.0
 v5(row) = dPlastic
v6(row) = pwi
v7(row) = cumwork
v8(row) =heatinc
v9(row) = cumheat
v10(row) = Trise
v11(row) = cumT
      else
        dPlastic = v4(row) - v4(row-1)
       pwi =dPlastic *  (v1(row)+v1(row-1))/2
       cumwork = v7(row-1) + pwi
       heatinc = pwi*1000000.0*beta
       cumheat = v9(row-1) + heatinc
       Trise = heatinc/(density*Cp)
       cumT = v11(row-1) + Trise
 v5(row) = dPlastic
v6(row) = pwi
v7(row) = cumwork
v8(row) =heatinc
v9(row) = cumheat
v10(row) = Trise
v11(row) = cumT
     endif
row = row +1
    End do
  End do
End do
end do
   !-----------------------------------------------
   ! Write final data
   !-----------------------------------------------


write(20,*) trim(header) // ',dPlastic,pwi,cumwork,heatinc,cumH,Trise,cumT'
   do i = 1, nrows
              write(20,'(F12.5,",",I8,",",I8,",",I8,",",11(1PE16.8,","))') &
             tval(i), noel(i), npt(i), ntens(i), &
             v1(i), v2(i), v3(i), v4(i),&
v5(i),v6(i),v7(i),v8(i),v9(i),v10(i),v11(i)
   
   end do
   close(20)
   print *, 'Finished. Output written to converged.csv'

END SUBROUTINE PDESolver
!------------------------------------------------------------------------------
