program blobInfo
implicit none

integer ( kind = 4 ), parameter :: N = 126
integer ( kind = 4 ), parameter :: xBin = 100

real*8 blobLimX(N), blobLimY(N), blobMidX(N), blobMidY(N), polyArea(N)
real*8 blobWidthX(N,10), blobWidthY(N,10), blobCenterX(N,10), blobCenterY(N,10), blobArea(N,10)
real*8 velocityX(N,N,10), velocityY(N,N,10), velocity(N,N,10)
real*8 diffWidthX(N,N,10,10), diffWidthY(N,N,10,10), diffCenterX(N,N,10,10), diffCenterY(N,N,10,10), diffArea(N,N,10,10)
real*8 blobLimXSum, blobLimYSum, polyAreaSum
real*8 dx, dy, vy, Sx, Sy, dA, diffCenterXBlob, diffCenterYBlob
real*8 xLeft, xRight, xWidth, deltaX, yTop, yBottom, simLimX, simLimY, simArea
real*8 corrPosition(N*N), corrPos(N*N), corrVal(N*N), corrValue(N*N), freqPos(N), freqVal(N), freqValue(N)
real*8 binCorr(xBin), binFreq(xBin),fFil(xBin)
integer i, frame(N), G(N), iter, iterBack, maxIter, nBlob, mBlob, maxBlob, time(N), correlation, count, flag
integer corrBin, freqBin, tempBin, corrMax, freqMax, tempMax
real*8 rho_s,C_s,microsec,second

rho_s   = 3.58*1e-4
C_s     = 6.3*1e4
microsec= 1.0d0
second  = microsec*10**6

xLeft   = 2.147
xRight  = 2.187
xWidth  = xRight - xLeft
deltaX  = xWidth/xBin
yTop    = 0.033
yBottom = -0.033
simLimX = xRight - xLeft
simLimY = yTop - yBottom
simArea = simLimX * simLimY

open(unit=10,file='blob_size.txt',status='unknown')
open(unit=20,file='blob_details.dat',status='unknown')

blobLimXSum = 0.0d0
blobLimYSum = 0.0d0
polyAreaSum = 0.0d0

do i = 1, N
  read(10,*) frame(i), G(i), blobLimX(i), blobLimY(i), blobMidX(i), blobMidY(i), polyArea(i)
  blobLimXSum = blobLimXSum + blobLimX(i)
  blobLimYSum = blobLimYSum + blobLimY(i)
  polyAreaSum = polyAreaSum + polyArea(i)
enddo

write(*,*) blobLimXSum/dfloat(N), blobLimYSum/dfloat(N), polyAreaSum/dfloat(N), polyAreaSum*100.0d0/(dfloat(N) * simArea)

!call exchange (blobLimX, blobLimX, N)
!call exchange (blobLimY, blobLimX, N)

!do i = N,N-5,-1
!  write(*,*) blobLimY(i)/rho_s
!enddo

i                        = 1
iter                     = 1
nBlob                    = 1
maxBlob                  = 1
time(iter)               = frame(i)
blobWidthX(iter, nBlob)  = blobLimX(i)
blobWidthY(iter, nBlob)  = blobLimY(i)
blobCenterX(iter, nBlob) = blobMidX(i)
blobCenterY(iter, nBlob) = blobMidY(i)
blobArea(iter, nBlob)    = polyArea(i)

do i = 2, N
  if (frame(i) - frame(i-1) /= 0) then
    iter                     = iter + 1
    nBlob                    = 1
    time(iter)               = frame(i)
    blobWidthX(iter, nBlob)  = blobLimX(i)
    blobWidthY(iter, nBlob)  = blobLimY(i)
    blobCenterX(iter, nBlob) = blobMidX(i)
    blobCenterY(iter, nBlob) = blobMidY(i)
    blobArea(iter, nBlob)    = polyArea(i)
  elseif (frame(i) - frame(i-1) == 0) then
    iter                     = iter
    nBlob                    = nBlob + 1
    blobWidthX(iter, nBlob)  = blobLimX(i)
    blobWidthY(iter, nBlob)  = blobLimY(i)
    blobCenterX(iter, nBlob) = blobMidX(i)
    blobCenterY(iter, nBlob) = blobMidY(i)
    blobArea(iter, nBlob)    = polyArea(i)
  endif
  if (nBlob > maxBlob) maxBlob = nBlob
!  print*, i, frame(i), iter, time(iter), nBlob, maxBlob
enddo

maxIter = iter

tempBin = 0
freqBin = 0

do iter = 2, maxIter
  do nBlob = 1, maxBlob
    correlation = 0
    flag = 0
    do iterBack = iter-1, 1, -1
      do mBlob = 1, maxBlob
        diffWidthX(iter, iterBack, nBlob, mBlob)  = blobWidthX(iter, nBlob)  - blobWidthX(iterBack, mBlob)
        diffWidthY(iter, iterBack, nBlob, mBlob)  = blobWidthY(iter, nBlob)  - blobWidthY(iterBack, mBlob)
        diffCenterX(iter, iterBack, nBlob, mBlob) = blobCenterX(iter, nBlob) - blobCenterX(iterBack, mBlob)
        diffCenterY(iter, iterBack, nBlob, mBlob) = blobCenterY(iter, nBlob) - blobCenterY(iterBack, mBlob)
        diffArea(iter, iterBack, nBlob, mBlob)    = blobArea(iter, nBlob)    - blobArea(iterBack, mBlob)
        diffCenterXBlob                           = abs(diffCenterX(iter, iterBack, nBlob, mBlob))*second/C_s
        diffCenterYBlob                           = abs(diffCenterY(iter, iterBack, nBlob, mBlob))*second/C_s
        if (diffCenterXBlob < 0.2 .and. diffCenterXBlob > 1e-10 .and. diffCenterYBlob < 0.2 .and. diffCenterYBlob > 1e-10) then
          correlation = correlation + 1
          flag = 1
          velocityX(iter, iterBack, nBlob) = abs(diffCenterX(iter, iterBack, nBlob, mBlob))*second/C_s
          velocityY(iter, iterBack, nBlob) = abs(diffCenterY(iter, iterBack, nBlob, mBlob))*second/C_s
          velocity(iter, iterBack, nBlob)  = dsqrt(velocityX(iter, iterBack, nBlob)*velocityX(iter, iterBack, nBlob) &
                                                 + velocityY(iter, iterBack, nBlob)*velocityY(iter, iterBack, nBlob))
          if (iterBack == iter-1) write(20,*) iter, nBlob, mBlob, blobWidthY(iter, nBlob)/rho_s, velocityX(iter, iterBack, nBlob)
          if (iterBack < iter-1 .and. correlation == iter - iterBack) then
            tempBin = tempBin + 1
            corrPos(tempBin) = blobCenterX(iter, nBlob)
            corrVal(tempBin) = dfloat(correlation)
            write(30,*) blobCenterX(iter, nBlob), dfloat(correlation)!, blobWidthX(iter, nBlob)/rho_s, velocityX(iter, iterBack, nBlob)
          endif
        endif
      enddo
    if (iterBack == iter-1 .and. flag == 0 .and. blobCenterX(iter, nBlob) > xLeft) then
      freqBin = freqBin + 1
      freqPos(freqBin) = blobCenterX(iter, nBlob)
      freqVal(freqBin) = time(iter)-time(1)
      !print*, blobCenterX(iter, nBlob), time(iter)-time(1)
    endif
    enddo
  enddo
enddo

tempMax = tempBin
freqMax = freqBin

freqValue(1) = freqVal(1)
corrValue(1) = corrVal(1)

corrBin = 0

do freqBin = 2, freqMax
  freqValue(freqBin) = freqVal(freqBin) - freqVal(freqBin-1)
  !print*, freqBin, freqPos(freqBin), freqValue(freqBin)
enddo

call exchange (freqPos, freqValue, freqMax)

do freqBin = 1, freqMax
  write(35,*) freqBin, freqPos(freqBin), freqValue(freqBin)
  i = int((freqPos(freqBin)-xLeft)/deltaX)
  binFreq(i) = binFreq(i) + freqValue(freqBin)
enddo

do tempBin = 2, tempMax
  if (corrPos(tempBin) /= corrPos(tempBin-1)) then
    corrBin = corrBin + 1
    corrPosition(corrBin) = corrPos(tempBin-1)
    corrValue(corrBin) = corrVal(tempBin-1)
    !print*, corrBin, tempBin, corrPosition(corrBin), corrValue(corrBin)
  endif
enddo

corrMax = corrBin

call exchange (corrPosition, corrValue, corrMax)

do corrBin = 1, corrMax
  write(40,*) corrBin, corrPosition(corrBin), corrValue(corrBin)
  i = int((corrPosition(corrBin)-xLeft)/deltaX)
  binCorr(i) = binCorr(i) + corrValue(corrBin)
enddo

do i = 1, xBin
  if (abs(binFreq(i)) > 1e-5 .and. abs(binFreq(i)) < 1e5) then
    fFil(i) = binCorr(i)/binFreq(i)
    write(45,*) xLeft+(dfloat(i)*deltaX), fFil(i), binCorr(i), 1.0d0/binFreq(i) ! This is auto-correlation time for blobs!
  endif
enddo

contains

subroutine exchange (b, s, n)

integer, intent (in) ::  n
real*8, dimension (n), intent (inout) :: b
real*8, dimension (n), intent (inout) :: s
real*8 :: c
integer*8 ::  i, j, d

do j = n-1, 1, -1
  do i = 1, j, 1
    if (b(i) >= b(i+1)) then
    c=b(i)
    b(i)=b(i+1)
    b(i+1)=c
    d = s(i)
    s(i) = s(i+1)
    s(i+1) = d
    endif
  end do
end do

end subroutine exchange

end program blobInfo
