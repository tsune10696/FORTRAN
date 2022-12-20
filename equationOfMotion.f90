program equationOfMotion
implicit none
!-------------------------------
!  equationOfMotion.f90
!  tunaProducts(2022.12.19)
!  (kは弾性）
!-------------------------------

!---文字の定義------
integer :: i, n, step
real(8) :: ydd_EqNS(20000), ydd_EqEW(20000), ydd_EqUD(20000) !地震波(gal)
real(8) :: y(20000), yd(20000), ydd(20000), time(20000) !変位(m)、速度(m/s)、加速度(m/s/s)、時間(s)
real(8) :: beta, dt, omega, T, h, m, k, c, pi, delta_p,epsilon, delta_y

!---ファイルを開く--
open(11,file="wave/JMAKobeWave.dat")
open(12,file="output/outputKobe.txt")

!--ファイルの読込-----
do i = 1,17996,1
read(11,"(3(f6.3))",end=9999) ydd_EqNS(i), ydd_EqEW(i), ydd_EqUD(i)
end do
9999 continue

!---数値の入力--------
beta = 0.25                       !β法(平均加速度法なので1/4)
dt = 0.02                         !地震波の時間刻み（単位はs）
h = 0.05                          !減衰定数(5%の減衰を与える)
epsilon = 1.0 * (10.0 ** (-5.0))  !収れん判定に用いる値（0に近い数値）
step = 17996                      !地震波のデータ数
pi = 2.0 * asin(1.0)              !円周率
m = 2.0 * 10.0**5                 !質量(単位はkg、平米1tの建物で200m2を想定)
k = 4.0 * pi**2 * m               !降伏前のバネ定数(単位はN/m、周期が1.0秒になるように値調整)
T = 2.0 * pi * (m / k)**0.5       !一次固有周期(s)
omega = 2 * pi / T                !固有振動数(1/s)
c = 2.0 * h * k / omega           !粘性係数(剛性比例型)
!c = 2.0 * h * omega * m          !粘性係数(質量比例型)

!---初期条件----------
ydd_EqNS(1) = 0.0             !時刻0.00秒の時の地動加速度(m/s/s)
ydd_EqEW(1) = 0.0             !時刻0.00秒の時の地動加速度(m/s/s)
ydd_EqUD(1) = 0.0             !時刻0.00秒の時の地動加速度(m/s/s)
y(1) = 0.0                    !時刻0.00秒の時の応答変位(m)
yd(1) = 0.0                   !時刻0.00秒の時の応答速度(m/s)
ydd(1) = 0.0                  !時刻0.00秒の時の応答加速度(m/s/s)

write(12,*)'# step,responseAcc(absolute)(m/s/s), inputWave(m/s/s)' !出力データのラベルを作成

!---次ステップの値を求める----------
do n = 1, step, 1
    y(n+1)   = y(n)                                                                                          !変位を仮定
    ydd(n+1) = (y(n+1) - y(n)) / (beta * dt**2) - yd(n) / (beta * dt) - ydd(n) * (1.0 / (2.0 * beta) - 1.0)  !加速度を求める
    yd(n+1)  = yd(n) + 0.50 * (ydd(n) + ydd(n+1)) * dt                                                       !速度を求める

    !---収れん判定（収れんしていない場合にはループの計算を続ける）---
    do
        delta_p = (-1.0) * m * ydd_EqNS(n+1) * 0.01 - (m*ydd(n+1) + c*yd(n+1) + k*y(n+1)) !不釣合力を計算
        
        if (abs(delta_p).lt.epsilon) exit
            delta_y = delta_p / (k + m / (beta * dt**2) + c /(2.0*beta*dt))
            y(n+1) = y(n+1) + delta_y                                                                               !変位修正
            ydd(n+1) = (y(n+1) - y(n)) / (beta * dt**2) - yd(n) / (beta * dt) - ydd(n) * (1.0 / (2.0 * beta) - 1.0) !加速度を求める
            yd(n+1)  = yd(n) + 0.50 * (ydd(n) + ydd(n+1)) * dt                                                      !速度を求める
    end do

!---出力-----------------
    write(12,*) n+1 , ydd(n+1) + ydd_EqNS(n+1) * 0.01 , ydd_EqNS(n+1) * 0.01
end do

!---ファイルを閉じる----
close(11)
close(12)

stop
end program equationOfMotion
