TITLE 'Marauder 2'
COORDINATES cartesian1 { coordinate system, 1D,2D,3D, etc }
VARIABLES { system variables }
rx(threshold=1e-5)
vx(threshold=1e-5)
ry(threshold=1e-5)
vy(threshold=1e-5)
SELECT { method controls }
ngrid = 11 !accuracy vs speed
DEFINITIONS { parameter definitions }
! Initial mass
mi=19.36
! Fuel mass
mfueli=3.547
! Motor firing time
tfuel=3.61
! Launch altitude
alt = 1400
! Launch angle
launchangle = 2
delta = (ry^2)/90000
theta = if (vy>0) then launchangle + delta else 0
! Main deploy altitude
altp = 396.24 !1300 ft
! Wind
ws1 = 8.9
alt1 = 1000
ws2 = 8.9
alt2 = 2000
ws3 = 8.9
alt3 = 3000
ws4 = 8.9
alt4 = 4000
! Parachute
Ddrogue = 3
Adrogue = 3.14*(Ddrogue/6.562)^2
CDdrogue = 0.40
Dmain = 13
Amain = 3.14*(Dmain/6.562)^2
CDmain = 0.39
Ayairframe = 0.0135
CDnose = 0.40
! X metrics
Areax=0.347 
CDx=0.38 !Clarify
! Force of thrust
Fthrust = if (t<0.2) then 2400 else if (t<1) then 2500 else if (t<2.6) then 2250 else if (t<2.8) then 2000 else if (t<3.3) then 750 else if (t<3.61) then 150 else 0 !Rough thrust curve
Fthrustx = -Fthrust*sin(theta*0.01745)
Fthrusty = Fthrust*cos(theta*0.01745)
! Mass of fuel over time
mfuelgone = if (t<tfuel) then t*mfueli/tfuel else mfueli
m = mi - mfuelgone
g=9.80665
v = sqrt(vx^2+vy^2)
vrelx = if (0<ry and ry<alt1) then vx-ws1 else if (alt1<ry and ry<alt2) then vx-ws2 else if (alt2<ry and ry<alt3) then vx-ws3 else if (alt3<ry and ry<alt4) then vx-ws4 else vx
vrely = vy
vrel = sqrt(vrelx^2+vrely^2)
! Air density
rho =  (-(((ry+alt)/1000)-44.3308)/42.2665)^(7418/1743)
Areay=if (vy>0) then Ayairframe else if (vy<0 and ry>altp) then Adrogue else if (vy<0 and ry<altp) then Amain else Ayairframe
CDy = if (vy>0) then CDnose else if (vy<0 and ry>altp) then CDdrogue else if (vy<0 and ry<altp) then CDmain else CDnose
Fdragy = -.5*rho*Areay*CDy*vrel^2*vrely/(vrel+1e-6)
Fdragx = -.5*rho*Areax*CDx*vrel^2*vrelx/(vrel+1e-6)
Fgrav = -m*g
ax = (Fdragx+Fthrustx)/m
ay = (Fthrusty+Fdragy)/m - g
INITIAL VALUES
vy = 0.01
vx = 0.01
EQUATIONS { PDE's, one for each variable }
rx: dt(rx) = vx
vx: dt(vx) = ax
ry: dt(ry) = vy
vy: dt(vy) = ay
BOUNDARIES { The domain definition }
REGION 1 START(0) LINE TO (1)
TIME 0 TO 2000 halt(ry<0) { if time dependent }
PLOTS { save result displays }
for t= 0 by endtime/50 to endtime
history(ry, vy, ay) at(0)
Export file="Output.csv"
history(rx, vx, ax) at (0)
history(ry) at (0) vs rx
history(Fthrustx, Fthrusty) at (0)
history(rho) at(0)
history(theta) at (0) vs rx
report t
summary
report val(rho,0)
report val(vy, 0)
report val(rx, 0)
END
