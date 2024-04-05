
FINAL PARAMETERS:

+ no. stations          :  109
+ dof                   :  215
+ varf                  :  32.3129474588
+ sigma0                :  5.68444785874

+ Estimated plate model :  1.50379 +/-  0.00417  1.18346 +/-  0.00401  1.20716 +/-  0.00370 mas
+ Estimated plate model :   0.6285 /Ma  32.2446 deg  45.0910 deg rate latitude longitude


Execute with the following command

./sinex2epoch.exe -e20:001 -s32.312947 -E1.503791:1.183464:1.207157 \
-q4.078805e-22:-1.281398e-22:9.127072e-23:3.774281e-22:-7.771154e-23:3.217797e-22 XVSOLFIN14_ver05.SNX.AUSONLY

+ project station coordinates to their individual mean epoch
+ reset velocities to plate model for those stations with a horizontal velocity difference > 0.0015 m/y
+ reset velocities to plate model for those stations with a horizontal velocity uncertainty > 0.0015 m/y
+ reset velocities uncertainties to 0.0005 m/y XYZ
+ project coordinates to the epoch of 20:001
+ reset velocities to plate model for those stations with a horizontal velocity difference > 0 m/y
+ reset velocities uncertainties to 0.0005 m/y XYZ

