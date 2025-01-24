## ===
# DESCRIPTION
# This file contains useful functions for common geodetic operations.
# 
# Functions:    - coordinate_transformation()
#               - plate_motion_model()
#               - get_velocities_APMM()
#               - xyz2llh()
#               - xyz2enu()
#               - xyz2enu_vcv()
#               - llh2xyz()
#               - enu2xyz()
#               - enu2xyz_vcv()
#               - enu2los()
#               - enu2los_vcv()
#               - doy2date()
#               - dd2ddm()
#               - dd2dms()
#               - dms2dd()
#               - vincentys_inverse()

"""
    coordinate_transformation(xi, yi, zi, convention="AGRS", tx=0.0, ty=0.0, tz=0.0, d=0.0, rx=0.0, ry=0.0, rz=0.0, t1=0.0, t0=0.0, tx_r=0.0, ty_r=0.0, tz_r=0.0, d_r=0.0, rx_r=0.0, ry_r=0.0, rz_r=0.0)

Coordinate transformation from one reference frame to another. 
- If no rates and epochs are given, it works as a 7 parameter transformation (in the propogation equations, the rate and elapsed time become zero).
- If rates and epochs are given, it works as a 14-parameter transformation. 
    
# NOTE:
- This function requires the rotation matrix sign convention to be set to "AGRS" (default) or "IERS"
- For AGRS convention, see GDA2020 Technical Manual (Section 2.2.1). 
- For IERS convention, see GNSS Handbook Eq. 2.29., or Chapter 4 of the IERS Conventions (https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html). 

# Arguments 
- `xi::Float`:          X component of input cartesian coordinates that are being transformed (metres).
- `yi::Float`:          Y component of input cartesian coordinates that are being transformed (metres).
- `zi::Float`:          Z component of input cartesian coordinates that are being transformed (metres).
- `convention::String`: Sign convention for roation matrix. "AGRS" (default) or "IERS".
- `tx::Float`:          Translation parameter of the X component (metres).
- `ty::Float`:          Translation parameter of the Y component (metres).
- `tz::Float`:          Translation parameter of the Z component (metres).
- `d::Float`:           Scale parameter (ppb).
- `rx::Float`:          Rotation parameter of the X component (mas).
- `ry::Float`:          Rotation parameter of the Y component (mas).
- `rz::Float`:          Rotation parameter of the Z component (mas).
- `t1::Float`:          Epoch of output TRF coordinates (decimal years).
- `t0::Float`:          Reference epoch of transformation (decimal years).
- `tx_r::Float`:        Rate of change for the translation parameter of the X component (metres per year).
- `ty_r::Float`:        Rate of change for the translation parameter of the Y component (metres per year).
- `tz_r::Float`:        Rate of change for the translation parameter of the Z component (metres per year).
- `d_r::Float`:         Rate of change for the scale parameter (ppb per year).
- `rx_r::Float`:        Rate of change for the rotation parameter of the X component (mas per year).
- `ry_r::Float`:        Rate of change for the rotation parameter of the Y component (mas per year).
- `rz_r::Float`:        Rate of change for the rotation parameter of the Z component (mas per year).

# References
- Dawson & Woods, 2010. 
- GDA2020 Technical Manual. 
- GNSS Handbook, 2017, Section 2.3.2.  
- IERS Conventions, 2010, Chapter 4 (https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html). 


"""
function coordinate_transformation(xi, yi, zi, convention="AGRS", tx=0.0, ty=0.0, tz=0.0, d=0.0, rx=0.0, ry=0.0, rz=0.0, t1=0.0, t0=0.0, tx_r=0.0, ty_r=0.0, tz_r=0.0, d_r=0.0, rx_r=0.0, ry_r=0.0, rz_r=0.0)

    # Initial coordinates
    Xi = xi
    Yi = yi
    Zi = zi

    # Elapsed time
    # - t0 = epoch of transformation reference (epoch of output TRF coordinates)
    # - t1 = epoch of input TRF coordinates.
    Δt = (t1 - t0)

    # Transformation parameters and associated units
    # - t = translation (metres)
    # - d = scale (ppb -> parts)
    # - r = rotation (millarcseconds -> degrees -> radians)

    # Parameters
    tx = tx
    ty = ty
    tz = tz
    d  = d * 10^-9
    rx = rx / 3600000 * (π/180) 
    ry = ry / 3600000 * (π/180)
    rz = rz / 3600000 * (π/180)

    # Rates
    tx_r = tx_r
    ty_r = ty_r
    tz_r = tz_r
    d_r  = d_r * 10^-9
    rx_r = rx_r / 3600000 * (π/180)
    ry_r = ry_r / 3600000 * (π/180)
    rz_r = rz_r / 3600000 * (π/180)

    # Transformation parameters with rates applied to the elasped time 
    # - If no time or rates given, 2nd terms are cancelled out by being equal to zero (and so a 7-parameter transformation)
    Tx = tx + tx_r*Δt
    Ty = ty + ty_r*Δt
    Tz = tz + tz_r*Δt
    D  =  d + d_r*Δt
    Rx = rx + rx_r*Δt
    Ry = ry + ry_r*Δt
    Rz = rz + rz_r*Δt

    # Transformation equation
    # - Equation 2.29 (GNSS Handbook):    
    #       
    #       𝐗₂ = 𝐗₁ + 𝚻 + D𝐗₁ + 𝐑𝐗₁
    #           
    #   Where:
    #               𝐗₂ is output TRF coordinate vector (X,Y,Z). 
    #               𝐗₁ is input TRF coordinate vector (X,Y,Z).
    #               𝚻  is translation parameter vector (Tx, Ty, Tz). 
    #               D  is scale parameter (scalar). 
    #               𝐑  is rotation parameter matrix (Rx, Ry, Rz). 
    #                  - Condition for AGRS/IERS sign convention. 
            
    if convention=="AGRS"

        # AGRS sign convention
        XYZ = [Xi Yi Zi]' + [Tx Ty Tz]' + D*[Xi Yi Zi]' + [0 Rz -Ry; -Rz 0 Rx; Ry -Rx 0]*[Xi Yi Zi]'
    
    end

    if convention=="IERS"

        # IERS sign convention
        XYZ = [Xi Yi Zi]' + [Tx Ty Tz]' + D*[Xi Yi Zi]' + [0 -Rz Ry; Rz 0 -Rx; -Ry Rx 0]*[Xi Yi Zi]'

    end
    
    # Return 
    return XYZ[1], XYZ[2], XYZ[3] 

end

"""
    plate_motion_model(xi, yi, zi, convention="AGRS", t1, t2, rx_r=0.00150379, ry_r=0.00118346, rz_r=0.00120716)

Propagate coordinates from one epoch to another on the same reference frame using a plate motion model.
- The default rotation rates are for the Australian Plate Motion Model (APMM). 
- This function requires the rotation matrix sign convention to be set to "AGRS" (default) or "IERS"

# Arguments
- `xi::Float`:              X component of cartesian coordinates (metres).
- `yi::Float`:              Y component of cartesian coordinates (metres).
- `zi::Float`:              Z component of cartesian coordinates (metres).
- `convention::String`:     Sign convention for roation matrix. "AGRS" (default) or "IERS".
- `t1::Float`:              Epoch of input coordinates (decimal years). Epoch to propagate from.
- `t2::Float`:              Epoch of output coordinates (decimal years). Epoch to propagate to.
- `rx_r::Float`:            Rate of change for the rotation parameter of the X component (mas per year).
- `ry_r::Float`:            Rate of change for the rotation parameter of the Y component (mas per year).
- `rz_r::Float`:            Rate of change for the rotation parameter of the Z component (mas per year).

# References
- GDA2020 Technical Manual. 

"""
function plate_motion_model(xi, yi, zi, t1, t2, convention="AGRS", rx_r=1.50379, ry_r=1.18346, rz_r=1.20716)

    # Initial coordinates
    Xi = xi
    Yi = yi
    Zi = zi

    # Propagation time
    # - t1 = epoch of input TRF coordinates.
    # - t2 = epoch to propagate to. 
    Δt = (t1 - t2)

    # Organise units
    # - _r = rotation rates (millarcseconds -> degrees -> radians)
    rx_r = rx_r / 3600000 * (π/180)
    ry_r = ry_r / 3600000 * (π/180)
    rz_r = rz_r / 3600000 * (π/180)

    # Apply rates to propagtion time delta
    Rx = rx_r*Δt
    Ry = ry_r*Δt
    Rz = rz_r*Δt

    # Transformation equation
    # - Equation 2.29 (GNSS Handbook):    
    #       
    #       𝐗₂ = 𝐗₁ + 𝐑𝐗₁
    #           
    #   Where:
    #               𝐑  is rotation parameter matrix (Rx, Ry, Rz). 
    #                  - Condition for AGRS/IERS sign convention. 

    if convention=="AGRS"

        # AGRS sign convention
        XYZ = [Xi Yi Zi]' + [0 Rz -Ry; -Rz 0 Rx; Ry -Rx 0]*[Xi Yi Zi]'
    
    end

    if convention=="IERS"

        # IERS sign convention
        XYZ = [Xi Yi Zi]' + [0 -Rz Ry; Rz 0 -Rx; -Ry Rx 0]*[Xi Yi Zi]'

    end
    
    # Return 
    return XYZ[1], XYZ[2], XYZ[3]

end

"""
    get_velocities_APMM(xi, yi, zi)

Obtain velocity at coordinates based on Australian Plate Motion Model (APMM). 
- The result of this function is a velocity (metres per year) for each component. 

# Arguments
- `xi::Float`:    X component of cartesian coordinates (metres).
- `yi::Float`:    Y component of cartesian coordinates (metres).
- `zi::Float`:    Z component of cartesian coordinates (metres).

# Example
- From the GDA2020 Technical Manal (section 3.3), this example can be used to check that the function is working properly.  
- Correct results are below in velocity (metres per year) for each component.

## Execute function
`get_velocities_APMM(-4052052.6588, 4212835.9938, -2545104.6946)`

## Results
X=-0.03925825664338662, Y=-0.005159255804548842, Z= 0.053962976434117876

"""
function get_velocities_APMM(xi, yi, zi)

    # Initial coordinates
    x = xi
    y = yi
    z = zi

    # Rotation rate (arcseconds per year)
    rot_rate = [0.00150379, 0.00118346, 0.00120716]

    # Rotation rate (radians per year) 
    ω = rot_rate*(π/(180*3600))

    # Rotation matrix
    R = [0 z -y; -z 0 x; y -x 0]

    # Velocity (metres per year)
    V = R*ω

    # Assign componants
    vX = V[1]
    vY = V[2]
    vZ = V[3]

    # Return
    return vX, vY, vZ

end

"""
    xyz2llh(x, y, z, a, e2)

Convert global cartesian coordinates to geodetic coordinates.
- X, Y and Z (cartesian), to longitude, latitude and ellipsoid height (geodetic). 
- Requires semi-major axis (a), and first eccentricity (e2) values of reference ellipsoid. 
- GRS80 values will be used as default if ellipsoidal parameters are not given. 
- Longitude, latitude are returned as decimal degrees. Ellipsoid height is in metres.
- Return three seperate components (i.e. lon, lat, h).
- Reference: Teunissen & Montenbruck 2017, Section 2.2, p. 32.

# Arguments
- `x::Float`:    X component of Earth-centred Earth-fixed (ECEF) cartesian system (m).
- `y::Float`:    Y component of Earth-centred Earth-fixed (ECEF) cartesian system (m).
- `z::Float`:    Z component of Earth-centred Earth-fixed (ECEF) cartesian system (m).
- `a::Float`:    Semi-major axis of reference ellipsoid (m).
- `e2::Float`:   First eccentricity of reference ellipsoid.

## Execute function
`xyz2llh(-4052052.623236821, 4212836.107268005, -2545104.5753482883)`

## Results
Lon=133.88552, Lat=-23.67011, H=603.2478

"""
function xyz2llh(x, y, z, a=6378137, e2=0.00669438002290)

    # Input cartesian coordinates (X, Y, Z)
    X = x
    Y = y
    Z = z
        
    # Input reference ellipsoid parameters (default is GRS80)
    # - a: semi-major axis (m)
    # - e2: first eccentricity
    a = a
    e2 = e2

    # Longitude (Equation 2.15, Teunissen & Montenbruck 2017, Section 2.2, p. 31)
    λ = atan(Y, X)*180/π

    # Latitude (ϕ) and ellipsoid height (h), computed iteratively 
    # - initial value of latitude needed (ϕ₀), and height (h=0)
    # - equation 2.21, Teunissen & Montenbruck 2017, Section 2.2, p. 32.
    ϕ₀ = atan(Z/sqrt(x^2 + y^2) * (1 + e2/(1 - e2)))
    h = 0

    for i=1:100

        # Radius of curvature of the prime vertical (Equation 2.19, Teunissen & Montenbruck 2017, Section 2,2, p. 32)
        N = a/(sqrt(1 - e2*sin(ϕ₀)^2))
        
        # Height (Equation 2.22, Teunissen & Montenbruck 2017)
        h = sqrt(X^2 + Y^2)*cos(ϕ₀) + Z*sin(ϕ₀) - a*sqrt(1 - e2*(sin(ϕ₀)^2))

        # Latitude (Equation 2.20, Teunissen & Montenbruck 2017)
        ϕ = atan(Z/sqrt(x^2 + y^2) * (1 + e2*N*sin(ϕ₀)/Z))

        # Iteration
        ϕ₀ = ϕ

    end

    # Latitude
    ϕ = ϕ₀ * 180/π

    # Return 
    return λ, ϕ, h

end

"""
    xyz2enu(Δx, Δy, Δz, λ, ϕ)

Convert global geocentric cartesian coordinate vector to local right-handed cartesian coordinate system.
- X, Y, Z (global), to E, N, U (local).
- Return three seperate components (i.e. ΔE, ΔN, ΔU).
- Coordinates of point Q with respect to local origin point P.
- Reference: Teunissen & Montenbruck 2017, Section 2.2, p. 33.  

# Arguments
- `Δx::Float`:    delta X component, coordinate difference (metres).
- `Δy::Float`:    delta Y component, coordinate difference (metres).
- `Δz::Float`:    delta Z component, coordinate difference (metres).
- `λ::Float`:     longitude of local origin point P (decimal degrees).
- `ϕ::Float`:     latitude of local origin point P (decimal degrees).

## Execute function
`xyz2enu(-4372.27392023982, 7382.3983550289, 9211.8395301752, 145.3521, -30.2912)`

## Results
E=-3587.437523474417, N=11885.443961169742, U=2083.5273208438184

"""
function xyz2enu(Δx, Δy, Δz, λ, ϕ)

    # Inputs
    ΔX = Δx
    ΔY = Δy
    ΔZ = Δz
    λ = λ * π/180
    ϕ = ϕ * π/180

    # Vector of coordinates
    X = [ΔX ΔY ΔZ]'

    # Rotation matrix (global --> local)
    R = [-sin(λ) cos(λ) 0; 
         -sin(ϕ)*cos(λ) -sin(ϕ)*sin(λ) cos(ϕ); 
          cos(ϕ)*cos(λ)  cos(ϕ)*sin(λ) sin(ϕ)]

    # Coordinate transformation (global --> local)
    ENU = R*X

    # Extract components 
    ΔE = ENU[1]
    ΔN = ENU[2]
    ΔU = ENU[3]

    # Return
    return ΔE, ΔN, ΔU

end

"""
    xyz2enu_vcv(xx, yy, zz, xy, xz, yz, λ, ϕ)

Convert variance-covariance of global geocentric cartesian coordinates to local cartesian right-handed coordinate system.
- XX, YY, ZZ,  (global), to EE, NN, UU (local).
- Return individual VCV parameters.
- Variance and covariance for coordinates of point Q with respect to local origin point P.
- Reference: Teunissen & Montenbruck 2017, Section 2.2, p. 33. 

# Arguments
- `xx::Float`:    variance of X component (metres).
- `yy::Float`:    variance of Y component (metres).
- `zz::Float`:    variance of Y component (metres).
- `xy::Float`:    covariance of X and Y components (metres).
- `xz::Float`:    covariance of X and Z components (metres).
- `xy::Float`:    covariance of Y and Z components (metres).
- `λ::Float`:     longitude of local origin point P (decimal degrees).
- `ϕ::Float`:     latitude of local origin point P (decimal degrees).

## Execute function
`xyz2enu_vcv(0.000038292676, 0.00014426068, 0.000041634281, -0.000013097656, 0.000002896884, -0.000023503130, 96.8339, -12.1883)`

## Results
EE=0.000036698200, NN=0.000036505950, UU=0.000150983489, EN=-0.000000034120, EU=0.000000219055, NU=-0.000000062007

"""
function xyz2enu_vcv(xx, yy, zz, xy, xz, yz, λ, ϕ)

    # Assign parameters
    xx = xx
    yy = yy
    zz = zz
    xy = xy
    xz = xz
    yz = yz
    λ = λ * π/180
    ϕ = ϕ * π/180

    # VCV matrix
    Q = [xx xy xz; xy yy yz; xz yz zz]

    # Rotation matrix (global --> local)
    R = [-sin(λ) cos(λ) 0; 
         -sin(ϕ)*cos(λ) -sin(ϕ)*sin(λ) cos(ϕ); 
          cos(ϕ)*cos(λ) cos(ϕ)*sin(λ) sin(ϕ)]

    # Compute transformation (global --> local)
    vcv = R*Q*R'

    # Return
    return vcv[1,1], vcv[2,2], vcv[3,3], vcv[1,2], vcv[1,3], vcv[2,3]

end

"""
    llh2xyz(λ, ϕ, h, a, e2)

Convert geodetic coordinates to cartesian coordinates. 
- Longitude, latitude and ellipsoidal height (geodetic), to X, Y and Z (cartesian). 
- Requires semi-major axis (a), and first eccentricity (e2) values of reference ellipsoid. 
- GRS80 values will be used as default if ellipsoidal parameters are not given. 
- Return XYZ coordinates in three seperate components (i.e. X, Y, Z).
- Reference: Teunissen & Montenbruck 2017, Section 2.2, p. 32.

# Arguments
- `λ::Float`:    Geodetic Longitude (decimal degrees).
- `ϕ::Float`:    Geodetic Latitude (decimal degrees).
- `h::Float`:      Ellipsoid height (m).
- `a::Float`:      Semi-major axis of reference ellipsoid (m).
- `e2::Float`:     First eccentricity of reference ellipsoid

## Execute function
`llh2xyz(133.88552, -23.67011, 603.2478)`

## Results
X=-4052052.623236821, Y=4212836.107268005, Z=-2545104.5753482883
 
"""
function llh2xyz(λ, ϕ, h, a=6378137, e2=0.00669438002290)

    # Input geodetic coordinates (longitude, latitude, ellipsoid height)
    λ = λ * π/180
    ϕ = ϕ * π/180
    h = h

    # Input ellipsoid parameters (default is GRS80)
    # - a: semi-major axis (m)
    # - e2: first eccentricity 
    a = a
    e2 = e2

    # Radius of curvature of the prime vertical (Equation 2.19, Teunissen & Montenbruck 2017)
    N = a/(sqrt(1 - e2*sin(ϕ)^2))

    # Compute geodetic to cartesian (Equation set 2.18, Teunissen & Montenbruck 2017)
    X = (N + h)*cos(ϕ)*cos(λ)
    Y = (N + h)*cos(ϕ)*sin(λ)
    Z = ((1 - e2)*N + h)*sin(ϕ)

    # Return
    return X, Y, Z

end

"""
    enu2xyz(Δe, Δn, Δu, λ, ϕ)

Convert local right-handed cartesian coordinates to global cartesian coordinate system.
- E, N, U (local), to X, Y, Z (global).
- Return three seperate components (i.e. ΔX, ΔY, ΔZ).
- Coordinates of point Q with respect to local origin point P.
- Reference: Teunissen & Montenbruck 2017, Section 2.2, p. 33.

# Arguments
- `Δe::Float`:    delta E component, coordinate difference (metres).
- `Δn::Float`:    delta N component, coordinate difference (metres).
- `Δu::Float`:    delta U component, coordinate difference (metres).
- `λ::Float`:     longitude of local origin point P (decimal degrees).
- `ϕ::Float`:     latitude of local origin point P (decimal degrees).

## Execute function
`enu2xyz(-3587.437523474417, 11885.443961169742, 2083.5273208438184, 145.3521, -30.2912)`

## Results
X=-4372.27392023982, Y=7382.3983550289, Z=9211.8395301752

"""
function enu2xyz(Δe, Δn, Δu, λ, ϕ)

    # Inputs
    ΔE = Δe
    ΔN = Δn
    ΔU = Δu
    λ = λ*π/180
    ϕ = ϕ*π/180

    # Vector of coordinates
    X = [ΔE ΔN ΔU]'

    # Rotation matrix (local --> global)   
    R = [-sin(λ) -sin(ϕ)*cos(λ) cos(ϕ)*cos(λ);
          cos(λ) -sin(ϕ)*sin(λ) cos(ϕ)*sin(λ);
          0 cos(ϕ) sin(ϕ)]

    # Compute transformation (local --> global)
    XYZ = R*X

    # Extract components
    ΔX = XYZ[1]
    ΔY = XYZ[2]
    ΔZ = XYZ[3]

    # Return
    return ΔX, ΔY, ΔZ

end


"""
    enu2xyz_vcv(ee, nn, uu, en, eu, nu, λ, ϕ)

Convert variance-covariance of local right-handed cartesian coordinates to global cartesian coordinate system.
- EE, NN, UU,  (local), to XX, YY, ZZ (global).
- Return a single vcv matrix (i.e. [xx xy xz; xy yy yz; xz yz zz]).
- Variance and covariance for coordinates of point Q with respect to local origin point P.
- Reference: Teunissen & Montenbruck 2017, Section 2.2, p. 33. 

# Arguments
- `ee::Float`:    variance of east component (metres).
- `nn::Float`:    variance of north component (metres).
- `uu::Float`:    variance of up component (metres).
- `en::Float`:    covariance of east and north components (metres).
- `eu::Float`:    covariance of east and up components (metres).
- `nu::Float`:    covariance of nort and u components (metres).
- `lon::Float`:   longitude of local origin point P (decimal degrees).
- `lat::Float`:   latitude of local origin point P (decimal degrees).

## Execute function
`enu2xyz_vcv(0.000036698200, 0.000036505950, 0.000150983489, -0.000000034120, 0.000000219055, -0.000000062007, 96.8339, -12.1883)`

## Results
XX=0.000038292676, YY=0.00014426068, ZZ=0.000041634281, XY=-0.000013097656, XZ=0.000002896884, YZ=-0.000023503130
 
"""
function enu2xyz_vcv(ee, nn, uu, en, eu, nu, λ, ϕ)

    # Assign parameters
    ee = ee
    nn = nn
    uu = uu
    en = en
    eu = eu
    nu = nu
    λ = λ*π/180
    ϕ = ϕ*π/180

    # VCV matrix
    Q = [ee en eu; en nn nu; eu nu uu]

    # Rotation matrix (local --> global)   
    R = [-sin(λ) -sin(ϕ)*cos(λ) cos(ϕ)*cos(λ);
          cos(λ) -sin(ϕ)*sin(λ) cos(ϕ)*sin(λ);
          0 cos(ϕ) sin(ϕ)]

    # Compute transformation (local --> global)
    xyzVCV = R*Q*R'

    # Return
    return xyzVCV

end

"""
    enu2los(α, θ, e, n, u)

Project east, north, up (ENU) displacement vector to line-of-sight (LoS) using incidence and azimuth angles of satellite.
- LoS look vector has origin at target (positive value is towards satellite).

# Arguments
- `α::Float`:     Azimuth angle of target-to-satellite (radians).
- `θ ::Float`:    Incidence angle of satellite (radians).
- `e::Float`:     East componant of displacement vector (metres).
- `n::Float`:     East componant of displacement vector (metres).
- `u::Float`:     East componant of displacement vector (metres).

# References
Equation 2.6, An analysis of the InSAR displacement vector decomposition (Wietske Brouwer, 2021).

    dₗₒₛ =  dₑSin(θ)Sin(α) + dₙSin(θ)Cos(α) + dᵤCos(θ)

    Where:     dₗₒₛ is line-of-sight displacement vector.
               dₑ, dₙ, dᵤ are the east, north and up discplacement.
               θ is the incidence angle between ellipsoid normal and line-of-sight vector.
               α is azimuth angle of line-of-sight vector from target to satellite. 

- See section 2.2 of Brouwer (2021) for other viewing geometry conventions. 
- Also, see equation 5.1.1 of Hansson (2001) for another perspective. 

"""
function enu2los(α, θ, e, n, u)

    # ENU to LOS projection
    # - Equation 2.6/2.7, Brouwer, (2021)
    # - Multiplication of line-of-sight projector and discplacement vector.  
    LOS = [sin(θ)sin(α) sin(θ)cos(α) cos(θ)]*[e n u]'
    #LOS = e*sin(θ)sin(α) + n*sin(θ)cos(α) + u*cos(θ)
    
    # Return
    return LOS[1]

end

"""
    enu2los_vcv(α, θ, ee, nn, uu, en, eu, nu)

Project east, north, up uncertainties to line-of-sight (LoS) vector.
- See the function enu2los() for more details. 

# Arguments
- `α::Float`:    Azimuth angle of target-to-satellite (radians).
- `θ::Float`:    Incidence angle of satellite (radians).
- `ee::Float`:   Variance for east componant of displacement vector (metres).
- `nn::Float`:   Variance for north componant of displacement vector (metres).
- `uu::Float`:   Variance for up componant of displacement vector (metres).
- `en::Float`:   Covariance for east componant of displacement vector (metres).
- `eu::Float`:   Covariance for north componant of displacement vector (metres).
- `nu::Float`:   Covariance for up componant of displacement vector (metres).

"""
function enu2los_vcv(α, θ, ee, nn, uu, en, eu, nu)

    # VCV matrix
    Q = [ee en eu; en nn nu; eu nu uu]

    # Unit vector for projector (ENU --> LoS)
    A = [sin(θ)sin(α) sin(θ)cos(α) cos(θ)]
    
    # Propagation of uncertainties
    LOS_VCV = A*Q*A'

    # Return
    return LOS_VCV[1]

end

"""
    doy2date(yyyy, doy)

Retreive calendar date (Year, Month, Day), from given year and day-of-year (DOY). 
- Operates by converting to modified Julian date (MJD), then adding the DOY, before then converting to calendar date.
- See this link for explaination: http://www.mat.uc.pt/~efemast/netscape/help/en/tem_dj.htm.
- Returns DateTime object.

# Arguments
- `yyyy::Integer`:    Year (e.g. 2021).
- `doy::Integer`:    Day of year (e.g. 190)

"""
function doy2date(yyyy, doy)

    # --- Clarify year variable
    if length(string(yyyy))==4
        yyyy = yyyy
    elseif length(string(yyyy))<=2
        if yyyy < 68
            yyyy = 2000 + yyyy
        elseif yyyy > 68
            yyyy = 1900 + yyyy
        end
    end

    ## --- Convert from YYYY-DOY to Modified Julian Date (MJD)

    # Initialise day (dd) and month (mm) for given year
    dd = 0
    mm = 1

    # Account for Julian year increment being at 14 months
    if mm > 2
        y = yyyy
        m = mm
    else
        y = yyyy - 1
        m = mm + 12
    end

    # Set date conversion parameters for Julian date equation
    Y = y
    M = m
    D = dd
    A = floor(Y/100)
    B = 2 - A + floor(A/4)
    E = floor(365.25*(Y+4716))
    F = floor(30.6001*(M+1))

    # Equation for Julian date (JD)
    JD = B + D + E + F - 1524.5

    # Modified Julian date (MJD) 
    MJD = JD - 2400000.5 + doy

    ## --- Convert from MJD to Calendar date (YYY, MM, DD)

    # Julian date with day-of-year added
    JD = MJD + 2400000.5

    # Set date conversion parameters
    Z = floor(JD + 0.5)
    F = JD - Z

    # Account for calendar adjustment (4th of October 1582)
    if Z < 2299161
        A = Z
    else
        G = floor((Z-1867216.25)/36524.25)
        A = Z + 1 + G - floor(G/4)
    end

    B = A + 1524
    C = floor((B-122.1)/365.25)
    D = floor(365.25*C)
    E = floor((B-D)/30.6001)

    # Day (DD), without fraction of a day. If fraction of day needed, do: + F
    DD = B - D - floor(30.6001*E)

    # Account for Julian year increment being at 14 months
    # Month
    if E < 14
        MM = E - 1
    else
        MM = E - 13
    end
    # Year
    if MM > 2
        YYYY = C - 4716
    else
        YYYY = C - 4715
    end

    # Return
    return Date(Int(YYYY), Int(MM), Int(DD))

end

"""
   dd2ddm(dd)

Convert coordinate in decimal-degrees (dd) format to degrees-decimal-minutes format (ddm).

# Arguments
- `dd::Float`:	decimal degrees.

"""
function dd2ddm(dd)
	
	# Degree calculation
	d = floor(abs(dd))

	# Minute calculation
	m = floor((abs(dd) - d) * 60) 

	# Second calculation
	s = (abs(dd) - d - m/60) * 3600 
	
	# Decimal minute calculation
	dm = m + s/60

	if dd >= 0

		# Return
		return Int(d), dm 
	else

		# Return 
	    return Int(-d), dm
	end
end

"""
   dd2dms(dd)

Convert coordinate in decimal-degrees (dd) format to degrees-minute-seconds format (dms).

# Arguments
- `dd::Float`:	decimal degrees.

"""
function dd2dms(dd)
	
	# Degree calculation
	d = floor(abs(dd))

	# Minute calculation
	m = floor((abs(dd) - d) * 60) 

	# Second calculation
	s = (abs(dd) - d - m/60) * 3600 

	# Account for sign
	if dd >= 0
		return Int(d), Int(m), s
	else
	    return Int(-d), Int(m), s
	end
end

"""
   dms2dd(d, m, s)

Convert coordinate in degrees-minute-seconds (dms) format to decimal degrees format (dd).

# Arguments
- `d::Int`:	    degrees.
- `m::Int`:     mintues.
- `s::Float`:   seconds.

"""
function dms2dd(d, m, s)

	# Calculation
	dd = abs(d) + m/60 + s/3600

	# Account for sign
	if d >= 0
		return dd
	else
		return -dd
	end
end

"""
   vincentys_inverse(λ₁,ϕ₁,λ₂,ϕ₂,a,f)

Given the geodetic coordinates of two points, calculate the ellipsoidal arc distance and forward/reverse geodetic azimuths.
- Equation numbering from AGRS Compendium.

# Arguments
- `λ₁::Float`:	Longitude of point 1 (decimal degrees).
- `ϕ₁::Float`:	Latitude of point 1 (decimal degrees).
- `λ₂::Float`:	Longitude of point 2 (decimal degrees).
- `ϕ₂::Float`:	Latitude of point 2 (decimal degrees).
- `a::Float`:   Semi-major axis of the ellipsoid.
- `f::Float`:   Flattening of the ellipsoid.

# References
- AGRS Compendium, Section 12.1.3, https://www.icsm.gov.au/publications/australian-geospatial-reference-system-compendium.
- GeodePy, https://github.com/GeoscienceAustralia/GeodePy/tree/master.
- Vincenty 1975, https://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf.
- Deakin 2010, http://www.mygeodesy.id.au/documents/Trav_Comp_V2.1.pdf.
"""
function vincentys_inverse(λ₁,ϕ₁,λ₂,ϕ₂,a,f)

    # Derived quantities
    # - b, Semi-minor axis of ellipsoid.
    b = a*(1-f)

    # Calculate reduced latitude
    # - Eq. 71 and 72.
    u1 = atan((1 - f) * tan(ϕ₁))
    u2 = atan((1 - f) * tan(ϕ₂))
    
    # Calculate initial approximation of longitude difference
    # - Difference in geodetic longitude on a sphere (λ)
    # - Eq. 72.
    λ = λ₂ - λ₁
    ω = λ

    # Iteratively estimate the longitude difference
    # - Iterate until no significant change (< 1e-12), 
    #   or at 1000 iterations. 

    # Set initial values
    α = 0
    σ = 0
    cos2σ𝑚 = 0

    # Iterate
    for i in 1:1000

        # Eq. 74
        sinσ = sqrt((cos(u2)sin(λ))^2 + (cos(u1)sin(u2) - sin(u1)cos(u2)cos(λ))^2)

        # Eq. 75
        cosσ = sin(u1)sin(u2) + cos(u1)cos(u2)cos(λ)

        # Eq. 76
        σ = atan(sinσ/cosσ)

        # Eq. 77
        α = asin(cos(u1)*cos(u2)*sin(λ)/sin(σ))

        # Eq. 78
        cos2σ𝑚 = cosσ - (2sin(u1)sin(u2)/cos(α)^2)

        # Eq. 79
        C = (f/16)cos(α)^2 * (4 + f*(4 - 3cos(α)^2))

        # Eq. 80
        λ̂ = ω + (1 - C) * f * sin(α) * (σ + C * sin(σ) * (cos2σ𝑚 + C * cos(σ) * (-1 + 2cos2σ𝑚^2)))

        # Update
        Δλ = λ̂ - λ 
        λ = λ̂
        
        # Stop condition
        if abs(Δλ) < 1e-12
            break
        end

    end
    
    # Eq. 81
    u² = cos(α)^2 * (a^2 - b^2)/b^2

    # Eq. 82
    A = 1 + (u²/16384) * (4096 + u² * (-768 + u² * (320 - 175u²)))

    # Eq. 83
    B = (u²/1024) * (256 + u² * (-128 + u² * (74 - 47u²)))
 
    # Eq. 84
    Δσ = B*sin(σ) * (cos2σ𝑚 + (B/4) * (cos(σ) * (-1 + 2cos2σ𝑚^2) - (B/6) * cos2σ𝑚 * (-3 + 4sin(σ)^2) * (-3 + 4cos2σ𝑚^2)))
    
    # === 
    # FINAL CALCULATIONS
    # - Geodesic, the elliposidal arc distance (Eq. 85)
    # - Geodetic azimuth (Eq. 86 and 87)
    
    # Eq. 85
    # - Geodesic (s)
    s = b*A*(σ - Δσ)

    # Eq. 86 and 87
    # - Point 1 to Point 2 (α₁₋₂).
    # - Point 2 to point 1 (α₂₋₁).
    α₁₋₂ = atan((cos(u2)sin(λ)) / (cos(u1)sin(u2) - sin(u1)cos(u2)cos(λ)))
    α₂₋₁ = atan((cos(u1)sin(λ)) / (-sin(u1)cos(u2) + cos(u1)sin(u2)cos(λ)))

    # Radians to degrees
    α₁₋₂ = α₁₋₂ * 180/π
    α₂₋₁ = α₂₋₁ * 180/π

    # Sign condition
    if α₁₋₂ < 0
        α₁₋₂ = α₁₋₂ + 360
    end

    # Reverse azimuth
    α₂₋₁ = α₂₋₁ + 180

    # Return
    return s, α₁₋₂, α₂₋₁

end
