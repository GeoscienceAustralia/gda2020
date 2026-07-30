## ===
# DESCRIPTION
# This file contains some useful functions to handle DynAdjust output. 
#
# Functions:    - read_DynAdjustAdjustmentSummary_BLOCK()
#               - dynxyz2dataframe()
#               - dynadj2dataframe()
#               - dyncor2dataframe()
#               - dynm2s2dataframe()
#               - dyndst2dataframe()
#               - dynrename2dataframe()
#               - dynnear2dataframe()
#               - dynstnxml2dataframe() | In Development.

"""
    read_DynAdjustAdjustmentSummary_BLOCK(fp)

"""
function read_DynAdjustAdjustmentSummary_BLOCK(fp)

    # Read block
    lines = readlines(fp)
    block_start = findfirst(contains("ITERATION"), lines)[1]
    block_end  = findfirst(contains("Chi-Square test"), lines)[1]
    adjustment_summary_block = lines[block_start:block_end]
    dash_lines_index = findall(contains("--------------------------------------------------------------------------------"), adjustment_summary_block)
    deleteat!(adjustment_summary_block, dash_lines_index)

    # Return
    return adjustment_summary_block

end

"""
    dynxyz2dataframe(fp)

Read in XYZ coordinates from DynAdjust file (.xyz). 

# Arguments
`f::String`:    File path to .xyz DynAdjust file.

"""
function dynxyz2dataframe(fp)

    # Read lines
    lines = readlines(fp)
    lines = deleteat!(lines, 1:20)
    lines = filter(!isempty, lines)

    # Extract columns (single-line block)
    site = SubString.(lines, 1, 20)
    constraint = SubString.(lines, 21, 23)
    easting = SubString.(lines, 27, 39)
    northing = SubString.(lines, 41, 55)
    zone = SubString.(lines, 59, 62)
    lon = SubString.(lines, 77, 91)
    lat = SubString.(lines, 63, 76)
    H = SubString.(lines, 92, 102)
    h = SubString.(lines, 103, 113)
    x = SubString.(lines, 114, 128)
    y = SubString.(lines, 129, 143)
    z = SubString.(lines, 144, 158)
    sigmaE = SubString.(lines, 159, 170)
    sigmaN = SubString.(lines, 171, 180)
    sigmaU = SubString.(lines, 181, 190)

    # Organise format
    site = String.(rstrip.(site))
    constraint = String.(strip.(constraint))
    easting = parse.(Float64, easting)
    northing = parse.(Float64, northing)
    zone = parse.(Int, zone)
    lon = parse.(Float64, lon)
    lat = parse.(Float64, lat)
    H = parse.(Float64, H)
    h = parse.(Float64, h)
    x = parse.(Float64, x)
    y = parse.(Float64, y)
    z = parse.(Float64, z)
    sigmaE = parse.(Float64, sigmaE)
    sigmaN = parse.(Float64, sigmaN)
    sigmaU = parse.(Float64, sigmaU)

    # Form dataframe
    df = DataFrame(
                    soln=site,
                    site=site,
                    cons=constraint,
                    lon=lon,
                    lat=lat,
                    H=H,
                    h=h,
                    x=x, 
                    y=y, 
                    z=z, 
                    e=easting,
                    n=northing,
                    zone=zone,
                    sigmaE=sigmaE, 
                    sigmaN=sigmaN,
                    sigmaU=sigmaU,
    )
    sort!(df)
    
    # Return
    return df

end

"""
    dynadj2dataframe(fp)

- Currently forming dataframe of adjusted measurements only. 
- ToDo: Add seperate dataframe of adjusted coordinates when needed.
- ToDo: Refine handling of measurement types.

"""
function dynadj2dataframe(fp)

    # Obtain adjusted measurements block
    lines = readlines(fp)
    block_start = findfirst(contains("Adjusted Measurements"), lines)[1]
    block_end = findfirst(contains("Adjusted Coordinates"), lines)[1]
    adjusted_msr_block = lines[block_start+5:block_end-3]
    lines = nothing

    # Seperate measurement types that occur over multiple lines (D)
    # - Single line cartesian measurement types (S, E, C, M, R, H, L, Y, X)
    # - Single line angle measurement types (P, Q, B, A, Z, I, J, K, V)
    # - Multiple line measurement types (D)
    # - ToDo: refine the list here.
    # - ToDo: convert DMS to DD for convenience (currently DMS are strings)

    # Single line cartesian measurements
    block = filter(
                    r->SubString.(r, 1, 1)!="D"  
                    && SubString.(r, 1, 1)!=" "
                    && SubString.(r, 1, 1)!="P"
                    && SubString.(r, 1, 1)!="Q"
                    && SubString.(r, 1, 1)!="B"
                    && SubString.(r, 1, 1)!="A"
                    && SubString.(r, 1, 1)!="Z"
                    && SubString.(r, 1, 1)!="I"
                    && SubString.(r, 1, 1)!="J"
                    && SubString.(r, 1, 1)!="K"
                    && SubString.(r, 1, 1)!="V", 
                    adjusted_msr_block
    )
    # Single line angle measurement types
    block_angles = filter(
                            r->SubString.(r, 1, 1)!="D"  
                            && SubString.(r, 1, 1)!=" "
                            && SubString.(r, 1, 1)!="S"
                            && SubString.(r, 1, 1)!="E"
                            && SubString.(r, 1, 1)!="C"
                            && SubString.(r, 1, 1)!="G"
                            && SubString.(r, 1, 1)!="M"
                            && SubString.(r, 1, 1)!="R"
                            && SubString.(r, 1, 1)!="H"
                            && SubString.(r, 1, 1)!="L"
                            && SubString.(r, 1, 1)!="Y"
                            && SubString.(r, 1, 1)!="X", 
                            adjusted_msr_block
    )
    # Multi-line measurement types
    block_multi = filter(
                        r->SubString.(r, 1, 1)=="D" 
                        || SubString.(r, 1, 1)==" ", 
                        adjusted_msr_block
    )
    adjusted_msr_block = nothing

    # Extract columns (single-line block)
    type = SubString.(block, 1, 1)
    station1 = SubString.(block, 3, 22)
    station2 = SubString.(block, 23, 42)
    station3 = SubString.(block, 43, 63)
    par = SubString.(block, 66, 66)
    msr = SubString.(block, 68, 86)
    adj = SubString.(block, 88, 105)
    cor = SubString.(block, 107, 117)
    sigmaMSR = SubString.(block, 119, 130)
    sigmaADJ = SubString.(block, 132, 143)
    sigmaCOR = SubString.(block, 145, 156)
    nstat = SubString.(block, 158, 167)
    pelzerRel = SubString.(block, 169, 179) 
    preAdjCor = SubString.(block, 181, 193) 

    # Organise format (single-line block)
    type = String.(type)
    station1 = String.(rstrip.(station1))
    station2 = String.(rstrip.(station2))
    station3 = String.(rstrip.(station3))
    par = String.(par)
    cor = parse.(Float64, cor)
    adj = parse.(Float64, adj)
    msr = parse.(Float64, msr)
    sigmaMSR = parse.(Float64, sigmaMSR)
    sigmaADJ = parse.(Float64, sigmaADJ)
    sigmaCOR = parse.(Float64, sigmaCOR)
    nstat = parse.(Float64, nstat)
    pelzerRel = parse.(Float64, pelzerRel)
    preAdjCor = parse.(Float64, preAdjCor)

    # Setup temporary dataframe (single-line block)
    df = DataFrame(
                    type=type,
                    station1=station1,
                    station2=station2,
                    station3=station3,
                    par=par,
                    msr=msr,
                    adj=adj,
                    cor=cor,
                    sigmaMSR=sigmaMSR, 
                    sigmaADJ=sigmaADJ, 
                    sigmaCOR=sigmaCOR, 
                    nstat=nstat, 
                    pelzerRel=pelzerRel, 
                    preAdjCor=preAdjCor, 
    )

    if isempty(block_angles)==false

        # Extract columns (single-line block angles)
        type = SubString.(block_angles, 1, 1)
        station1 = SubString.(block_angles, 3, 22)
        station2 = SubString.(block_angles, 23, 42)
        station3 = SubString.(block_angles, 43, 63)
        msr = SubString.(block_angles, 68, 86)
        adj = SubString.(block_angles, 88, 105)
        cor = SubString.(block_angles, 107, 117)
        sigmaMSR = SubString.(block_angles, 119, 130)
        sigmaADJ = SubString.(block_angles, 132, 143)
        sigmaCOR = SubString.(block_angles, 145, 156)
        nstat = SubString.(block_angles, 158, 167)
        pelzerRel = SubString.(block_angles, 169, 179) 
        preAdjCor = SubString.(block_angles, 181, 193) 

        # Organise format (single-line block angles)
        type = String.(type)
        station1 = String.(rstrip.(station1))
        station2 = String.(rstrip.(station2))
        station3 = String.(rstrip.(station3))
        par = fill("N/A", length(block_angles))
        if length(split(msr[2]))==3
            msr = split.(msr)
            msr_dd = []
            for i in eachindex(msr)
                push!(msr_dd, dms2dd(parse(Int, msr[i][1]), parse(Int, msr[i][2]), parse(Float64, msr[i][3])))
            end
            msr = msr_dd
        else
            msr = parse.(Float64, strip.(msr))
        end
        if length(split(adj[2]))==3
            adj = split.(adj)
            adj_dd = []
            for i in eachindex(adj)
                push!(adj_dd, dms2dd(parse(Int, adj[i][1]), parse(Int, adj[i][2]), parse(Float64, adj[i][3])))
            end
            adj = adj_dd
        else
            adj = parse.(Float64, strip.(adj))
        end
        cor = parse.(Float64, cor)
        sigmaMSR = parse.(Float64, sigmaMSR)
        sigmaADJ = parse.(Float64, sigmaADJ)
        sigmaCOR = parse.(Float64, sigmaCOR)
        nstat = parse.(Float64, nstat)
        pelzerRel = parse.(Float64, pelzerRel)
        preAdjCor = parse.(Float64, preAdjCor)

        # Setup temporary dataframe (single-line block angles)
        df1 = DataFrame(
                        type=type,
                        station1=station1,
                        station2=station2,
                        station3=station3,
                        par=par,
                        msr=msr,
                        adj=adj,
                        cor=cor,
                        sigmaMSR=sigmaMSR, 
                        sigmaADJ=sigmaADJ, 
                        sigmaCOR=sigmaCOR, 
                        nstat=nstat, 
                        pelzerRel=pelzerRel, 
                        preAdjCor=preAdjCor, 
        )

        # Combine into one dataframe
        df = vcat(df, df1)

    end


    if isempty(block_multi)==false

        # Extract columns (multi-line block)
        type = []
        station1 = []
        station2 = []
        station3 = []
        par = []
        msr = []
        adj = []
        cor = []
        sigmaMSR = []
        sigmaADJ = []
        sigmaCOR = []
        nstat = []
        pelzerRel = []
        preAdjCor = []
        for l in eachindex(block_multi)

            line = block_multi[l]
            if length(line)==67
            
                # Extract column
                type_l = SubString(line, 1, 1)
                station1_l = SubString(line, 3, 22)
                station2_l = SubString(line, 23, 42)
                par_l =  SubString(line, 66, 66)

                # Oragnise format
                type_l = String(type_l)
                station1_l = String(rstrip(station1_l))
                station2_l = String(rstrip(station2_l))

                # Add to vector
                push!(type, type_l)
                push!(station1, station1_l)
                push!(station2, station2_l)
                push!(station3, "-")
                push!(par, par_l)
                push!(msr, NaN)
                push!(adj, NaN)
                push!(cor, NaN)
                push!(sigmaMSR, NaN)
                push!(sigmaADJ, NaN)
                push!(sigmaCOR, NaN)
                push!(nstat, NaN)
                push!(pelzerRel, NaN)
                push!(preAdjCor, NaN)

            else
            
                # Extract Column from line
                station3_l = SubString(line, 43, 63)
                msr_l = SubString(line, 68, 86)
                adj_l = SubString(line, 88, 105)
                cor_l = SubString(line, 107, 117)
                sigmaMSR_l = SubString(line, 119, 130)
                sigmaADJ_l = SubString(line, 132, 143)
                sigmaCOR_l = SubString(line, 145, 156)
                nstat_l = SubString(line, 158, 167)
                pelzerRel_l = SubString(line, 169, 179)
                preAdjCor_l = SubString(line, 181, 193)

                # Oragnise format
                station3_l = String(rstrip(station3_l))
                if length(split(msr_l))==3
                    msr_l = split(msr_l)      
                    msr_l = dms2dd(parse(Int, msr_l[1]), parse(Int, msr_l[2]), parse(Float64, msr_l[3]))
                else
                    msr_l = parse(Float64, msr_l)
                end
                if length(split(adj_l))==3
                    adj_l = split(adj_l)      
                    adj_l = dms2dd(parse(Int, adj_l[1]), parse(Int, adj_l[2]), parse(Float64, adj_l[3]))
                else
                    adj_l = parse(Float64, adj_l)
                end
                cor_l = parse(Float64, cor_l)
                sigmaMSR_l = parse(Float64, sigmaMSR_l)
                sigmaADJ_l = parse(Float64, sigmaADJ_l)
                sigmaCOR_l = parse(Float64, sigmaCOR_l)
                nstat_l = parse(Float64, nstat_l)
                pelzerRel_l = parse(Float64, pelzerRel_l)
                preAdjCor_l = parse(Float64, preAdjCor_l)

                # Add to vector
                push!(station3, station3_l)
                push!(msr, msr_l)
                push!(adj, adj_l)
                push!(cor, cor_l)
                push!(sigmaMSR, sigmaMSR_l)
                push!(sigmaADJ, sigmaADJ_l)
                push!(sigmaCOR, sigmaCOR_l)
                push!(nstat, nstat_l)
                push!(pelzerRel, pelzerRel_l) 
                push!(preAdjCor, preAdjCor_l)

                # Fill in preceding columns
                push!(type, type[end])
                push!(station1, station1[end])
                push!(station2, station2[end])
                push!(par, par[end])

            end

        end
        
        # Setup temporary dataframe (single-line block)
        df2 = DataFrame(
                        type=type,
                        station1=station1,
                        station2=station2,
                        station3=station3,
                        par=par,
                        msr=msr,
                        adj=adj,
                        cor=cor,
                        sigmaMSR=sigmaMSR, 
                        sigmaADJ=sigmaADJ, 
                        sigmaCOR=sigmaCOR, 
                        nstat=nstat, 
                        pelzerRel=pelzerRel, 
                        preAdjCor=preAdjCor, 
        )

        # Combine into one dataframe
        df = vcat(df, df2)
        
    end

    # Return
    return df

end

"""
    dyncor2dataframe(fp)

"""
function dyncor2dataframe(fp)

    lines = readlines(fp)
    block_start = findfirst(contains("Station    "), lines)[1]
    coord_correction_block = lines[block_start+2:end]
    coord_correction_block = split.(coord_correction_block, " ", keepempty=false)
    site = string.(getindex.(coord_correction_block, 1))
    corrE = parse.(Float64, getindex.(coord_correction_block, 10))
    corrN = parse.(Float64, getindex.(coord_correction_block, 11))
    corrU = parse.(Float64,getindex.(coord_correction_block, 12))
    df = DataFrame(
                    site=site,
                    corrE=corrE*1000,
                    corrN=corrN*1000,
                    corrU=corrU*1000,
    )
    df = insertcols(df, :corr3D => sqrt.(df.corrE.^2 + df.corrN.^2 + df.corrU.^2))

    # Return
    return df

end

"""
    dynm2s2dataframe(fp)

"""
function dynm2s2dataframe(fp)

    # Read lines
    lines = readlines(fp)

    # m2s block
    block_start = findfirst(startswith("Station    "), lines)[1]
    block_end = findfirst(startswith("Totals    "), lines)[1]
    m2s_block = lines[block_start+2:block_end-2]
    # Extract columns as substrings
    site = String.(rstrip.(SubString.(m2s_block, 1, 20)))
    A = tryparse.(Int, SubString.(m2s_block, 21, 28))
    B = tryparse.(Int, SubString.(m2s_block, 29, 36))    
    C = tryparse.(Int, SubString.(m2s_block, 37, 44))      
    D = tryparse.(Int, SubString.(m2s_block, 45, 52))      
    E = tryparse.(Int, SubString.(m2s_block, 53, 60))      
    G = tryparse.(Int, SubString.(m2s_block, 61, 68))      
    H = tryparse.(Int, SubString.(m2s_block, 69, 76))       
    I = tryparse.(Int, SubString.(m2s_block, 77, 84))      
    J = tryparse.(Int, SubString.(m2s_block, 85, 92))     
    K = tryparse.(Int, SubString.(m2s_block, 93, 100))      
    L = tryparse.(Int, SubString.(m2s_block, 101, 108))      
    M = tryparse.(Int, SubString.(m2s_block, 109, 116))      
    P = tryparse.(Int, SubString.(m2s_block, 117, 124))      
    Q = tryparse.(Int, SubString.(m2s_block, 125, 132))      
    R = tryparse.(Int, SubString.(m2s_block, 133, 140))       
    S = tryparse.(Int, SubString.(m2s_block, 141, 148))      
    V = tryparse.(Int, SubString.(m2s_block, 149, 156))      
    X = tryparse.(Int, SubString.(m2s_block, 157, 164))      
    Y = tryparse.(Int, SubString.(m2s_block, 165, 172))       
    Z = tryparse.(Int, SubString.(m2s_block, 173, 180))
    total = tryparse.(Int, SubString.(m2s_block, 181, 191))
    # Organise into DataFrame
    df_m2s = DataFrame(
                        site=site,
                        A=A,    
                        B=B,       
                        C=C,       
                        D=D,       
                        E=E,       
                        G=G,       
                        H=H,       
                        I=I,
                        J=J,       
                        K=K,       
                        L=L,       
                        M=M,       
                        P=P,       
                        Q=Q,       
                        R=R,       
                        S=S,       
                        V=V,       
                        X=X,       
                        Y=Y,       
                        Z=Z,
                        total=total,
    )

    # m2s warning block
    block_start = findfirst(startswith("WARNING: "), lines)[1]
    m2s_warnings_block = lines[block_start+4:end]
    site = String.(rstrip.(SubString.(m2s_warnings_block, 1, 20)))
    measurements = String.(rstrip.(SubString.(m2s_warnings_block, 21, 37)))
    count = parse.(Int, SubString.(m2s_warnings_block, 51))
    # Organise into DataFrame
    df_m2s_warnings = DataFrame(
                        site=site,
                        measurements=measurements,
                        count=count,
    )               
    # Return
    return df_m2s, df_m2s_warnings

end

"""
    dyndst2dataframe(fp)

"""
function dyndst2dataframe(fp)

    # Duplicate lines
    lines = readlines(fp)
    near_switch = false
    if findfirst(startswith("Nearby"), lines) != nothing
        
        # Set switch to form near station dataframe too
        near_switch = true

    end

    block_start = findfirst(startswith("Removed "), lines)[1]
    if near_switch == true
        block_end = findfirst(startswith("Nearby"), lines)[1]
        block = lines[block_start+1:block_end-1]
    else
        block = lines[block_start+1:end]
    end
    block = filter(!isempty, block)

    # Extract duplicate site names
    site = []
    for i in eachindex(block)
        temp = SubString(block[i], 5, lastindex(block[i]))
        temp = String(temp) 
        push!(site, temp)
    end

    # Form dataframe
    df_dup = DataFrame(site=site)

    if near_switch == true

        # Nearby lines
        block = lines[block_end+5:end]
        block = filter(!isempty, block)
        lines = nothing

        # Extract nearby site names
        #site = []
        #nearSite = []
        #hrzDelta = []
        #vrtDelta = []
        #for i in eachindex(block)
        #    temp = SubString(block[i], 5, lastindex(block[i]))
        #    temp = String(temp) 
        #    push!(site, temp)
        #end
        site = String.(rstrip.(SubString.(block, 1, 20)))
        nearSite = String.(rstrip.(SubString.(block, 21, 40)))
        hrzDelta = tryparse.(Float64, SubString.(block, 41, 60))
        vrtDelta = tryparse.(Float64, SubString.(block, 61, 80))

        # Form dataframe
        df_near = DataFrame(
                            site=site,
                            nearSite=nearSite,
                            hrzDelta=hrzDelta,
                            vrtDelta=vrtDelta,
        )

    end


    if near_switch == true
        # Return
        return df_dup, df_near
    else
        # Return
        return df_dup 
    end 

end

"""
    dynrename2dataframe(fp)

"""
function dynrename2dataframe(fp)

    # Read lines
    lines = readlines(fp)
    if startswith(lines[1], "!#=")
        block = lines[2:end]
    else
        block = lines[1:end]
    end
    lines = nothing

    # Extract site names
    rename = []
    original = []
    for i in eachindex(block)
        temp_rename = SubString(block[i], 1, 20)
        temp_rename = String(rstrip(temp_rename)) 
        push!(rename, temp_rename)
        temp_original = SubString(block[i], 21, lastindex(block[i]))
        temp_original = String(rstrip(temp_original))
        push!(original, temp_original)
    end

    # Form dataframe
    df = DataFrame(rename=rename, original=original)

    # Return
    return df

end

"""
    dynnear2dataframe(fp)

"""
function dynnear2dataframe(fp)

    # Read lines
    lines = readlines(fp)
    lines = filter(!isempty, lines)

    # Extract site names
    # - site (column 1)
    # - nearSite (column 2)

    site = []
    nearSite = []
    for i in eachindex(lines)
        temp_site = SubString(lines[i], 1, 20)
        temp_site = String(rstrip(temp_site)) 
        push!(site, temp_site)
        temp_nearSite = SubString(lines[i], 21, lastindex(lines[i]))
        temp_nearSite = String(rstrip(temp_nearSite))
        push!(nearSite, temp_nearSite)
    end

    # Form dataframe
    df = DataFrame(site=site, nearSite=nearSite)

    # Return
    return df

end

"""
    dynstnxml2dataframe(fp)

# STATUS: IN DEVELOPMENT

"""
function dynstnxml2dataframe(fp)

    # Read xml document
    doc = readxml(fp)

    # Assign root of xml document
    stn = doc.root

    # Assign vectors
    sites = nodecontent.(findall("//StationCoord/Name", stn))
    constraints = nodecontent.(findall("//Constraints", stn))
    x = nodecontent.(findall("//StationCoord/XAxis", stn))
    y = nodecontent.(findall("//StationCoord/YAxis", stn))
    z = nodecontent.(findall("//StationCoord/Height", stn))
    description = nodecontent.(findall("//Description", stn))

    # Form dataframe
    df = DataFrame(
                    site=sites,
                    constraint=constraints,
                    x=x,
                    y=y,
                    z=z,
                    description=description,
    )

    # Return
    return df

end
