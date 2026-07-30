## ===
# DESCRIPTION
# This file contains multiple functions that relate to common SINEX operations. 
#
# Functions:    - readSINEX_HEADER()
#               - readSINEX_CUSTOM()
#               - readSINEX_FILE_COMMENT()
#               - readSINEX_SITE_ID()
#               - readSINEX_SOLUTION_EPOCHS()
#               - readSINEX_SOLUTION_ESTIMATE()
#               - readSINEX_SOLUTION_MATRIX_ESTIMATE()
#               - sinex2dataframe_SiteIDBlock()
#               - sinex2dataframe_SolutionEpochsBlock()
#               - sinex2dataframe_SolutionEstimateBlock()
#               - sinex2dataframe_SolutionMatrixEstimateBlock()
#               - sinex2dataframe_SolutionDiscontinuityBlock()
#               - dataframe2sinex_SolutionEstimateBlock()
#               - dataframe2sinex_SolutionMatrixEstimateBlock()
#               - df2df_solMatEst_to_vcv()
#               - writeSINEX()
#               - sinex2dataframe()
#               - subset_sites_POLYGON()



"""
    readSINEX_HEADER(fp)

Read the header line of SINEX to a string variable. 

# Arguments
`fp::String`:   Full or relative path to SINEX file. 

"""
function readSINEX_HEADER(fp)
    
    ## ===
    # SINEX HEADER
    snx_header_line = readlines(fp)[1]

    # Return
    return snx_header_line

end

"""
    readSINEX_CUSTOM(fp, start, end_line)

Read custome line range from SINEX.

# Arguments
`fp::String`:   Full or relative path to SINEX file.
`start_line::Int`:   Beginning line number. 
`end_line::Int`:     Finishing line number.

"""
function readSINEX_CUSTOM(fp, start_line, end_line)

    # Range
    i = start_line
    j = end_line

    # Custom range
    custom = readlines(fp)[i:j]

    # Return
    return custom

end

"""
    readSINEX_FILE_COMMENT(fp)

Read the first FILE/COMMENT block of SINEX to a vector of strings. 

# Arguments
`fp::String`:   Full or relative path to SINEX file. 

"""
function readSINEX_FILE_COMMENT(fp)
    
    ## ===
    # FILE/COMMENT BLOCK

    # Read file
    lines_SNX = readlines(fp)

    # Extract lines
    block_start = findfirst(contains("+FILE/COMMENT"), lines_SNX)[1]
    block_end  = findfirst(contains("-FILE/COMMENT"), lines_SNX)[1]
    snx_block_FILE_COMMENT = lines_SNX[block_start:block_end]
    lines_SNX = nothing

    # Return
    return snx_block_FILE_COMMENT

end


"""
    readSINEX_SITE_ID(fp)

Read the SITE/ID block of SINEX to a vector of strings. 

# Arguments
`fp::String`:   Full or relative path to SINEX file. 

"""
function readSINEX_SITE_ID(fp)

    ## ===
    # SITE/ID BLOCK

    # Read file
    lines = readlines(fp)

    # Extract lines
    block_start = findall(contains("+SITE/ID"), lines)[1]
    block_end  = findall(contains("-SITE/ID"), lines)[1]
    snx_block_SITE_ID = lines[block_start:block_end]
    lines = nothing

    # Return
    return snx_block_SITE_ID

end

"""
    readSINEX_SOLUTION_EPOCHS(fp)

Read the SOLUTION/EPOCHS block of SINEX to a vector of strings. 

# Arguments
`fp::String`:   Full or relative path to SINEX file. 

"""
function readSINEX_SOLUTION_EPOCHS(fp)

    ## ===
    # SOLUTION/EPOCHS BLOCK

    # Read file
    lines_SNX = readlines(fp)

    # Extract lines
    block_start = findall(contains("+SOLUTION/EPOCHS"), lines_SNX)[1]
    block_end  = findall(contains("-SOLUTION/EPOCHS"), lines_SNX)[1]
    snx_block_SOLUTION_EPOCHS = lines_SNX[block_start:block_end]
    lines_SNX = nothing

    # Return
    return snx_block_SOLUTION_EPOCHS

end


"""
    readSINEX_SOLUTION_ESTIMATE(fp)

Read the SOLUTION/ESTIMATE block of SINEX to a vector of strings. 

# Arguments
`fp::String`:   Full or relative path to SINEX file. 

"""
function readSINEX_SOLUTION_ESTIMATE(fp)
    
    ## ===
    # SOLUTION/ESTIMATE BLOCK

    # Read file
    lines_SNX = readlines(fp)

    # Extract lines
    block_start = findall(contains("+SOLUTION/ESTIMATE"), lines_SNX)[1]
    block_end  = findall(contains("-SOLUTION/ESTIMATE"), lines_SNX)[1]
    snx_block_SOLUTION_ESTIMATE = lines_SNX[block_start:block_end]
    lines_SNX = nothing

    # Return
    return snx_block_SOLUTION_ESTIMATE

end

"""
    readSINEX_SOLUTION_MATRIX_ESTIMATE(fp)

Read the SOLUTION/MATRIX_ESTIMATE block of SINEX to a vector of strings. 

# Arguments
`fp::String`:   Full or relative path to SINEX file. 

"""
function readSINEX_SOLUTION_MATRIX_ESTIMATE(fp)
    
    ## ===
    # SOLUTION/MATRIX_ESTIMATE BLOCK

    # Read file
    lines_SNX = readlines(fp)

    # Extract lines
    block_start = findall(contains("+SOLUTION/MATRIX_ESTIMATE"), lines_SNX)[1]
    block_end  = findall(contains("-SOLUTION/MATRIX_ESTIMATE"), lines_SNX)[1]
    snx_block_SOLUTION_MATRIX_ESTIMATE = lines_SNX[block_start:block_end]
    lines_SNX = nothing

    # Return
    return snx_block_SOLUTION_MATRIX_ESTIMATE

end

"""
sinex2dataframe_SiteIDBlock(fp)

Convert the SITE/ID block of SINEX file to dataframe. 

# Arguments
`fp::String`:   Full or relative path to SINEX file.

"""
function sinex2dataframe_SiteIDBlock(fp)

    # Read block
    lines = readSINEX_SITE_ID(fp)

    # Handling sinex format variations
    if startswith(lines[2], "*")
        lines = lines[3:end-1]
    else
        lines = lines[2:end-1]
    end

    # Extract columns as substrings
    site = SubString.(lines, 2, 5)
    pt = SubString.(lines, 8, 8)
    domes = SubString.(lines, 10, 18)
    T = SubString.(lines, 20, 20)
    lonDeg = SubString.(lines, 45, 47)
    lonMin = SubString.(lines, 49, 50)
    lonSec = SubString.(lines, 52, 55)
    latDeg = SubString.(lines, 58, 59)
    latMin = SubString.(lines, 61, 62)
    latSec = SubString.(lines, 64, 67)
    ellHgt = SubString.(lines, 69, 75)

    # Correct format
    site = String.(site)
    pt = String.(pt)
    domes = String.(domes)
    T = String.(T)
    lonDeg = parse.(Int, lonDeg[:])
    lonMin = parse.(Int, lonMin[:])
    lonSec = parse.(Float64, lonSec[:])
    latDeg = parse.(Int, latDeg[:])
    latMin = parse.(Int, latMin[:])
    latSec = parse.(Float64, latSec[:])
    ellHgt = parse.(Float64, ellHgt[:])

    # Organise into DataFrame
    df = DataFrame(
                    site=site,
                    pt=pt,
                    domes=domes,
                    T=T,
                    lonDeg=lonDeg,
                    lonMin=lonMin,
                    lonSec=lonSec,
                    latDeg=latDeg,
                    latMin=latMin,
                    latSec=latSec,
                    ellHgt=ellHgt
    )

    # Return
    return df

end

"""
sinex2dataframe_SolutionEpochsBlock(fp)

Convert the SOLUTION/EPOCHS block of SINEX file to dataframe. 

# Arguments
`fp::String`:   Full or relative path to SINEX file.

"""
function sinex2dataframe_SolutionEpochsBlock(fp)

    # Read block
    lines = readSINEX_SOLUTION_EPOCHS(fp)

    # Handling sinex format variations
    if startswith(lines[2], "*")
        lines = lines[3:end-1]
    else
        lines = lines[2:end-1]
    end
    lines = split.(lines, " ", keepempty=false)

    # Isolate into vectors
    site = getindex.(lines,1)
    solnNum = getindex.(lines, 3)
    startEpoch = getindex.(lines,5)
    endEpoch = getindex.(lines,6)
    meanEpoch = getindex.(lines,7)

    # Correct format
    site = String.(site)
    solnNum = String.(solnNum)

    # Convert to DateTime
    new_startEpoch = []
    for t in startEpoch
        strings = split(t, ":")
        i = parse.(Int64, strings)
        epoch = doy2date(i[1], i[2])
        push!(new_startEpoch, epoch)
    end
    new_meanEpoch = []
    for t in meanEpoch
        strings = split(t, ":")
        i = parse.(Int64, strings)
        epoch = doy2date(i[1], i[2])
        push!(new_meanEpoch, epoch)
    end
    new_endEpoch = []
    for t in endEpoch
        strings = split(t, ":")
        i = parse.(Int64, strings)
        epoch = doy2date(i[1], i[2])
        push!(new_endEpoch, epoch)
    end

    # Create time delta
    duration = []
    for t in eachindex(new_startEpoch)
        timeDelta = round(Dates.value(new_endEpoch[t] - new_startEpoch[t])/365.25, digits=2)
        push!(duration, timeDelta) 
    end

    # Organise into DataFrame
    df = DataFrame(
                    site=site,
                    solnNum=solnNum,
                    startEpoch=new_startEpoch,
                    endEpoch=new_endEpoch,
                    meanEpoch=new_meanEpoch,
                    duration=duration
    )

    # Return
    return df

end

"""
    sinex2dataframe_SolutionEstimateBlock(fp)

Convert the SOLUTION/ESTIMATE block of SINEX file to dataframe. 

# Arguments
`fp::String`:   Full or relative path to SINEX file. 

"""
function sinex2dataframe_SolutionEstimateBlock(fp)

    ## ===
    # SOLUTION/ESTIMATE BLOCK

    # Get lines
    lines = readSINEX_SOLUTION_ESTIMATE(fp)

    # Remove non-data lines
    ind = findall(r -> startswith(r, "*") | startswith(r, "+") | startswith(r, "-"), lines)
    lines = deleteat!(lines, ind)

    # Split by column
    lines = split.(lines, " ", keepempty=false)

    # Isolate into vectors
    rowNum = getindex.(lines,1)
    type = getindex.(lines,2)
    site = getindex.(lines,3)
    pt = getindex.(lines,4)
    solnNum = getindex.(lines,5)
    refEpoch = string.(getindex.(lines,6))
    unit = getindex.(lines,7)
    s = getindex.(lines,8)
    est = getindex.(lines,9)
    std = getindex.(lines,10)

    # Correct format
    rowNum = parse.(Int, rowNum[:])
    type = string.(type)
    site = string.(site)
    pt = string.(pt)
    solnNum = string.(solnNum)
    unit = string.(unit)
    s = string.(s)
    est = parse.(Float64, est[:])
    std = parse.(Float64, std[:])

    # Convert reference epoch to DateTime
    new_refEpoch = []
    for t in refEpoch
        strings = split(t, ":")
        i = parse.(Int64, strings)
        epoch = doy2date(i[1], i[2])
        push!(new_refEpoch, epoch)
    end
    refEpoch = new_refEpoch

    # Create site/solution column
    siteSoln = string.(site[:], "_", solnNum[:])

    # Create column for number of solutions (numSolns)
    # - Identify number of parameters
    # - Count the occurences of the site code

    # Number of paramters
    if type[4]=="STAX"
        numPar = 3
    elseif  type[4]=="VELX"
        numPar = 6
    end

    # Number of solutions
    numSolns = []
    for s in site
        number_of_solutions = round(Int, count(==("$s"), site)/numPar)
        push!(numSolns, number_of_solutions)
    end

    # Organise into DataFrame
    df = DataFrame(
                                siteSoln=siteSoln,
                                rowNum=rowNum,
                                par=type,
                                site=site,
                                pt=pt,
                                solnNum=solnNum,
                                numSolns=numSolns,
                                refEpoch=refEpoch,
                                unit=unit,
                                s=s,
                                est=est,
                                sigma=std,
    )

    # Return
    return df

end

"""
    sinex2dataframe_SolutionMatrixEstimateBlock(fp)

Convert the SOLUTION/MATRIX_ESTIMATE block of SINEX file to dataframe. 

# Arguments
`fp::String`:   Full or relative path to SINEX file. 

"""
function sinex2dataframe_SolutionMatrixEstimateBlock(fp)

    ## ===
    # SOLUTION/MATRIX_ESTIMATE BLOCK

    # Read lines
    lines = readSINEX_SOLUTION_MATRIX_ESTIMATE(fp)

    # Identify if upper or lower triangular matrix
    if occursin("U COVA", lines[1])
        triangular = "upper"
    elseif occursin("L COVA", lines[1])
        triangular = "lower"
    end

    @info("$(now()) - SINEX VCV MATRIX IS:                        $triangular triangular")

    # Handling SIENX file variations
    if startswith(lines[2], "*")
        lines = lines[3:end-1]
    else
        lines = lines[2:end-1]
    end
    lines = split.(lines, " ", keepempty=false)

    # Pad out lines with NaN
    for i in eachindex(lines)
        if length(lines[i])==5
            continue
        end
        if length(lines[i]) == 4
            push!(lines[i], "NaN")
        end
        if length(lines[i]) == 3
            push!(lines[i], "NaN")
            push!(lines[i], "NaN")
        end
    end

    # Isolate into vectors
    row = getindex.(lines,1)
    col = getindex.(lines,2)
    q1 = getindex.(lines,3)
    q2 = getindex.(lines,4)
    q3 = getindex.(lines,5)

    # Correct format
    row = parse.(Int, row[:])
    col = parse.(Int, col[:])
    q1 = parse.(Float64, q1[:])
    q2 = parse.(Float64, q2[:])
    q3 = parse.(Float64, q3[:])

    # Create triangularity column
    triangle = fill(triangular, length(q1))

    # Organise into DataFrame
    df = DataFrame(
                        row=row,
                        col=col,
                        q1=q1,
                        q2=q2,
                        q3=q3,
                        triangle=triangle,
    )

    # Return
    return df

end

"""
    sinex2dataframe_SolutionDiscontinuityBlock(fp)

Convert the SOLUTION/DISCONTINUITY block of SINEX file to dataframe. 

# Arguments
`fp::String`:   Full or relative path to SINEX discontinuity file. 

"""
function sinex2dataframe_SolutionDiscontinuityBlock(fp)
    
    # Disconts dataframe
    l = readlines(fp)
    block_start = findfirst(contains("+SOLUTION/DISCONTINUITY"), l)[1]
    block_end = findfirst(contains("-SOLUTION/DISCONTINUITY"), l)[1]
    l = l[block_start+1:block_end-1]
    l = deleteat!(l, findall(contains("*"), l))
    l = split.(l, " ", keepempty=false)    
    site = String.(getindex.(l,1))
    A = String.(getindex.(l,2))
    soln = String.(getindex.(l,3))
    P = String.(getindex.(l,4))
    start = String.(getindex.(l,5))
    finish = string.(getindex.(l,6))
    par = string.(getindex.(l,7))
    df = DataFrame(
                    site=site,
                    A=A,
                    solnNum=soln,
                    start=start,
                    finish=finish,
                    par=par,
    )

    # Return
    return df

end

"""
    dataframe2sinex_SolutionEstimateBlock(df)

Convert dataframe of SINEX SOLUTION/ESTIMATE block to vector of strings ready for writing to SINEX. 

# Arguments
`df::DataFRame`:   Dataframe of SINEX SOLUTION/ESTIMATE block. 

"""
function dataframe2sinex_SolutionEstimateBlock(df)
    
    lines = ["+SOLUTION/ESTIMATE"]
    push!(lines, "*INDEX TYPE__ CODE PT SOLN _REF_EPOCH__ UNIT S __ESTIMATED VALUE____ _STD_DEV___")

    for i in 1:length(df.index)

        l = @sprintf(" %5i %s   %s  %s    %s %s %-3s  %s %.14e %.5e", df.index[i], df.type[i], df.code[i], df.pt[i], df.soln[i], df.refEpoch[i], df.unit[i], df.s[i], df.est[i], df.std[i])
        push!(lines, l)
        
    end

    push!(lines, "-SOLUTION/ESTIMATE")
    
    # Return
    return lines

end

"""
    dataframe2sinex_SolutionMatrixEstimateBlock(df)

Convert dataframe of SINEX SOLUTION/MATRIX_ESTIMATE block to vector of strings ready for writing to SINEX. 

# Arguments
`df::DataFRame`:   Dataframe of SINEX SOLUTION/ESTIMATE block. 

"""
function dataframe2sinex_SolutionMatrixEstimateBlock(df)

    lines = ["+SOLUTION/MATRIX_ESTIMATE L COVA"]
    push!(lines, "*PARA1 PARA2 ____PARA2+0__________ ____PARA2+1__________ ____PARA2+2__________")

    for i in eachindex(df.row)

        l = @sprintf(" %5i %5s %21.14e %21.14e %21.14e", 
                     df.row[i], 
                     df.col[i],
                     df.q1[i],
                     df.q2[i],
                     df.q3[i],
        )
        push!(lines, l)
        
    end
    push!(lines, "-SOLUTION/MATRIX_ESTIMATE L COVA")

    # Remove NaNs
    for i in eachindex(lines)

        lines[i] = replace(lines[i], " NaN" => "")

    end

    # Return
    return lines

end 


"""
    df2df_solMatEst_to_vcv(df_solMatEst, df_solEst)

Convert dataframe of SINEX SOLUTION/MATRIX_ESTIMATE dataframe of site based vcv. 

# Arguments
`df_solMatEst::DataFrame`:   Dataframe of SINEX SOLUTION/MATRIX_ESTIMATE block (sinex2dataframe_SolutionMatrixEstimate()).
`df_solEst::DataFrame`:      Dataframe of SINEX SOLUTION/ESTIMATE block (sinex2dataframe_SolutionEstimate()). 

"""
function df2df_solMatEst_to_vcv(df_solMatEst, df_solEst)

    # Identify triangularity
    triangular = df_solMatEst.triangle[1]

    # Remove velocity parameters, and create new DataFrame
    # - Shouldn't need this anymore.
    #df_solEstXYZ = filter(row -> row.par == "STAX" || row.par == "STAY" || row.par == "STAZ", df_solEstVel)
    #df_solEstVel = nothing

    # Declare empty dataframe
    df = DataFrame(
                    xx = Float64[], 
                    yy = Float64[],
                    zz = Float64[],
                    xy = Float64[],
                    xz = Float64[],
                    yz = Float64[],
    )

    # Convert dataframe to table for speed (due to non-typed columns in dataframes)
    tb_vcv = columntable(df_solMatEst)
    df_solMatEst = nothing

    # Extract vcv for each site solution
    for s in unique(df_solEst.siteSoln)

        # Parameter indices
        x_ind, y_ind, z_ind = filter(name -> name.siteSoln == s, df_solEst)[:,"rowNum"]

        # VCV matrix triangularity
        if triangular == "upper"

            # VCV matrix parameters (variances)
            xx = tb_vcv.q1[findfirst((tb_vcv.row.==x_ind) .& (tb_vcv.col.==x_ind))]
            yy = tb_vcv.q1[findfirst((tb_vcv.row.==y_ind) .& (tb_vcv.col.==y_ind))]
            zz = tb_vcv.q1[findfirst((tb_vcv.row.==z_ind) .& (tb_vcv.col.==z_ind))]
            # VCV matrix parameters (covariances)
            xy = tb_vcv.q2[findfirst((tb_vcv.row.==x_ind) .& (tb_vcv.col.==x_ind))]
            xz = tb_vcv.q3[findfirst((tb_vcv.row.==x_ind) .& (tb_vcv.col.==x_ind))]
            yz = tb_vcv.q2[findfirst((tb_vcv.row.==y_ind) .& (tb_vcv.col.==y_ind))]
        
        elseif triangular == "lower"

            # VCV matrix parameters (variances)
            xx = tb_vcv.q1[findfirst((tb_vcv.row.==x_ind) .& (tb_vcv.col.==x_ind))]
            yy = tb_vcv.q2[findfirst((tb_vcv.row.==y_ind) .& (tb_vcv.col.==x_ind))]
            zz = tb_vcv.q3[findfirst((tb_vcv.row.==z_ind) .& (tb_vcv.col.==x_ind))]
            # VCV matrix parameters (covariances)
            xy = tb_vcv.q1[findfirst((tb_vcv.row.==y_ind) .& (tb_vcv.col.==x_ind))]
            xz = tb_vcv.q1[findfirst((tb_vcv.row.==z_ind) .& (tb_vcv.col.==x_ind))]
            yz = tb_vcv.q2[findfirst((tb_vcv.row.==z_ind) .& (tb_vcv.col.==x_ind))]
        
        end

        # Output
        new_row = [xx yy zz xy xz yz]
        push!(df, new_row)

    end

    # Insert site column
    df = insertcols(df, 1, :site=>unique(df_solEst.siteSoln))
    df = insertcols(df, 2, :siteSoln=>unique(df_solEst.siteSoln))

    # Return
    return df

end  

"""
    writeSINEX(fp, header, comment, SiteID, SolutionEpochs, SolutionEstimate, SolutionMatrixEstimate)

Documentation... 

# Arguments
- 

"""
function writeSINEX(fp, header, comment, SiteID, SolutionEpochs, SolutionEstimate, SolutionMatrixEstimate)

    # Open File
    f = open(fp, "w")

    # Header
    @printf(f, "%s\n", header)

    print(f, "*-------------------------------------------------------------------------------\n")

    # Comment
    @printf(f, "%s\n", comment[1])

    print(f, "*-------------------------------------------------------------------------------\n")

    #SITE/ID
    for i in eachindex(SiteID)

        @printf(f, "%s\n", SiteID[i])

    end

    print(f, "*-------------------------------------------------------------------------------\n")

    # SOLUTION/EPOCHS
    for i in eachindex(SolutionEpochs)

        @printf(f, "%s\n", SolutionEpochs[i])

    end

    print(f, "*-------------------------------------------------------------------------------\n")

    # SOLUTION/ESTIMATE
    for i in eachindex(SolutionEstimate)

        @printf(f, "%s\n", SolutionEstimate[i])

    end

    print(f, "*-------------------------------------------------------------------------------\n")

    # SOLUTION/MATRIX_ESTIMATE
    for i in eachindex(SolutionMatrixEstimate)

        @printf(f, "%s\n", SolutionMatrixEstimate[i])

    end

    # End Line
    print(f, "%ENDSNX")

    close(f)

end

"""
    sinex2dataframe(fp, out_dp=false, numPar=false, dsc_fp=false)

Read in SINEX file and return a dataframe.

# Arguments
`fp::String`:      Full or relative path to SINEX file. 
`out_do::String`:      Output directory path for saving CSV file. Default is to not save CSV.
`p::Int64`:        Custom input for number of SINEX parameters. Default is to identify it from header automatically.  
`dsc_fp::String`:    Full or relative path to discontinuity file (optional). Default is not to add discontinuity column (SITE_YYYYDOY ). 

"""
function sinex2dataframe(fp, out_dp=false, numPar=false, dsc_fp=false)

    @info("$(now()) - BEGIN SINEX CONVERSION:                     sinex2dataframe()")

    ## ===
    # SETUP

    # Input path
    fp = fp

    # Sinex file
    snx_file = split(fp, "/")[end]

    # Disconts file
    if dsc_fp != false
        dsc_file = split(dsc_fp, "/")[end]
    else
        dsc_file = "Not provided"
    end

    # Output path
    if out_dp != false
        csv_dp = out_dp
    else
        csv_dp = "Not provided"
    end

    # SINEX types
    # - SOL (default)
    # - CRD
    CRD = false
    type = split(snx_file, ".")[1]
    type = split(type, "_")[end]
    if type == "CRD"
        CRD = true
    end

    @info("$(now()) - INPUT SINEX FILE:                           $snx_file")   
    @info("$(now()) - INPUT DISCONTS FILE:                        $dsc_file")
    @info("$(now()) - OUTPUT DIRECTORY:                           $csv_dp")

    # Read file header
    snx_header = readSINEX_HEADER(fp)

    # Number of parameters
    if numPar!==false
        numPar = numPar
        @info("$(now()) - CUSTOM SET NUMBER OF PARAMETERS:            $numPar")
    else
        header = snx_header
        header = split(header, " ", keepempty=false)
        if header[end]=="V"
            numPar = 6
            @info("$(now()) - NUMBER OF IDENTIFIED PARAMETERS:            $numPar")
        elseif header[end]=="X"
            numPar = 3
            @info("$(now()) - NUMBER OF IDENTIFIED PARAMETERS:            $numPar")
        elseif header[end]=="S"
            numPar = 6
            @info("$(now()) - WARNING:                                    From SINEX header, solution type = 'S', will assume 6 parameters.")
        end
    end
    
    ## ===
    # SITE/ID BLOCK
    
    @info("$(now()) - PROCESS SITE/ID BLOCK.")
    #df_siteID = sinex2dataframe_SiteIDBlock(fp)

    ## ===
    # SOLUTION/EPOCHS BLOCK

    @info("$(now()) - PROCESS SOLUTION/EPOCHS BLOCK.")

    df_epochs = sinex2dataframe_SolutionEpochsBlock(fp)

    ## ===
    # SOLUTION/ESTIMATE BLOCK

    @info("$(now()) - PROCESS SOLUTION/ESTIMATE BLOCK.")

    df_solEst = sinex2dataframe_SolutionEstimateBlock(fp)

    ## ===
    # SOLUTION/MATRIX_ESTIMATE BLOCK

    if CRD == false

        @info("$(now()) - PROCESS SOLUTION/MATRIX_ESTIMATE BLOCK.")

        df_solMatEst = sinex2dataframe_SolutionMatrixEstimateBlock(fp)

    end

    ## ===
    # PRODUCE SITE BASED VCV DATA FRAME

    if CRD == false

        @info("$(now()) - PRODUCE SITE BASED VCV DATAFRAME.")

        df_vcv = df2df_solMatEst_to_vcv(df_solMatEst, df_solEst)

    end

    ## ===
    # PRODUCE FINAL DATAFRAME

    @info("$(now()) - PRODUCE FINAL DATAFRAME.")

    # Dataframe columns:
    # - site code (site), 
    # - solution (soln),
    # - Combined site solution name (siteSoln), 
    # - reference epoch (refEpoch),
    # - start epoch (startEpoch),
    # - end epoch (endEpoch),
    # - mean epoch (meanEpoch), 
    # - duration in years (years),
    # - solution number (solnNum),
    # - number of solutions (numSolns)
    # - ToDo: document the rest of columns here.
    df = DataFrame(site = filter(row -> row.par == "STAX", df_solEst)[:,"site"])
    df = insertcols(df, 2, :soln => unique(df_solEst.siteSoln))
    df = insertcols(df, 3, :refEpoch => filter(row -> row.par == "STAX", df_solEst)[:,"refEpoch"])
    df = insertcols(df, 4, :startEpoch => df_epochs.startEpoch)
    df = insertcols(df, 5, :endEpoch => df_epochs.endEpoch)
    df = insertcols(df, 6, :meanEpoch => df_epochs.meanEpoch)
    df = insertcols(df, 7, :years => df_epochs.duration)
    df = insertcols(df, 8, :solnNum => filter(row -> row.par == "STAX", df_solEst)[:,"solnNum"])
    df = insertcols(df, 9, :numSolns => filter(row -> row.par == "STAX", df_solEst)[:,"numSolns"])
    
    # Cartesian coordinates and sigmas (XYZ)
    xSol = filter(row -> row.par == "STAX", df_solEst)
    ySol = filter(row -> row.par == "STAY", df_solEst)
    zSol = filter(row -> row.par == "STAZ", df_solEst)
    df = insertcols(df, 10, :x => xSol.est)
    df = insertcols(df, 11, :y => ySol.est)
    df = insertcols(df, 12, :z => zSol.est)
    df = insertcols(df, 13, :sigmaX => xSol.sigma)
    df = insertcols(df, 14, :sigmaY => ySol.sigma)
    df = insertcols(df, 15, :sigmaZ => zSol.sigma)

    # Geographic coordinates (longitude, latitude, and ellipsoid height)
    tb_llh = columntable(xyz2llh.(df.x, df.y, df.z))
    lonDD = tb_llh[1]
    latDD = tb_llh[2]
    ellHgt = tb_llh[3]
    df = insertcols(df, 16, :lon => lonDD)
    df = insertcols(df, 17, :lat => latDD)
    df = insertcols(df, 18, :ellHgt => ellHgt)

    if CRD == false

        ## ===
        # ADD IN VCV and SIGMAS FOR ENU 

        # VCV columns
        df = insertcols(df, 19, :xx => df_vcv.xx)
        df = insertcols(df, 20, :yy => df_vcv.yy)
        df = insertcols(df, 21, :zz => df_vcv.zz)
        df = insertcols(df, 22, :xy => df_vcv.xy)
        df = insertcols(df, 23, :xz => df_vcv.xz)
        df = insertcols(df, 24, :yz => df_vcv.yz)

        # Empty vectors
        ee = Vector{Float64}()
        nn = Vector{Float64}()
        uu = Vector{Float64}()
        en = Vector{Float64}()
        eu = Vector{Float64}()
        nu = Vector{Float64}()

        for i in 1:length(df.xx)

            varEE, varNN, varUU, covarEN, covarEU, covarNU = xyz2enu_vcv(df.xx[i], df.yy[i], df.zz[i], df.xy[i], df.xz[i], df.yz[i], df.lon[i], df.lat[i])

            # Add to vector
            push!(ee, varEE)
            push!(nn, varNN)
            push!(uu, varUU)
            push!(en, covarEN)
            push!(eu, covarEU)
            push!(nu, covarNU)

        end

        # Add in ENU VCV coloumns
        df = insertcols(df, 25, :ee => ee)
        df = insertcols(df, 26, :nn => nn)
        df = insertcols(df, 27, :uu => uu)
        df = insertcols(df, 28, :en => en)
        df = insertcols(df, 29, :eu => eu)
        df = insertcols(df, 30, :nu => nu)

        # Add in ENU standard deviation columns
        # ToDo: investigate if this should be done from df_solEst 
        df = insertcols(df, 31, :sigmaE => sqrt.(ee))
        df = insertcols(df, 32, :sigmaN => sqrt.(nn))
        df = insertcols(df, 33, :sigmaU => sqrt.(uu))
    
    end

    ## ===
    # ORDER
    df = sort(df, order(:soln))

    ## ===
    # ADD DISCONTINUITY COLUMN (optional)
    # - i.e. CODE_YYYDOY

    # Discontinuity dataframe
    if dsc_fp != false

        @info("$(now()) - ADD DISCONTINUITY COLUMN.")

        # Discontinuity dataframe
        df_dsc = sinex2dataframe_SolutionDiscontinuityBlock(dsc_fp)

        # Remove velocity rows
        df_dsc = filter(r -> r.par=="P", df_dsc)

        # Create column for number of solutions (numSolns)
        numSolns = []
        for s in df_dsc.site
            number_of_solutions = count(==("$s"), df_dsc.site)
            push!(numSolns, number_of_solutions)
        end
        df_dsc = insertcols(df_dsc, 4, :numSolns=>string.(numSolns))

        # Add discontinuity column (CODE_YYYYYDOY)
        disconts =  []
        for i in eachindex(df_dsc.site)

            start = df_dsc.start[i]

            # Isolate year, DOY, seconds, and site code
            YY = split(start, ":")[1]
            DOY = split(start, ":")[2]
            CODE = df_dsc.site[i]

            # YY -> YYYY
            if start == "00:000:00000"
                name = string(CODE, "_1900001")
            elseif parse(Int, YY) < 94
                yy = "20"
                YYYY = string(yy, YY)
                name = string(CODE, "_", YYYY, DOY)
            else 
                yy = "19"
                YYYY = string(yy, YY)
                name = string(CODE, "_", YYYY, DOY)
            end
            
            push!(disconts, name)

        end
        
        # Add column
        df_dsc = insertcols(df_dsc, 1, :soln=>disconts)
 
        # Add discontinuity column to main dataframe
        disconts = []
        for s in unique(df.site)

            # Isolate particular station in both sinex and disconts dataframes
            sinex_df = filter(r -> r.site==s, df)
            disconts_df = filter(r -> r.site==s, df_dsc)

            # Apply discont naming to correct solution
            # - not every discont solution will be in sinex solution file
            # - Only apply discont naming to sites that have multiple solutions

            if length(disconts_df.site)==0

                for i in eachindex(sinex_df.soln)
            
                    discont_name = sinex_df.site[i]
                    push!(disconts, discont_name)
                end
            else
                for i in eachindex(sinex_df.site)
                    
                    if disconts_df.numSolns[1]=="1"

                        discont_name = sinex_df.site[i]
                        push!(disconts, discont_name)

                    else

                        sol = sinex_df.solnNum[i]
                        discont_name = filter(r -> r.solnNum==sol, disconts_df)[1, "soln"]
                        push!(disconts, discont_name)

                    end

                end
            end
        end
        
        # Insert solution discontinuity column
        df = rename(df, :soln => :siteSoln)
        df = insertcols(df, 2, :soln => disconts)         

    end

    ## ===
    # CSV FILE
    if out_dp != false
        csv_name = split(snx_file, ".")[1]
        csv_fp = "$out_dp/$(csv_name)_SNX.csv"
        @info("$(now()) - SAVE TO CSV:                              $csv_fp")
        CSV.write("$csv_fp", df)
    end

    @info("$(now()) - FINISH SINEX CONVERSION:                    sinex2dataframe().")

    # Return
    return df

end

"""
    sinex2dataframe_NGCA(dp, JUR, marker_names)

"""
function sinex2dataframe_NGCA(dp, JUR, marker_names)

    # Setup
    dp = dp
    JUR = JUR
    jur = lowercase(JUR)

    # Read each SINEX file to dataframe
    snx_files = filter(endswith(".SNX.$(JUR).NGCA"), readdir("$dp/snx"))
    for i in eachindex(snx_files)

        # Cluster code
        cluster = split(snx_files[i], ".")[1]
        cluster = "$(cluster)_ls"

        snx_name = snx_files[i]
        if i == 1
            df = sinex2dataframe("$dp/snx/$snx_name", false, 3, false)
            df.cluster = fill(cluster, length(df.soln))
        else
            df_temp = sinex2dataframe("$dp/snx/$snx_name", false, 3, false)
            df_temp.cluster = fill(cluster, length(df_temp.soln))
            df = vcat(df, df_temp)
        end
    end

    # Remove non-jurisdiction sites
    ind = []
    for i in eachindex(marker_names)
        marker_name = marker_names[i]
        j = findall(row -> startswith(row, "$(marker_name)"), df.site)
        push!(ind, j)
    end
    ind = reduce(vcat, ind)
    df = df[ind,:]

    # Add: marker number
    tt_file = filter(contains("$(jur)TransTable"), readdir("$dp/other"))[1]
    df_tt = CSV.read("$dp/other/$(tt_file)", DataFrame, delim=',', types=String, header=false, ignorerepeated=false)
    marker_number = String[]
    for i in eachindex(df.site)
        j = findfirst(isequal(df.site[i]), df_tt[:,1])
        if isnothing(j)
            number = df.site[i]
        else
            number = df_tt[j,2]
        end
        push!(marker_number, number)
    end
    df = insertcols(df, 2, :number => marker_number)

    # Refine: soln column
    solnCol = []
    for s in unique(df.soln)

        df_temp = filter(r -> r.soln=="$s", df)
        for i in eachindex(df_temp.soln)
            suffix = df_temp.refEpoch[i]
            solnName = "$(df_temp.site[1])_$suffix"
            push!(solnCol, solnName)
        end
    end
    df.soln = solnCol

    # Refine: soln column again
    # - To account for multiple occupations per day
    df = sort(df, order(:soln))
    solnCol = []
    letters = collect('A':'Z')
    for s in unique(df.soln)

        # Per solution dataframe
        df_temp = filter(r->r.soln=="$s", df)

        if length(df_temp.soln)==1

            push!(solnCol, df_temp.soln[1])

        else

            for i in eachindex(df_temp.soln)

                solution = df_temp.soln[i]
                letter = letters[i]
                new_soln_name = "$(solution)_$(letter)"
                push!(solnCol, new_soln_name)
            end   
        end
    end 
    df.soln = solnCol

    # Return
    return df

end

"""
    subset_sites_POLYGON(df, lonBounds, latBounds)

# STATUS: IN DEVELOPMENT

"""
function subset_sites_POLYGON(df, lonBounds, latBounds)

    # Create polygon nodes
    polygon = [lonBounds latBounds]

    # Create data points
    data = [df.lon df.lat]

    # Check for inbounds condition
    inBounds = inpoly2(data, polygon)

    # Get indices
    indices_inbound = findall((inBounds[:,1] .== true) .|| (inBounds[:,2] .== true))

    # Subset dataframe to inbound rows only
    df = df[indices_inbound, :]

    # Return
    return df

end