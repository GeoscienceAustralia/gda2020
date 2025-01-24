## ===
# DESCRIPTION
# This file contains useful functions to aid QA/QC procedures for GDA2020.
# 
# Functions:    - coordinate_delta_SINEX()
#               - uncertainty_delta_SINEX()
#               - coordinate_delta_DYNADJUST()
#               - coordinate_delta_to_RVS()
#               - plot_delta_PDF()
#               - plot_delta_BAR()
#               - disconts2dataframe()




"""
    coordinate_delta_SINEX(df1, df2)

Compute coordinate differences between two dataframes. 
- Originally designed to be used for dataframes formed from sinex2dataframe(). 
- But it can work with any dataframe that contains the correct fields. 
- Dataframe needs at minimum these fields: `site`, `code`, `lon`, `lat`, `x`, `y`, `z`.

# Arguments
`df1::DataFrame`:    Dataframe of dataset one.
`df2::DataFrame`:    Dataframe of dataset two.

# Returns
- Dataframe of coordinate differences in ENU system (mm).

"""
function coordinate_delta_SINEX(df1, df2)

    @info("$(now()) - BEGIN COORDINATE COMPARISON:                coordinate_delta_SINEX()")     

    # ===
    # LIMIT TO COMMON SITES ONLY
    # - Remove the added site solutions (from dataset 1)
    # - Remove the removed site solutons (from dataset 2)
    additional = setdiff(df1[:, "soln"], df2[:, "soln"])
    removed = setdiff(df2[:, "soln"], df1[:, "soln"])
    for s in additional
        df1 = filter(r -> r.soln != s, df1)
    end
    for s in removed
        df2 = filter(r -> r.soln != s, df2)
    end

    @info("$(now()) - COMMON SOLUTIONS:                           $(length(df1.soln))")
    @info("$(now()) - NON-COMMON SOLUTIONS (EXCLUDED):            $(length(additional))")

    # Order
    sort!(df1, order(:soln))
    sort!(df2, order(:soln))

    ## ===
    # COMPUTE CARTESIAN COORDINATE DIFFERENCES  
    # - account for multiple solutions per site.
    # - also create comparison solution column

    # Coordinate delta
    dX = df1.x .- df2.x
    dY = df1.y .- df2.y
    dZ = df1.z .- df2.z

    # Solution used for comparison (siteDelta)
    solnCompare = df2.soln

    ## ===
    # TRANSFORM CARTESIAN COORDINATE DIFFERENCES TO LOCAL TOPOCENTRIC
    # - XYZ -> ENU.
    tb_enuDelta = columntable(xyz2enu.(dX, dY, dZ, df2.lon, df2.lat))
    dE = tb_enuDelta[1]
    dN = tb_enuDelta[2]
    dU = tb_enuDelta[3]

    ## ===
    # TRANSFORM DATASET(2) CARTESIAN UNCERTAINTIES TO LOCAL TOPOCENTRIC
    # - XYZ -> ENU
    # - Convert to 2-sigma
    covar=zeros(length(df2.sigmaX))
    tb_1sigma_s2_ENU = columntable(xyz2enu_vcv.(df2.sigmaX, df2.sigmaY, df2.sigmaZ, covar, covar, covar, df2.lon, df2.lat))
    d2_std95_E = 1.96*tb_1sigma_s2_ENU[1]
    d2_std95_N = 1.96*tb_1sigma_s2_ENU[2]
    d2_std95_U = 1.96*tb_1sigma_s2_ENU[3]

    ## ===
    # UNCERTAINTY PROPAGATION OF COORDINATE DIFFERENCES

    # Dataset(1) coordinate (1sigma)
    d1_stdX = df1.sigmaX
    d1_stdY = df1.sigmaY
    d1_stdZ = df1.sigmaZ

    # Dataset(2) coordinate (1sigma)
    d2_stdX = df2.sigmaX
    d2_stdY = df2.sigmaY
    d2_stdZ = df2.sigmaZ

    # Propagation of coordinate uncertainties to coordinate delta uncertainties
    # - also converted to 2-Sigma (95% CI)
    delta_uncX = 1.96*sqrt.(d2_stdX.^2 + d1_stdX.^2)
    delta_uncY = 1.96*sqrt.(d2_stdY.^2 + d1_stdY.^2)
    delta_uncZ = 1.96*sqrt.(d2_stdZ.^2 + d1_stdZ.^2)

    # Transform to local topocentric coordinate uncertainties
    # - no covariance here (could add later). 
    covar=zeros(length(delta_uncX))
    tb_delta_std95 = columntable(xyz2enu_vcv.(delta_uncX, delta_uncY, delta_uncZ, covar, covar, covar, df2.lon, df2.lat))
    delta_std95_E = tb_delta_std95[1]
    delta_std95_N = tb_delta_std95[2]
    delta_std95_U = tb_delta_std95[3]

    ## ===
    # CONDITION OF STATISTICAL SIGNIFICANCE
    # - coordinate deltas outside the bounds of 95% CI.
    # - outside bounds = true, inside bounds = false.

    # East
    outlierE = []
    for i in eachindex(dE)
        c = abs(dE[i]) - delta_std95_E[i]
        if c <= d2_std95_E[i]
            push!(outlierE, false)
        else
            push!(outlierE, true)
        end
    end

    # North
    outlierN = []
    for i in eachindex(dN)
        c = abs(dN[i]) - delta_std95_N[i]
        if c <= d2_std95_N[i]
            push!(outlierN, false)
        else
            push!(outlierN, true)
        end
    end

    # Up
    outlierU = []
    for i in eachindex(dU)
        c = abs(dU[i]) - delta_std95_U[i]
        if c <= d2_std95_U[i]
            push!(outlierU, false)
        else
            push!(outlierU, true)
        end
    end

    ## ===
    # FORM DATAFRAME OF COORDINATE DIFFERENCES FROM OTHER SOLUTION
    # - convert m --> mm.
    # - Add in distance column.

    # Compute distance
    distance = sqrt.(dE.^2 + dN.^2 + dU.^2)

    # Dataframe
    df = DataFrame(
                    soln=df1.soln, 
                    site=df1.site, 
                    solnCompare=solnCompare,
                    numSolns=df1.numSolns,
                    dE=dE*1000,
                    dN=dN*1000,
                    dU=dU*1000,
                    distance3D=distance*1000,
                    d2_std95_E=d2_std95_E*1000,
                    d2_std95_N=d2_std95_N*1000,
                    d2_std95_U=d2_std95_U*1000,
                    delta_std95_E=delta_std95_E*1000,
                    delta_std95_N=delta_std95_N*1000,
                    delta_std95_U=delta_std95_U*1000,
                    outlierE=outlierE,
                    outlierN=outlierN,
                    outlierU=outlierU,
    )

    # AUSPOS cluster column
    if "cluster" in names(df1)
        df.cluster = df1.cluster
    end
    
    ## ===
    # RETURN
    @info("$(now()) - FINISH COORDINATE COMPARISON:               coordinate_delta_SINEX()") 

    return df

end

"""
    uncertainty_delta_SINEX(df1, df2)

Compute uncertainty differences between two dataframes. 
- Originally designed to be used for dataframes formed from sinex2dataframe(). 
- But it can work with any dataframe that contains the correct fields. 
- Dataframe needs at minimum these fields: `site`, `code`, `sigmaE`, `sigmaN`, `sigmaU`.

# Arguments
`df1::DataFrame`:    Dataframe of dataset one.
`df2::DataFrame`:    Dataframe of dataset two.

# Returns
- Dataframe of uncertainty differences in ENU system (mm).

"""
function uncertainty_delta_SINEX(df1, df2)

    @info("$(now()) - BEGIN UNCERTAINTY COMPARISON:               uncertainty_delta_SINEX()") 

    # ===
    # LIMIT TO COMMON SITES ONLY
    # - Remove the added sites (from dataset 1)
    # - Remove the removed sites (from dataset 2)
    additional = setdiff(df1[:, "soln"], df2[:, "soln"])
    for s in additional
        df1= filter(r -> r.soln != s, df1)
    end
    removed = setdiff(df2[:, "soln"], df1[:, "soln"])
    for s in removed
        df2= filter(r -> r.soln != s, df2)
    end

    @info("$(now()) - COMMON SOLUTIONS:                           $(length(df1.soln))")
    @info("$(now()) - NON-COMMON SOLUTIONS (EXCLUDED):            $(length(additional))")

    sort!(df1, order(:soln))
    sort!(df2, order(:soln))

    ## ===
    # COMPUTE TOPOCENTRIC UNCERTAINTY DIFFERENCES  
    # - account for multiple solutions per site. 
    dE = df1.sigmaE .- df2.sigmaE
    dN = df1.sigmaN .- df2.sigmaN
    dU = df1.sigmaU .- df2.sigmaU

    ## ===
    # FORM DATAFRAME OF UNCERTAINTY DIFFERENCES FROM OTHER NADJ
    # - convert m --> mm.
    # - Add in distance column.

    # Compute distance
    distance = sqrt.(dE.^2 + dN.^2 + dU.^2)

    # Dataframe
    df = DataFrame(
                    soln=df1.soln,
                    site=df1.site,
                    dE=dE*1000,
                    dN=dN*1000,
                    dU=dU*1000,
                    distance3D=distance*1000,
    )

    ## ===
    # RETURN
    @info("$(now()) - FINISH UNCERTAINTY COMPARISON:              uncertainty_delta_SINEX()")

    return df 

end

"""
    coordinate_delta_DYNADJUST(df1, df2)

    Compute the coordinate differences between two national adjustment (NADJ) dataframes.
    - Created to work on dataframes formed from NADJ solutions.
    
    # Arguments
    `df1::DataFrame`:      DataFrame of solution 1 (NADJ). 
    `df2::DataFrame`:      DataFrame of solution 2 (NADJ). 
    
    # Returns
    - Dataframe of coordinate differences in ENU system (mm).

"""
function coordinate_delta_DYNADJUST(df1, df2)

    @info("$(now()) - BEGIN COORDINATE COMPARISON:                coordinate_delta_DYNADJUST()") 

    # ===
    # LIMIT TO COMMON SITES ONLY
    # - Remove the added site solutions (from dataset 1)
    # - Remove the removed site solutons (from dataset 2)
    additional = setdiff(df1[:, "soln"], df2[:, "soln"])
    for s in additional
        df1= filter(r -> r.soln != s, df1)
    end
    removed = setdiff(df2[:, "soln"], df1[:, "soln"])
    for s in removed
        df2= filter(r -> r.soln != s, df2)
    end

    @info("$(now()) - COMMON SOLUTIONS:                           $(length(df1.soln))")
    @info("$(now()) - NON-COMMON SOLUTIONS (EXCLUDED):            $(length(additional))")

    # Order
    sort!(df1, order(:soln))
    sort!(df2, order(:soln))

    ## ===
    # COMPUTE CARTESIAN COORDINATE DIFFERENCES  
    # - account for multiple solutions per site. 
    dX = df1.x .- df2.x
    dY = df1.y .- df2.y
    dZ = df1.z .- df2.z

    ## ===
    # TRANSFORM CARTESIAN COORDINATE DIFFERENCES TO LOCAL TOPOCENTRIC
    # - XYZ -> ENU.
    tb_enuDelta = columntable(xyz2enu.(dX, dY, dZ, df2.lon, df2.lat))
    dE = tb_enuDelta[1]
    dN = tb_enuDelta[2]
    dU = tb_enuDelta[3]

    ## ===
    # CONVERT NADJ-2 1-SIGMA to 2-SIGMA (95% CI)
    d2_std95_E = 1.96*df2.sigmaE
    d2_std95_N = 1.96*df2.sigmaN
    d2_std95_U = 1.96*df2.sigmaU

    ## ===
    # UNCERTAINTY PROPAGATION OF COORDINATE DIFFERENCES

    # SINEX-1 coordinate (1sigma)
    d1_stdE = df1.sigmaE
    d1_stdN = df1.sigmaN
    d1_stdU = df1.sigmaU

    # SINEX-2 coordinate (1sigma)
    d2_stdE = df2.sigmaE
    d2_stdN = df2.sigmaN
    d2_stdU = df2.sigmaU

    # Propagation of coordinate uncertainties to coordinate delta uncertainties
    # - also convert to 2-Sigma (95% CI)
    delta_std95_E = 1.96*sqrt.(d2_stdE.^2 + d1_stdE.^2)
    delta_std95_N = 1.96*sqrt.(d2_stdN.^2 + d1_stdN.^2)
    delta_std95_U = 1.96*sqrt.(d2_stdU.^2 + d1_stdU.^2)

    ## ===
    # CONDITION OF STATISTICAL SIGNIFICANCE
    # - coordinate deltas outside the bounds of 95% CI.
    # - outside bounds = true, inside bounds = false.

    # East
    outlierE = []
    for i in eachindex(dE)
        c = abs(dE[i]) - delta_std95_E[i]
        if c <= d2_std95_E[i]
            push!(outlierE, false)
        else
            push!(outlierE, true)
        end
    end

    # North
    outlierN = []
    for i in eachindex(dN)
        c = abs(dN[i]) - delta_std95_N[i]
        if c <= d2_std95_N[i]
            push!(outlierN, false)
        else
            push!(outlierN, true)
        end
    end

    # Up
    outlierU = []
    for i in eachindex(dU)
        c = abs(dU[i]) - delta_std95_U[i]
        if c <= d2_std95_U[i]
            push!(outlierU, false)
        else
            push!(outlierU, true)
        end
    end

    ## ===
    # FORM DATAFRAME OF COORDINATE DIFFERENCES FROM OTHER NADJ
    # - convert m --> mm.
    # - Add in distance column.

    # Compute distance
    distance = sqrt.(dE.^2 + dN.^2 + dU.^2)

    # Dataframe
    df = DataFrame(
                    soln=df1.soln,
                    site=df1.site,
                    dE=dE*1000,
                    dN=dN*1000,
                    dU=dU*1000,
                    distance3D=distance*1000,
                    d2_std95_E=d2_std95_E*1000,
                    d2_std95_N=d2_std95_N*1000,
                    d2_std95_U=d2_std95_U*1000,
                    delta_std95_E=delta_std95_E*1000,
                    delta_std95_N=delta_std95_N*1000,
                    delta_std95_U=delta_std95_U*1000,
                    outlierE=outlierE,
                    outlierN=outlierN,
                    outlierU=outlierU,
    )

    ## ===
    # RETURN
    @info("$(now()) - FINISH COORDINATE COMPARISON:               coordinate_delta_DYNADJUST()")

    return df

end

"""
    coordinate_delta_to_RVS(df, df_rvs)

Compute coordinate differences between dataframe and Recognised Value Standard (RVS). 
- Originally designed to be used for dataframes formed from sinex2dataframe(). 
- But it can work with any dataframe that contains the correct fields. 
- Dataframe needs at minimum these fields: `site`, `code`, `lon`, `lat`, `x`, `y`, `z`.

# Arguments
`df::DataFrame`:     Dataframe of dataset one.
`rvs::DataFrame`:    Dataframe of RVS.

# Returns
- Dataframe of coordinate differences in ENU system (mm).

# To Do
- Think about the numSolns requirement in dataframe between APREF/NADJ summary.

"""
function coordinate_delta_to_RVS(df, df_rvs)

    @info("$(now()) - BEGIN COORDINATE COMPARISON:                coordinate_delta_to_RVS()")

    # Assign dataframes
    # - df1 (d1) = dataframe of interest
    # - df2 (d2) = dataframe of comparison (rvs)
    df1 = df
    df2 = df_rvs

    ## ===
    # REMOVE NON-COMMON SITES
    # - Remove sites from NADJ dataframe that are not present in RVS
    ind = []
    for s in df2.site
        i = findall(r -> startswith(r, "$(s)_") || isequal(r, "$(s)"), df1.soln)
        push!(ind, i)
    end
    ind = reduce(vcat, ind)
    df1 = df1[ind,:]

    # Order
    sort!(df1, order(:site))
    sort!(df2, order(:site))

    ## ===
    # COMPUTE COORDINATE DIFFERENCES
    # - Account for multiple solutions per site
    dX = []
    dY = []
    dZ = []
    for d2_site in df2.site
        
        i = findall(r -> r[1:4]==d2_site, df1.soln)

        for index in i

            # X component
            x_delta = df1.x[index] - df2[df2.site.==d2_site, :].X[1]
            push!(dX, x_delta)
            # Y component
            y_delta = df1.y[index] - df2[df2.site.==d2_site, :].Y[1]
            push!(dY, y_delta)
            # Z component
            z_delta = df1.z[index] - df2[df2.site.==d2_site, :].Z[1]
            push!(dZ, z_delta)
        end
    end

    ## ===
    # GATHER RVS COORDINATES AND UNCERTAINTIES FOR OUTPUT DATAFRAME
    # - loops are accounting for dimension mismatch between RVS and NADJ dataframe.

    # Cartesian coordinates
    rvsX = []
    rvsY = []
    rvsZ = []
    for d2_site in df2.site

        i = findall(r -> r[1:4]==d2_site, df1.soln)

        count = 1
        while count <= length(i)

            # X component
            local x = df2[df2.site.==d2_site, :].X[1]
            push!(rvsX, x)
            # Y component
            local y = df2[df2.site.==d2_site, :].Y[1]
            push!(rvsY, y)
            # X component
            local z = df2[df2.site.==d2_site, :].Z[1]
            push!(rvsZ, z)

            count += 1

        end
    end

    # Cartesian coordinate uncertainties
    d2_uncX = []
    d2_uncY = []
    d2_uncZ = []
    for d2_site in df2.site

        i = findall(r -> r[1:4]==d2_site, df1.soln)

        count = 1
        while count <= length(i)
            
            # X component
            xUnc = df2[df2.site.==d2_site, :]."95% CI X"[1]
            push!(d2_uncX, xUnc)

            # Y component
            yUnc = df2[df2.site.==d2_site, :]."95% CI Y"[1]
            push!(d2_uncY, yUnc)

            # X component
            zUnc = df2[df2.site.==d2_site, :]."95% CI Z"[1]
            push!(d2_uncZ, zUnc)

            count += 1

        end
    end

    #  ## ===
    # TRANSFORM CARTESIAN COORDINATES TO LOCAL TOPOCENTRIC
    # - no covariance here (could add later). 
    # - need to calculate lat/lon for RVS.

    # XYZ -> LLH
    tb_d2_LLH = columntable(xyz2llh.(rvsX, rvsY, rvsZ))
    d2_λ = tb_d2_LLH[1]
    d2_ϕ = tb_d2_LLH[2]

    # XYZ -> ENU
    covar=zeros(length(d2_uncX))
    tb_2sigma_d2_ENU = columntable(xyz2enu_vcv.(d2_uncX, d2_uncY, d2_uncZ, covar, covar, covar, d2_λ, d2_ϕ))
    d2_std95_E = tb_2sigma_d2_ENU[1]
    d2_std95_N = tb_2sigma_d2_ENU[2]
    d2_std95_U = tb_2sigma_d2_ENU[3]

    ## ===
    # TRANSFORM CARTESIAN COORDINATE DIFFERENCES TO LOCAL TOPOCENTRIC
    # - XYZ -> ENU.
    tb_enuDelta = columntable(xyz2enu.(dX, dY, dZ, d2_λ, d2_ϕ))
    dE = tb_enuDelta[1]
    dN = tb_enuDelta[2]
    dU = tb_enuDelta[3]

    ## ===
    # UNCERTAINTY PROPAGATION OF COORDINATE DIFFERENCES

    # RVS coordinate (1sigma)
    d2_stdE = d2_std95_E/1.96
    d2_stdN = d2_std95_N/1.96
    d2_stdU = d2_std95_U/1.96

    # NADJ coordinate (1sigma)
    d1_stdE = df1.sigmaE
    d1_stdN = df1.sigmaN
    d1_stdU = df1.sigmaU

    # Propagation of coordinate uncertainties to coordinate delta uncertainties
    # - also converted to 2-Sigma (95% CI)
    delta_std95_E = 1.96*sqrt.(d2_stdE.^2 + d1_stdE.^2)
    delta_std95_N = 1.96*sqrt.(d2_stdN.^2 + d1_stdN.^2)
    delta_std95_U = 1.96*sqrt.(d2_stdU.^2 + d1_stdU.^2)

    ## ===
    # CONDITION OF STATISTICAL SIGNIFICANCE
    # - coordinate deltas outside the bounds of 95% CI.
    # - outside bounds = true, inside bounds = false.

    # East
    outlierE = []
    for i in eachindex(dE)
        c = abs(dE[i]) - delta_std95_E[i]
        if c <= d2_std95_E[i]
            push!(outlierE, false)
        else
            push!(outlierE, true)
        end
    end

    # North
    outlierN = []
    for i in eachindex(dN)
        c = abs(dN[i]) - delta_std95_N[i]
        if c <= d2_std95_N[i]
            push!(outlierN, false)
        else
            push!(outlierN, true)
        end
    end

    # Up
    outlierU = []
    for i in eachindex(dU)
        c = abs(dU[i]) - delta_std95_U[i]
        if c <= d2_std95_U[i]
            push!(outlierU, false)
        else
            push!(outlierU, true)
        end
    end

    # Account for dataframe with no numSolns column
    if !in("numSolns", names(df1))
        numSolns = fill(NaN, length(df1.soln))
        df1 = insertcols(df1, :numSolns => numSolns)
    end

    ## ===
    # FORM DATAFRAME OF COORDINATE DIFFERENCES FROM RVS
    # - Add in distance column.
    # - convert m --> mm.
    distance = sqrt.(dE.^2 + dN.^2 + dU.^2)
    df = DataFrame(
                    soln=df1.soln, 
                    site=df1.site,
                    numSolns=df1.numSolns,
                    dE=dE*1000,
                    dN=dN*1000,
                    dU=dU*1000,
                    distance3D=distance*1000,
                    d2_std95_E=d2_std95_E*1000,
                    d2_std95_N=d2_std95_N*1000,
                    d2_std95_U=d2_std95_U*1000,
                    delta_std95_E=delta_std95_E*1000,
                    delta_std95_N=delta_std95_N*1000,
                    delta_std95_U=delta_std95_U*1000,
                    outlierE=outlierE,
                    outlierN=outlierN,
                    outlierU=outlierU,
    )

    ## ===
    # RETURN
    @info("$(now()) - FINISH COORDINATE COMPARISON:               coordinate_delta_RVS()")

    return df

end

"""
    plot_delta_PDF(df, title)

Plot Probability Density Function (PDF) of coordinate or uncertainty differences. 

# Arguments
`df::DataFrame`:    DataFrame of coordinate or uncertainty differences.
`title::String`:    String of plot title.

# Returns
- Figure of PDF for differences in ENU system (mm).

"""
function plot_delta_PDF(df, title)

    # East componant
    p1 = density(df.dE,
    ylabel="Density",
    xlabel="\\DeltaE (mm)",
    )
    annotate!((0.7,0.95), text("Count:    $(length(df.dE))", 12, :left))
    annotate!((0.7,0.9), text(@sprintf("Min:       % .2f", minimum(df.dE)), 12, :left))
    annotate!((0.7,0.85), text(@sprintf("Max:      % .2f", maximum(df.dE)), 12, :left))
    annotate!((0.7,0.8), text(@sprintf("Median:  % .2f", median(df.dE)), 12, :left))
    annotate!((0.7,0.75), text(@sprintf("Mean:    % .2f", mean(df.dE)), 12, :left))
    annotate!((0.7,0.7), text(@sprintf("Sigma:   % .2f", std(df.dE)), 12, :left))
    # North componant
    p2 = density(df.dN,
            ylabel="Density",
            xlabel="\\DeltaN (mm)",
    )
    annotate!((0.7,0.95), text("Count:    $(length(df.dN))", 12, :left))
    annotate!((0.7,0.9), text(@sprintf("Min:       % .2f", minimum(df.dN)), 12, :left))
    annotate!((0.7,0.85), text(@sprintf("Max:      % .2f", maximum(df.dN)), 12, :left))
    annotate!((0.7,0.8), text(@sprintf("Median:  % .2f", median(df.dN)), 12, :left))
    annotate!((0.7,0.75), text(@sprintf("Mean:    % .2f", mean(df.dN)), 12, :left))
    annotate!((0.7,0.7), text(@sprintf("Sigma:   % .2f", std(df.dN)), 12, :left))
    # Up componant
    p3 = density(df.dU,
            ylabel="Density",
            xlabel="\\DeltaU (mm)",
    )
    annotate!((0.7,0.95), text("Count:    $(length(df.dU))", 12, :left))
    annotate!((0.7,0.9), text(@sprintf("Min:       % .2f", minimum(df.dU)), 12, :left))
    annotate!((0.7,0.85), text(@sprintf("Max:      % .2f", maximum(df.dU)), 12, :left))
    annotate!((0.7,0.8), text(@sprintf("Median:  % .2f", median(df.dU)), 12, :left))
    annotate!((0.7,0.75), text(@sprintf("Mean:    % .2f", mean(df.dU)), 12, :left))
    annotate!((0.7,0.7), text(@sprintf("Sigma:   % .2f", std(df.dU)), 12, :left))
    # Finalise figure
    fig = plot(p1, 
        p2, 
        p3, 
        layout=(1,3), 
        size=(1500, 500), 
        margin=15Plots.mm,
        plot_title=title,
        plot_titlefontsize=16,
        xlim=(-25.0, 25.0),
        #ylim=(-0.01, 0.2),
        xrotation=50,
        minorticks=5,
        tickfontsize=12,
        labelfontsize=16,
        linewidth=3,
        legend=false,
        grid=true,
        gridalpha=0.6
    )

    # Return
    return fig

end


"""
    plot_delta_BAR(df, title)

Plot bar graph of coordinate differences. 

# Arguments
`df::DataFrame`:    DataFrame of coordinate differences.
`title::String`:    String of plot title. 

"""
function plot_delta_BAR(df, title)
    
    # East component
    p1 = plot(df.soln, zeros(length(df.soln)),
                ribbon=df.d2_std95_E,    
                linestrokewidth=0,
                c=:grey,
                fillalpha=0.4,
                label=false,
    ) 
    plot!(df.soln, df.dE, yerr=df.delta_std95_E, 
                seriestype=:bar,
                c=:lightblue,
                markerstrokecolor=:tan,
                ylabel="\\Delta East (mm)",
                xlabel="Station (and solution number)",
                label=false,
    )
    annotate!((1.0,1.05), text(@sprintf("μ:  %.2f mm,  σ:  %.2f mm", mean((df.dE)), std((df.dE))), 9, :right))
    annotate!((0.01,1.05), text("GREY AREA: comparison dataset coordinate uncertainty (95% CI),", 9, :left))
    annotate!((0.4,1.05), text("ERROR BARS: coordinate difference uncertainty (95% CI)", 9, :left)) 
    # North component
    p2 = plot(df.soln, zeros(length(df.soln)),          
                fillrange=(-df.d2_std95_N, df.d2_std95_N),
                linestrokewidth=0,
                c=:grey,
                fillalpha=0.4,
                label=false,
    ) 
    plot!(df.soln, df.dN, yerr=df.delta_std95_N, 
                seriestype=:bar,
                c=:lightblue,
                markerstrokecolor=:tan,
                ylabel="\\Delta North (mm)",
                xlabel="Station (and solution number)",
                label=false,
    )
    annotate!((1.0,1.05), text(@sprintf("μ:  %.2f mm,  σ:  %.2f mm", mean((df.dN)), std((df.dN))), 9, :right))
    annotate!((0.01,1.05), text("GREY AREA: comparison dataset coordinate uncertainty (95% CI),", 9, :left))
    annotate!((0.4,1.05), text("ERROR BARS: coordinate difference uncertainty (95% CI)", 9, :left)) 
    # Up component
    p3 = plot(df.soln, zeros(length(df.soln)),          
                fillrange=(-df.d2_std95_U, df.d2_std95_U),
                linestrokewidth=0,
                c=:grey,
                fillalpha=0.4,
                label=false,
    ) 
    plot!(df.soln, df.dU, yerr=df.delta_std95_U, 
                seriestype=:bar,
                c=:lightblue,
                markerstrokecolor=:tan,
                ylabel="\\Delta Up (mm)",
                xlabel="Station (and solution number)",  
                label=false,
    )
    annotate!((1.0,1.05), text(@sprintf("μ:  %.2f mm,  σ:  %.2f mm", mean((df.dU)), std((df.dU))), 9, :right))
    annotate!((0.01,1.05), text("GREY AREA: comparison dataset coordinate uncertainty (95% CI),", 9, :left))
    annotate!((0.4,1.05), text("ERROR BARS: coordinate difference uncertainty (95% CI)", 9, :left))    
    # Finalise figure
    fig = plot(p1, 
                p2, 
                p3,
                layout=(3,1), 
                size=(1500, 1000), 
                top_margin=2Plots.mm,
                bottom_margin=7Plots.mm,
                left_margin=7Plots.mm,
                plot_title=title,
                plot_titlefontsize=16,
                dpi=300,
                xwiden=0.96,
                xtickfontsize=4,
                ytickfontsize=12,
                labelfontsize=16,
                xrotation=60,
                grid=true,
                gridalpha=0.4,
                xticks=:all,
                ylims=(-60, 60),
    )

    # Return
    return fig

end
