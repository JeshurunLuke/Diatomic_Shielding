using PlotlyJS

using Colors

function set_figure_style!(p; title = "", ylabel = "", xlabel ="")
    relayout!(p,
        title=title, 
        xaxis=attr(
            title=xlabel, 
            showline=true,
            linecolor="black",
            mirror=true,
            ticks="outside",
            gridcolor="lightgrey"
        ),
        yaxis=attr(
            title=ylabel, 
            showline=true,
            linecolor="black",
            mirror=true,
            ticks="outside",
            gridcolor="lightgrey"
        ),
        plot_bgcolor="white",
        width=1000,
        height=400
    )
end
plotScan(mol::moleculeInteraction, sol_vec; stateOI = [nothing, nothing, nothing], units = "K", maxDisplay = 300) = plotScan(PlotlyJS.plot(), mol::moleculeInteraction, sol_vec; stateOI = stateOI, units = units, maxDisplay = maxDisplay)
function plotScan(p::PlotlyJS.SyncPlot, mol::moleculeInteraction, sol_vec; stateOI = [nothing, nothing, nothing], units = "K", maxDisplay = 300)
    yConv = Dict("K"=> 1/1e6, "uK" => 1)

    xlabel = ["E-Field (V/cm)","R (a0)"]
    conv = [1e5,BohrRadius.val]
    
    x_var = hcat([[sol_i.Ef, sol_i.R] for sol_i in sol_vec]...)
    x_varOI = argmax([maximum(abs.(diff(x_var[i, :])[1])) for i in 1:size(x_var, 1)])
    
    AMNodes = endNode(mol.M.basisTree)
    basisUC_Full = getBasisUC(mol.M.basisTree)
    basisUC = [getBasis(node.spin) for node in AMNodes]
    kronStateOI = [zeros(Int64, length(basis)) for basis in basisUC]
    for (ind, basisNode) in enumerate(basisUC)
        if stateOI[ind] == nothing
    
            kronStateOI[ind] .= 1
        else 
            subIndOI = findfirst(isequal(stateOI[ind]), basisNode)
            kronStateOI[ind][subIndOI] = 1
        end
    end
    
    stateOIind = findall(kron(kronStateOI...) .== 1)
    function findStateVec(sol_end::sol, stateOIs::Vector{Int64})
        outputVec = zeros(Int64, length(stateOIs))
    
        vec_end = sol_end.vec
        arrayOI = zeros(Float64, size(vec_end, 1), length(stateOIs))
        for col in 1:size(arrayOI, 2)
            arrayOI[stateOIs[col], col] = 1
        end
    
        B = zeros(ComplexF64, size(vec_end))
        for j in 1:size(vec_end, 2)
            max_index = argmax(abs.(vec_end[:, j]))
            B[max_index, j] = vec_end[max_index, j]
        end
        maxOverlapMat = B'*arrayOI
        for col in 1:size(maxOverlapMat, 2)
            outputVec[col] = argmax(abs.(maxOverlapMat[:, col]))
        end
        unique(outputVec)
    end
    
    adiabaticInd = argmax([sol_i.R for sol_i in sol_vec])
    filteredStateOI = findStateVec(sol_vec[adiabaticInd], stateOIind)
    filteredStateOI = length(filteredStateOI) > maxDisplay ? filteredStateOI[1:maxDisplay] : filteredStateOI
    
    #p = PlotlyJS.plot()
    for state in filteredStateOI
        PlotlyJS.addtraces!(p, PlotlyJS.scatter(x = x_var[x_varOI, :]/conv[x_varOI], y = [sol_i.val[state] for sol_i in sol_vec]*yConv[units], text=[KetName(sol_i.vec[:, state], basisUC_Full, QMorder = [6, 5, 4, 3, 2, 1]) for sol_i in sol_vec], name = "State $state"))
    end
    set_figure_style!(p, title = "Energy Spectrum", xlabel = xlabel[x_varOI], ylabel = "Energy ($(units))")
    p
end