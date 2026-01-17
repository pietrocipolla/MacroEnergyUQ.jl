function name_cuts!(EP_master::Model, counter::Int64)
    for con in all_constraints(EP_master,include_variable_in_set_constraints=false)
        if name(con) == ""
            set_name(con,"BendersCut"*string(counter))
        end
        counter+=1
    end 
    return counter
end

function setup_mga_master_problem!(EP_master::Model,setup::Dict)
 

    @constraint(EP_master,cMGABudget, EP_master[:eObj] + sum(EP_master[:vTHETA]) == setup["MGABudget"])

end

function make_rand_vecs(iterations::Int64, TechTypes::Int64, n_lines::Int64, Zones::Int64, ag::Bool)
    gen_vecs = rand(Float64,(TechTypes,Zones,iterations))
    if ag == true
        gen_vecs = rand(Float64,(TechTypes,iterations))
    end
        
    line_vecs = rand(Float64,(n_lines,iterations))
    return gen_vecs, line_vecs
end

function make_capMM_vecs(iterations::Int64, TechTypes::Int64, n_lines::Int64,Zones::Int64)
    cap_vecs =  rand(-1:1,TechTypes,Zones,2*iterations)
    cap_vecs = check_it_a(cap_vecs,iterations)
    line_vecs = rand(-1:1,n_lines,iterations)
    #check_it_a_ag!(vecs,iterations)
    #cap_vecs = convert_ag_to_disag(vecs,Zones)
    return cap_vecs, line_vecs
end

function unique_int(points::AbstractArray)
    pointst = transpose(points)
    nrow, ncol = size(points)

    uniques = fill(-2, (nrow, ncol))
    counter=0
    for i in 1:ncol
        for k in 1:ncol
            if points[:,i]==uniques[:,k]
                break
            elseif k == ncol
                counter = counter + 1
                uniques[:,counter] = points[:,i]
            end
        end
    end
    uniques = uniques[1:end, 1:counter]
    uniquesT = transpose(uniques)
    println("Done with uniques")
    return uniques
end

function make_combo_vecs(iterations::Int64, TechTypes::Int64,n_lines::Int64, Zones::Int64, ratio::Float64)
    rand_vecs, r_line_vecs = make_rand_vecs(ceil(Int64,iterations*ratio),TechTypes,n_lines,Zones)
    cap_vecs, c_line_vecs = make_capMM_vecs(floor(Int64,iterations*(1-ratio)),TechTypes,n_lines,Zones)
    
    gen_vecs = cat(rand_vecs,cap_vecs,dims=3)
    gen_vecs = vecs[:,:,1:iterations]
    
    line_vecs = cat(r_line_vecs, c_line_vecs, dims = 2)
    return gen_vecs, line_vecs
end

function convert_ag_to_disag(ag_vecs::AbstractArray, Zones::Int64)
    (techs,iterations) = size(ag_vecs)
    vecs = Array{Float64,3}(undef,(techs,Zones,iterations))
    for i in 1:iterations
        for j in 1:techs
			vecs[j,:,i] .= ag_vecs[j,i]
        end
    end
    return vecs
end

function check_it_a_ag!(a::AbstractArray, iterations::Int64)
    (r,i) = size(a)
    if iterations < i
        a = a[1:r,1:iterations]
        return a
    else
        println("Error")
    end
end

function check_it_a(a::AbstractArray, iterations::Int64)
    (r,c,i) = size(a)
    if iterations < i
        a = a[1:r,1:c, 1:iterations]
        return a
    else
        println("Error")
    end
end

function find_ratio(setup::Dict)
    ratio = 0.0
    if "ComboRatio" in keys(setup)
        ratio = setup["ComboRatio"]
        if ratio < 1
            return ratio
        else
            throw(ErrorException("Ratio greater than 1"))
        end
    else
        ratio = 0.25
    end
    return ratio
end

function generate_vecs(inputs::Dict, setup::Dict)
    iterations = setup["ModelingToGenerateAlternativeIterations"]
    TechTypes = collect(eachindex(unique(inputs["RESOURCES"].resource_type)))[end]
    n_lines = length(inputs["EXPANSION_LINES"])
    zones = inputs["Z"]
    method = setup["MGAMethod"]
    cluster_vecs = setup["ClusterMGAVecs"]

    
    n_its = iterations
    
    if method == 0
        ratio = find_ratio(setup)
        mats, line_vecs = make_combo_vecs(iterations,TechTypes,n_lines,zones,ratio)
    elseif method == 1
        mats, line_vecs = make_rand_vecs(iterations,TechTypes,n_lines,zones)
    elseif method == 2
        mats, line_vecs = make_capMM_vecs(iterations,TechTypes,n_lines,zones)
    end
    println(size(mats))
    max_mats = -1.0 .* mats
    max_line_vecs = -1.0 .* line_vecs
    all_mats = cat(mats, max_mats, dims=3)
    all_line_vecs = cat(line_vecs, max_line_vecs, dims=2)
    
    println(size(all_mats))
    if cluster_vecs == 1
        nclusters= setup["NumMGACluster"]
        focus_cluster = setup["FocusCluster"]
        if focus_cluster == 1
            iterations = 320
            nclusters = 16
        end
        all_vecs = Vector{Vector{Float64}}(undef,0)
        (r,c,its) = size(all_mats)
        vec_leng = r*c
        for i in 1:its
            mat = all_mats[:,:,i]
            vec = reshape(mat, vec_leng)
            push!(all_vecs, vec)
        end
        all_vecs = mapreduce(permutedims, vcat, all_vecs)
        if focus_cluster == 0
            all_vecs = kmeanscluster_vecs(all_vecs, nclusters)
        else
            all_vecs = kmeansfocuscluster_vecs(all_vecs, nclusters, n_its)
        end
        
        for i in 1:n_its*2
            all_mats[:,:,i] = reshape(all_vecs[i,:], (r,c))
        end
    end
    return all_mats, all_line_vecs
end

function kmeanscluster_vecs(vecs::AbstractArray, nclusters::Int64)
    vecsT = (vecs')
    result= kmeans(vecsT,nclusters)
    clusters=Vector{Vector{Vector{Float64}}}(undef,0)
    final_clusters = Vector{Vector{Float64}}(undef,0)
    for i in 1:nclusters
        push!(clusters, Vector{Vector{Float64}}(undef,0))
    end
    assignments = result.assignments
    for i in 1:length(assignments)
        push!(clusters[assignments[i]], vecs[i,:])
    end
    for i in 1:nclusters
        append!(final_clusters, clusters[i])
    end
    vecs_out = mapreduce(permutedims, vcat, final_clusters)
    return vecs_out
end

function kmeansfocuscluster_vecs(vecs::AbstractArray, nclusters::Int64, n_its::Int64)
    vecsT = (vecs')
    result= kmeans(vecsT,nclusters)
    clusters=Vector{Vector{Vector{Float64}}}(undef,0)
    final_clusters = Vector{Vector{Float64}}(undef,0)
    for i in 1:nclusters
        push!(clusters, Vector{Vector{Float64}}(undef,0))
    end
    assignments = result.assignments
    for i in 1:length(assignments)
        push!(clusters[assignments[i]], vecs[i,:])
    end
    for i in 1:nclusters
        if length(clusters[i]) >= n_its*2
            final_clusters = clusters[i][1:n_its*2]
            break
        end
        #append!(final_clusters, clusters[i])
    end
    vecs_out = mapreduce(permutedims, vcat, final_clusters)
    return vecs_out
end