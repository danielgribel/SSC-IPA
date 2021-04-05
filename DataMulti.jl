using SimpleWeightedGraphs

include("Input.jl")

struct DataMulti
    # Attribute dataset
    X::Array{Float64}
    
    # Number of samples
    n::Int

    # Number of dimensions (attributes)
    d::Int

    # Number of clusters
    k::Int

    # Number of graphs (layers)
    L::Int

    # Graphs
    G::Array{SimpleWeightedGraph}

    # Samples degree
    degree::Array{Float64}

    # List of nodes in each graph with degree > 0
    nodes::Array{Any,1}

    # List of nodes with at least 1 connection
    annotated::Array{Int}

    # List of nodes with no connection
    unannotated::Array{Int}

    # Number of edges in each graph
    nb_edges::Array{Float64}

    # SBM constants
    C::Array{Float64}

    # Dataset name
    instance::String

    input::Input

    beta_must::Float64

    beta_cannot::Float64
end

function DataMulti(X, G, k, instance, input)
    n = size(X)[1]
    d = size(X)[2]
    L = length(G)

    degree = Array{Float64, 2}(undef, L, n)
    [ degree[l, :] = [ sum(G[l].weights[i,:]) for i = 1:n ] for l = 1:L ]

    nb_edges = Array{Float64, 1}(undef, L)
    [ nb_edges[l] = 0.5*sum(degree[l, :]) for l = 1:L ]

    C = Array{Float64, 1}(undef, L)
    [ C[l] = nb_edges[l]*(log(2.0*nb_edges[l]) - 1) for l = 1:L ]

    nodes = []
    [ push!(nodes, collect(findall(i -> i > 0, degree[l, :]))) for l = 1:L ]

    annotated = findall([ sum(degree[:, i]) for i = 1:n ] .> 0)
    
    unannotated = setdiff(1:n, annotated)

    beta_must   = 0.0
    beta_cannot = 0.0

    return DataMulti(X, n, d, k, L, G, degree, nodes, annotated, unannotated, nb_edges, C, instance, input, beta_must, beta_cannot)
end

