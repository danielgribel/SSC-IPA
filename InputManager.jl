using LightGraphs
using CSV
using DelimitedFiles

include("DataMulti.jl")

function load_data(dataset, graph_must, graph_cannot, input)
    DATA_FILE = "data/" * dataset
    INSTANCE = replace(dataset, ".data" => "")
    LABEL_FILE = "data/" * INSTANCE * ".label"
    
    # Load the dataset
    X = readdlm(DATA_FILE)

    # Load the labels file
    label = Int[]

    open(LABEL_FILE) do file
        [ push!(label, parse(Int, ln)) for ln in eachline(file) ]
    end

    # Number of clusters
    k = length(unique(collect(label)))

    # Number of layers (graphs)
    L = 2

    G = SimpleWeightedGraph[]

    # Load the graphs
    E = Matrix[]

    [ push!(G, SimpleWeightedGraph(size(X)[1])) for l = 1:L ]

    # [ push!(E, readdlm("data/" * graph * "-" *string(l)* ".link", Int)) for l = 1:L ]

    push!(E, readdlm("data/" * graph_must, Int))
    push!(E, readdlm("data/" * graph_cannot, Int))
    
    [ add_edge!(G[l], E[l][j,:][1], E[l][j,:][2], E[l][j,:][3]) for l = 1:L for j = 1:size(E[l])[1] ]

    # Create an instance of data
    data = DataMulti(X, G, k, INSTANCE, input)

    return data, label
end