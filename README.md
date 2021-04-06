# SSC-IPA
Source code of SSC-IPA, from "Semi-Supervised Clustering with Inaccurate Pairwise Annotations" (Gribel, Gendreau and Vidal, 2021)

## Run

To run the SSC-IPA algorithm, open the Julia terminal try the following commands:

```
julia> include("Optimizer.jl")

julia> in = Input(seed, max_it, supervision_flag, prior)

julia> main("dataset", "graph_prefix", in)
```

### Example

`julia> include("Optimizer.jl")`

`julia> in = Input(1234, 50, 1, 0.9)`

`julia> main("vertebral", "vertebral-graph", in)`
