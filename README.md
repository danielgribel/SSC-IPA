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

```
julia> include("Optimizer.jl")

julia> in = Input(1234, 50, 1, 0.9)

julia> main("vertebral", "vertebral-graph", in)
```

### Parameters of `Input`

`seed`: Numerical seed

`max_it`: Maximum number of iterations the algorithm will take.

`supervision_flag`: Determines if pairwise supervision is used (0: unsupervised algorithm, 1: semi-supervised algorithm).

`prior`: Prior estimation regarding the experts' accuracy (between 0 and 1; enter -1 for no priors)
