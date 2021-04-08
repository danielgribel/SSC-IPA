# SSC-IPA
Source code of "Semi-Supervised Clustering with Inaccurate Pairwise Annotations" (Gribel, Gendreau and Vidal, 2021)

## Run

To run the SSC-IPA algorithm, open the Julia terminal try the following commands:

```
julia> include("Optimizer.jl")

julia> in = Input(seed, max_it, supervision_flag, prior)

julia> main("dataset", "must_graph", "cannot_graph", in)
```

### Example

```
julia> include("Optimizer.jl")

julia> in = Input(1234, 50, 1, 0.9)

julia> main("vertebral.data", "vertebral-must.link", "vertebral-cannot.link", in)
```

### Parameters of `Input`

`seed`: Numerical seed

`max_it`: Maximum number of iterations the algorithm will take.

`supervision_flag`: Determines if pairwise supervision is used (0: unsupervised algorithm, 1: semi-supervised algorithm).

`prior`: Prior estimation regarding the experts' accuracy (between 0 and 1; enter -1 for no priors)

### Parameters of the `main` function

`dataset`: Dataset file. **Important:** The dataset file must be within the /data folder inside the project. You must provide a file with the .data extension along with a labels file. The labels file must have the .label extension. Example: For a dataset named vertebral.data, you must provide the vertebral.label file in the same folder.

`must_graph`: _Must-link_ graph. **Important:** The _must-link_ graph file must be within the /data folder inside the project.

`cannot_graph`: _Cannot-link_ graph. **Important:** The _cannot-link_ graph file must be within the /data folder inside the project.

### Data format

Dataset files. The dataset file has n rows and d columns, where n is the number of data samples and d is the number of features. Each line contains the values of the d features of a data sample, where x_ij correspond to the j-th feature of the i-th sample of the data. Each feature value is separated by a single space, as depicted in the scheme below:

| x<sub>11</sub> | x<sub>12</sub> | x<sub>13</sub> | ... | x<sub>1d</sub> |
| x<sub>21</sub> | x<sub>22</sub> | x<sub>23</sub> | ... | x<sub>2d</sub> |
| .... | .... | .... | ... | .... |
| x<sub>n1</sub> | x<sub>n2</sub> | x<sub>n3</sub> | ... | x<sub>nd</sub> |
