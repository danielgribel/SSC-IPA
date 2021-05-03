# SSC-IPA
Source code of "Semi-Supervised Clustering with Inaccurate Pairwise Annotations" (Gribel, Gendreau and Vidal, 2021)

<p align="center">
<img src="https://user-images.githubusercontent.com/4787247/116889763-efe70300-abfa-11eb-9810-5df99eaf65fc.png" width="400" height="225">
</p>

## Related Article

Semi-Supervised Clustering with Inaccurate Pairwise Annotations: https://arxiv.org/abs/2104.02146

## Run

To run the SSC-IPA algorithm, open the Julia terminal and try the following commands:

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

`dataset`: Dataset file. **Important:** You must provide a file with the `.data` extension along with a labels (ground-truth) file. The labels file must have the `.label` extension. Example: For a dataset named "vertebral.data", you must provide the "vertebral.label" file in the same folder.

`must_graph`: _Must-link_ graph file.

`cannot_graph`: _Cannot-link_ graph file.

**Important:** The dataset, labels, _must-link_ graph, and _cannot-link_ graph files must be within the `/data` folder inside the project.

### Data format

**Dataset files.** The dataset file has `N` rows and `D` columns, where `N` is the number of data samples and `D` is the number of features. Each line contains the values of the `D` features of a data sample, where x<sub>ij</sub> correspond to the j-th feature of the i-th sample of the data. Each feature value is separated by a single space, as depicted in the scheme below:

| x<sub>11</sub> | x<sub>12</sub> | x<sub>13</sub> | ... | x<sub>1d</sub> |
|------|------|------|-----|------|
| x<sub>21</sub> | x<sub>22</sub> | x<sub>23</sub> | ... | x<sub>2d</sub> |
| ... | ... | ... | ... | ... |
| x<sub>n1</sub> | x<sub>n2</sub> | x<sub>n3</sub> | ... | x<sub>nd</sub> |

**Important**: The dataset files must have the `.data` extension.

**Graph files.** A graph file (_must-link_ or _cannot-link_) has `m` rows and 3 columns, where `m` is the number of connections (links) in the graph. The first two columns represent the two data samples of an edge, whereas and third column represents the edge weight. The scheme below describes a graph file, where s<sub>i</sub> and t<sub>i</sub> are two connected samples, and w<sub>i</sub> is the corresponding edge weight:

| s<sub>1</sub> | t<sub>1</sub> | w<sub>1</sub> |
|------|------|------|
| s<sub>2</sub> | t<sub>2</sub> | w<sub>2</sub> |
| ... | ... | ... |
| s<sub>m</sub> | t<sub>m</sub> | w<sub>m</sub> |

**Labels files.** The content of a labels file exhibits the cluster of each sample of the dataset according to the ground-truth, where y<sub>i</sub> corresponds to the label of the i-th sample:

y<sub>1</sub>

y<sub>2</sub>

...

y<sub>n</sub>

**Important**: The labels files must have the `.label` extension.
