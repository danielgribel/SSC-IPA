using Random
using LinearAlgebra
using Discreet
using StatsBase
using Clustering
using Hungarian
using Printf
using Distances
using BenchmarkTools

include("InputManager.jl")
include("Solution.jl")
include("Model.jl")
include("Hash.jl")

# Tolerance epsilon
const Tol = 1e-4

# function profiler()
#     to = TimerOutput()
#     begin
#     @timeit to "title1" begin
#     # code to profile 1
#     end
#     @timeit to "title2" begin
#     # code to profile 2
#     end
#     reset_timer!(to)
#     show(to; allocations = false)
#     println("")
#     end
# end

function assign_rand_center(y1, y2, matching, data, Mt)
    # Resulting assignment
    y = Array{Int, 1}(undef, data.n)

    # Make the group assignment
    for i = 1:data.n
        if y1[i] == y2[i]
            y[i] = y1[i]
        else
            rdm = bitrand(Mt, 1)[1]
            y[i] = rdm*y1[i] + (1 - rdm)*matching[ y2[i] ]
        end
    end
    return y
end

function assign_closest_center(c, data)
    k = size(c)[1]
    y = Array{Int, 1}(undef, data.n)
    dist = Array{Float64, 1}(undef, data.n)
    for i = 1:data.n
        max_dist = Inf
        for r = 1:k
            dist_r = distance(data.X[i, :], c[r, :])
            if dist_r < max_dist
                dist[i] = dist_r
                y[i] = r
                max_dist = dist_r
            end
        end
    end
    return y, dist
end

function is_degenerated(y, k)
    return length(unique(y)) < k
end

function max_k(a, k)
    return partialsortperm(a, 1:k, rev = true)
end

function repair_degeneracy!(y, mu, dist, data)
    empty_clusters = setdiff(1:data.k, unique(y))
    nb_empty = length(empty_clusters) 
    most_distant = max_k(dist, nb_empty)
    for c = 1:nb_empty
        mu[empty_clusters[c], :] = data.X[most_distant[c], :]
    end
    y, dist = assign_closest_center(mu, data)
    return y, mu, dist
end

function crossover(data, mu1, mu2, y1, y2, Mt)
    # Number of clusters
    k = data.k
    
    # Cost matrix in the bipartite graph
    cost = Array{Float64, 2}(undef, k, k)
    [ cost[r, s] = distance(mu1[r, :], mu2[s, :]) for r = 1:k for s = 1:k ]

    # Solve the assignment problem
    matching1 = hungarian(cost)[1]

    # Random selector
    coins = bitrand(Mt, k)
    
    # Create the new centers coordinates
    coord = Array{Float64, 2}(undef, k, data.d)
    [ coord[r, :] = coins[r]*mu1[r, :] + (1 - coins[r])*mu2[ matching1[r], :] for r = 1:k ]
    
    # Assign samples to the closest center
    y, dist = assign_closest_center(coord, data)
    
    while is_degenerated(y, k)
        y, coord, dist = repair_degeneracy!(y, coord, dist, data)
    end

    return y, coord, dist
end

function eval_mu(sol, src, tgt)
    # Number of samples in source cluster
    elem_src = findall(sol.y .== src)

    # Number of samples in target cluster
    elem_tgt = findall(sol.y .== tgt)
    
    # Source center after relocation
    mu_s = sum(sol.data.X[i, :] for i in elem_src)/length(elem_src)
    
    # Target center after relocation
    mu_t = sum(sol.data.X[i, :] for i in elem_tgt)/length(elem_tgt)
    
    return mu_s, mu_t
end

function eval_dist(sol, src, tgt, mu_s, mu_t)
    # Get the samples in the source cluster: before relocation
    samples_src = findall(sol.y .== src)
        
    # Get the samples in the target cluster: before relocation
    samples_tgt = findall(sol.y .== tgt)

    # Sum of distances after relocation
    dist = copy(sol.dist)
    [ dist[e] = distance(sol.data.X[e, :], mu_s) for e in samples_src ]
    [ dist[e] = distance(sol.data.X[e, :], mu_t) for e in samples_tgt ]

    return dist
end

# Evaluate relocation: relocate sample if new likelihood is better than current likelihood
function eval_relocate(sol, p, tgt)
    # Source cluster
    src = sol.y[p]

    # Copy current SBM parameters
    m_ = copy(sol.m)

    # Current log-likelihood of the SBM model
    sbm_cost = sol.data.input.ANNOTATION*sol.llsbm

    # Get the updated values of SBM parameters -- usually cheap to obtain
    if sum(sol.data.degree[:, p]) > 0
        m_ = estimate_SBM(sol, m_, p, src, tgt)
        # m_ = @btime estimate_SBM($sol, $m_, $p, $src, $tgt)
        sol.counter[src] -= 1
        sol.counter[tgt] += 1
        beta = calc_beta(sol.data, sol.counter)
        w = get_omega_prior(sol.data, m_, sol.counter, beta)
        sbm_cost = sol.data.input.ANNOTATION*ll_SBM_fixed_prior(sol.data, m_, sol.counter, w, beta)
        # w = @btime get_omega_prior($sol.data, $m_, $sol.counter, $beta)
        # sbm_cost = @btime $sol.data.input.ANNOTATION*ll_SBM_fixed_prior($sol.data, $m_, $sol.counter, $w, $beta)
        sol.counter[src] += 1
        sol.counter[tgt] -= 1
    end

    dist_tgt = 0.0
    dist_src = 0.0
    
    # Get the updated values of sigma, if the source or target variance is zero (cluster size <= 1)
    if sol.sigma2[src] == 0 || sol.sigma2[tgt] == 0
        sig = zeros(Float64, sol.data.k)
        [ sig[ sol.y[i] ] += sol.dist[i] for i = 1:sol.data.n ]
        sig[src] -= sol.dist[p]
        distance_tgt = distance(sol.data.X[p, :], sol.mu[tgt, :])
        sig[tgt] += distance_tgt
        sig[src] = sig[src]/((sol.counter[src] - 1) * sol.data.d)
        sig[tgt] = sig[tgt]/((sol.counter[tgt] + 1) * sol.data.d)
        dist_src = sol.dist[p]/(2.0*sig[src]) + sol.data.d*log(sqrt(sig[src]))
        dist_tgt = distance_tgt/(2.0*sig[tgt]) + sol.data.d*log(sqrt(sig[tgt]))
    else
        dist_src = sol.dist[p]/(2.0*sol.sigma2[src]) + sol.data.d*log(sqrt(sol.sigma2[src]))
        dist_tgt = distance(sol.data.X[p, :], sol.mu[tgt, :])/(2.0*sol.sigma2[tgt]) + sol.data.d*log(sqrt(sol.sigma2[tgt]))
    end

    if (dist_tgt + sbm_cost) < (dist_src + sol.data.input.ANNOTATION*sol.llsbm)

        # @time begin
        # Update samples membership
        update_assignment(sol, p, src, tgt)

        # Update centers
        sol.mu[src, :], sol.mu[tgt, :] = eval_mu(sol, src, tgt)
        
        # Update distances
        sol.dist = eval_dist(sol, src, tgt, sol.mu[src, :], sol.mu[tgt, :])
        
        # Update variances
        [ sol.sigma2[r] = 0.0 for r = 1:sol.data.k ]
        [ sol.sigma2[ sol.y[i] ] += sol.dist[i] for i = 1:sol.data.n ]
        [ sol.sigma2[r] = sol.sigma2[r]/(sol.counter[r] * sol.data.d) for r = 1:sol.data.k ]
        
        # Update SBM parameters
        update_sbm_param(sol, m_)

        # Log-likelihood of the GMM model
        gmm_cost = ll_GMM(sol.data, sol.sigma2, sol.counter)

        # Update likelihood
        update_ll(sol, gmm_cost, sbm_cost)
        # end

        return true
    
    end
    return false
end

function estimate_SBM(sol, m, p, src, tgt)
    # Update the number of edges between groups and the sum of degrees
    for l = 1:sol.data.L
        for v in neighbors(sol.data.G[l], p)
            e = sol.data.G[l].weights[p, v]
            if p == v
                m[l, src, src] -= e
                m[l, tgt, tgt] += e
            else
                m[l, sol.y[v], src] -= e
                m[l, src, sol.y[v]] -= e
                m[l, tgt, sol.y[v]] += e
                m[l, sol.y[v], tgt] += e
            end
        end
    end
    return m
end

function localsearch(sol, Mt)
    if sol.data.input.ANNOTATION == 0
        # K-means
        assign_unsupervised(sol)
    else
        has_relocated = true
        while has_relocated
            has_relocated = relocate(sol, Mt)
            if length(sol.data.unannotated) > 0
                assign(sol)
            end
        end
    end
end

function relocate(sol, Mt)
    prev_ll, curr_ll = Inf, sol.ll
    it_changed, has_relocated = true, false
    while (prev_ll - curr_ll) > Tol && it_changed
        it_changed = false
        samples = randperm(Mt, length(sol.data.annotated))
        for i1 in samples
            i = sol.data.annotated[i1]
            prev_c = sol.y[i]
            clusters = randperm(Mt, sol.data.k)
            for c in clusters
                if prev_c != c && sol.counter[ sol.y[i] ] > 1
                    if eval_relocate(sol, i, c)
                        prev_ll, curr_ll = curr_ll, sol.ll
                        it_changed, has_relocated = true, true
                    end
                end
            end
        end
    end
    return has_relocated
end

function closest_center(sol, i)
    return argmin( [ distance(sol.data.X[i, :], sol.mu[r, :])/(2.0*sol.sigma2[r]) + sol.data.d*log(sqrt(sol.sigma2[r])) for r = 1:sol.data.k ] )
end

function assign(sol)
    it_changed = true
    while it_changed
        it_changed = false
        for i in sol.data.unannotated
            src = sol.y[i]
            tgt = closest_center(sol, i)
            if tgt != src && sol.counter[src] > 2
                it_changed = true
                update_assignment(sol, i, src, tgt)
            end
        end
        sol.mu = update_mu(sol.data, sol.y)
        [ sol.dist[i] = distance(sol.data.X[i, :], sol.mu[sol.y[i], :]) for i = 1:sol.data.n ]
        [ sol.sigma2[r] = 0.0 for r = 1:sol.data.k ]
        [ sol.sigma2[ sol.y[i] ] += sol.dist[i] for i = 1:sol.data.n ]
        [ sol.sigma2[r] = sol.sigma2[r]/(sol.counter[r] * sol.data.d) for r = 1:sol.data.k ]
    end
    beta = calc_beta(sol.data, sol.counter)
    w = get_omega_prior(sol.data, sol.m, sol.counter, beta)
    sbm_cost = sol.data.input.ANNOTATION*ll_SBM_fixed_prior(sol.data, sol.m, sol.counter, w, beta)
    gmm_cost = ll_GMM(sol.data, sol.sigma2, sol.counter)
    update_ll(sol, gmm_cost, sbm_cost)
end

function assign_unsupervised(sol)
    it_changed = true
    while it_changed
        it_changed = false
        for i = 1:sol.data.n
            src = sol.y[i]
            tgt = closest_center(sol, i)
            if tgt != src && sol.counter[src] > 2
                it_changed = true
                update_assignment(sol, i, src, tgt)
            end
        end
        sol.mu = update_mu(sol.data, sol.y)
        [ sol.dist[i] = distance(sol.data.X[i, :], sol.mu[sol.y[i], :]) for i = 1:sol.data.n ]
        [ sol.sigma2[r] = 0.0 for r = 1:sol.data.k ]
        [ sol.sigma2[ sol.y[i] ] += sol.dist[i] for i = 1:sol.data.n ]
        [ sol.sigma2[r] = sol.sigma2[r]/(sol.counter[r] * sol.data.d) for r = 1:sol.data.k ]
    end
    gmm_cost = ll_GMM(sol.data, sol.sigma2, sol.counter)
    update_ll(sol, gmm_cost, 0.0)
end

function initial_assignment(data, Mt)
    # Random permutation
    rdm_order = randperm(Mt, data.n)
    # Create equaly-sized clusters from the random permutation 
    y = Array{Int, 1}(undef, data.n)
    [ y[ rdm_order[i] ] = ceil(data.k*i/data.n) for i = 1:data.n ]
    return y
end

function initial_population(data, pi1, Mt)
    population = Solution[]
    for i = 1:pi1
        # Obtain a solution with the k-means algorithm
        km_sol = kmeans(transpose(data.X), data.k; maxiter = 200, init = :kmpp, display = :none)
        @assert nclusters(km_sol) == data.k
        # Create a new solution with the assignment given by k-means
        sol = Solution(data, assignments(km_sol), transpose(km_sol.centers))
        # y = initial_assignment(data, Mt)
        # c = update_mu(data, y)
        # sol = Solution(data, y, c)
        # Apply the local search
        localsearch(sol, Mt)
        # Add solution to the population
        push!(population, sol)
    end
    # Sort population by likelihood
    sort!(population, by = v -> v.ll)
    return population
end

# Select two distinct solutions
function select_solutions(n)
    p1 = p2 = 1
    while p1 == p2
        p1 = rand(1:n)
        p2 = rand(1:n)
    end
    return p1, p2
end

function select_center(dist, data)
    c = rand(1:data.n)
    # items = Array((1:data.n))
    # c = sample(items, Weights(dist))
    return c
end

function remove_center!(y, mu, dist, c, data)
    active_centers = Array(1:data.k)[1:end .!= c]
    samples_c = findall(y .== c)
    for e in samples_c
        mindist = Inf
        for r in active_centers
            dist_e = distance(data.X[e, :], mu[r, :])
            if dist_e < mindist
                dist[e] = dist_e
                y[e] = r
            end
        end
    end
    return y, dist
end

function reinsert_center!(y, mu, dist, c, data)
    for i = 1:data.n
        dist_c = distance(data.X[i, :], mu[c, :])
        if dist_c < dist[i]
            y[i] = c
            dist[i] = dist_c
        end
    end
    return y, dist
end

function mutate(y, coord, dist, data, Mt)
    c = rand(1:data.k)
    y, dist = remove_center!(y, coord, dist, c, data)
    p = select_center(dist, data)
    coord[c, :] = data.X[p, :]
    y, dist = reinsert_center!(y, coord, dist, c, data)

    while is_degenerated(y, data.k)
        y, coord, dist = repair_degeneracy!(y, coord, dist, data)
    end
    
    mm = Solution(data, y, coord)
    return mm
end

function general_loop!(population, pi1, pi2, data, Mt, label)
    for it = 1:data.input.MAX_IT
        for i = 1:(pi2 - pi1)
            # Select two parent solutions
            p1, p2 = select_solutions(length(population))
            # Aplly the crossover operator
            y, coord, dist = crossover(data, population[p1].mu, population[p2].mu, population[p1].y, population[p2].y, Mt)
            # Mutate the crossover solution
            mm = mutate(y, coord, dist, data, Mt)
            # Apply the local search
            localsearch(mm, Mt)
            # Add solution to the population
            push!(population, mm)
        end
        # Select solutions for the next generation
        population = select_survivors(data, population, pi1, pi2)
    end
    sort!(population, by = v -> v.ll)
    t = 1
    return population[1:t]
end

function create_output_file(data, seed)
    OUTPUT_FILE = "out/"  * data.instance * "-" * string(seed) * "-" * string(data.input.GAMMA) * ".txt"
    return OUTPUT_FILE
end

function perc_edges_violations(y, data)
    must_viol = 0
    cann_viol = 0
    for e in collect(edges(data.G[1]))
        if y[e.src] != y[e.dst]
            must_viol += e.weight
        end
    end
    for e in collect(edges(data.G[2]))
        if y[e.src] == y[e.dst]
            cann_viol += e.weight
        end
    end

    total_edges = 0
    [ total_edges = total_edges + data.nb_edges[r] for r = 1:data.L ]
    total_viol = (must_viol + cann_viol)/total_edges

    return total_viol
end

function optimize(pi1, pi2, data, label, Mt)
    # Get the seed
    seed = signed(Mt.seed[1])

    # Create the output file
    OUTPUT_FILE = create_output_file(data, seed)

    # Start to measure CPU time
    t1 = time_ns()
    
    # Create the set of initial solutions
    population = initial_population(data, pi1, Mt)

    # Optimize with the general GA loop
    solutions = general_loop!(population, pi1, pi2, data, Mt, label)
    
    # Elapsed time of the algorithm
    cputime = (time_ns() - t1)/1.0e9

    # Ground-truth centroids
    g = update_mu(data, label)
    sol_truth = Solution(data, label, g)
    
    print_result(solutions[1], sol_truth, cputime, OUTPUT_FILE)
end

function print_result(sol, truth, cputime, OUTPUT_FILE)
    # Centroid index
    ci = centroid_index(sol.mu, truth.mu, sol.data.k)

    # Normalized Mutual Information
    nmi = mutual_information(sol.y, truth.y; normalize = true)

    # Percentage of violated annotations
    total_viol = perc_edges_violations(sol.y, sol.data)

    mu1 = sol.mu
    mu2 = truth.mu

    sig1 = sol.sigma2
    sig2 = truth.sigma2

    pi1 = zeros(Float64, sol.data.k)
    pi2 = zeros(Float64, sol.data.k)

    for r = 1:sol.data.k
        pi1[r] = sol.counter[r]/sol.data.n
        pi2[r] = truth.counter[r]/sol.data.n
    end

    # KL divergence
    kl = kl_matching(mu1, mu2, sig1, sig2, pi1, pi2, sol.data.k, sol.data.d)

    line = sol.data.instance * " "
    line *= string(sol.data.k) * " "
    line *= string(sol.data.input.ANNOTATION) * " "
    line *= string(sol.data.input.BETA) * " "
    line *= @sprintf("%.4f", truth.ll) * " "
    line *= @sprintf("%.4f", sol.ll) * " "
    line *= @sprintf("%.4f", sol.llgmm) * " "
    line *= @sprintf("%.4f", sol.llsbm) * " "
    line *= @sprintf("%.4f", nmi) * " "
    line *= string(ci) * " "
    line *= @sprintf("%.4f", kl) * " "
    line *= @sprintf("%.4f", total_viol) * " "
    line *= @sprintf("%.4f", cputime) * " "
    println(line)
    println(sol.y)
    # line *= string(sol.y) * " "
    # line *= "\n"
    # write_output(line, OUTPUT_FILE)
end

function write_output(line, OUTPUT_FILE)
    io = open(OUTPUT_FILE, "a")
    write(io, line)
    close(io)
end

function print_omega(w)
    k = size(w)[1]
    x = ""
    for r = 1:k
        for s = 1:k
            x *= @sprintf("%.4f", w[r, s]) * " "
        end
        x *= @sprintf("\n")
    end
    println(x)
end

function main(dataset, graph, input)
    # Random number generators
    Mt = MersenneTwister(input.SEED)
    Random.seed!(input.SEED)

    # Create an instance of data
    data, label = load_data(dataset, graph, input)

    # Initial population size
    pi1 = 10

    # Maximum population size
    pi2 = 20

    optimize(pi1, pi2, data, label, Mt)
end