include("DataMulti.jl")
include("Model.jl")

mutable struct Solution
    # Total log-likelihood
    ll::Float64

    # GMM log-likelihood
    llgmm::Float64

    # SBM log-likelihood
    llsbm::Float64

    # Dataset
    data::DataMulti
    
    # Sample-cluster assignment -- class representation
    y::Array{Int}
    
    # Sample-cluster assignment -- binary indicator
    z::Array{Int}
    
    # Number of edges for a pair of clusters
    m::Array{Float64}
    
    # Sum of sample degrees on each cluster
    kappa::Array{Float64}
    
    # Mean point of clusters
    mu::Array{Float64}
    
    # Standard deviation
    sigma2::Array{Float64}
    
    # Cluster cardinality -- number of samples inside the cluster
    counter::Array{Int}
    
    # Distance to center (Squared Euclidean)
    dist::Array{Float64}
end

function init_variables(data)
    z = zeros(Int, data.n, data.k)
    counter = Array{Int, 1}(undef, data.k)
    m = Array{Float64, 3}(undef, data.L, data.k, data.k)
    kappa = Array{Float64, 2}(undef, data.L, data.k)
    dist = Array{Float64, 1}(undef, data.n)
    sigma2 = zeros(Float64, data.k)
    return z, counter, m, kappa, dist, sigma2
end

function Solution(data, y, centers)
    # Create variables
    z, counter, m, kappa, dist, sigma2 = init_variables(data)

    # Initialize cluster binary-indicator variable
    [ z[i, y[i]] = 1 for i = 1:data.n ]

    # Initialize cluster cardinality
    [ counter[c] = sum(z[:, c]) for c = 1:data.k ]
    
    # Initialize the number of edges for a pair of clusters
    [ m[l, :, :] = update_m(data, data.G[l], y, data.nodes[l]) for l = 1:data.L ]
    
    # Initialize the sum of degrees of clusters
    [ kappa[l, :] = update_kappa(data, data.degree[l,:], y, data.nodes[l]) for l = 1:data.L ]
    
    # Initialize the cluster mean-point
    mu = copy(centers)
    
    # Initialize the sample-cluster distance
    [ dist[i] = distance(data.X[i, :], mu[ y[i], :]) for i = 1:data.n ]
    
    # Calculate the standard deviation
    [ sigma2[ y[i] ] += dist[i] for i = 1:data.n ]
    [ sigma2[r] = sigma2[r]/(counter[r] * data.d) for r = 1:data.k ]

    # Calculate the GMM log-likelihood
    # llgmm = data.n * data.d * (.5 + log(sqrt(sigma2)))
    llgmm = ll_GMM(data, sigma2, counter)
    
    beta = calc_beta(data, counter)

    # Calculate SBM probability matrix
    w = get_omega_prior(data, m, counter, beta)

    llsbm = data.input.ANNOTATION*ll_SBM_fixed_prior(data, m, counter, w, beta)

    ll = data.input.GAMMA*llgmm + llsbm

    # Create the solution
    return Solution(ll, llgmm, llsbm, data, y, z, m, kappa, mu, sigma2, counter, dist)
end

# Update the log-likelihood values
function update_ll(sol, llgmm, llsbm)
    sol.llgmm = llgmm
    sol.llsbm = sol.data.input.ANNOTATION*llsbm
    sol.ll = sol.data.input.GAMMA*llgmm + sol.data.input.ANNOTATION*llsbm
end

# Update the GMM parameters
function update_gmm_param(sol, sigma2, src, tgt, mu_src, mu_tgt)
    sol.sigma2 = sigma2
    sol.mu[src, :] = mu_src
    sol.mu[tgt, :] = mu_tgt
end

# Update the SBM parameters
function update_sbm_param(sol, m)
    sol.m = m
end

# Update the solution assignment
function update_assignment(sol, i, src, tgt)
    sol.z[i, src] = 0
    sol.z[i, tgt] = 1
    sol.y[i] = tgt
    sol.counter[src] -= 1
    sol.counter[tgt] += 1
end

# Update the samples-cluster distances
function update_distance(sol, dist)
    sol.dist = dist
end

# Compute the centroid-index metric
function centroid_index(c, g, k)
    orphan1 = trues(k)
    ci1 = k
    for r = 1:k
        cmin = argmin( [ distance(c[r, :], g[s, :]) for s = 1:k ] )
        if orphan1[cmin] == true
            ci1 -= 1
        end
        orphan1[cmin] = false
    end

    orphan2 = trues(k)
    ci2 = k
    for r = 1:k
        cmin = argmin( [ distance(g[r, :], c[s, :]) for s = 1:k ] )
        if orphan2[cmin] == true
            ci2 -= 1
        end
        orphan2[cmin] = false
    end
    return max(ci1, ci2)
end

# Jeffreys-divergence for two gaussians
function jeffreys_divergence(mu1, mu2, sig1, sig2, d)
    return kl_divergence(mu1, mu2, sig1, sig2, d) + kl_divergence(mu2, mu1, sig2, sig1, d)
end

# KL-divergence for two gaussians
function kl_divergence(mu1, mu2, sig1, sig2, d)
    kl = 0.5* ( d*(log(sig2) - log(sig1)) - d + distance(mu1, mu2)/sig2 + d*(sig1/sig2) )
    return kl
end

function kl_matching(mu1, mu2, sig1, sig2, pi1, pi2, k, d)
    sum = 0.0
    for r = 1:k
        min_s = minimum( [ kl_divergence(mu1[r, :], mu2[s, :], sig1[r], sig2[s], d) + log(pi1[r]/pi2[s]) for s = 1:k ] )
        sum += pi1[r] * min_s
    end
    return sum
end
