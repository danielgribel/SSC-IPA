# Constrained-problem flag
const CONSTRAINED = true

function distance(x, y)
    return sqeuclidean(x, y)
end

function distance(x, y, d)
    s = 0.
    for i = 1:d
        s += (x[i] - y[i])^2
    end
    return s
end

# Update means mu, given y
function update_mu(data, y)
    mu = Array{Float64, 2}(undef, data.k, data.d)
    for r = 1:data.k
        elem_r = findall(y .== r)
        mu[r, :] = sum(data.X[i, :] for i in elem_r)/length(elem_r)
    end
    return mu
end

# Update matrix of the number of edges going from cluster r to s
function update_m(data, G, y, nodes)
    m = zeros(Float64, data.k, data.k)
    for i = 1:length(nodes)
        m[ y[ nodes[i] ], y[ nodes[i] ] ] += Float64(G.weights[ nodes[i], nodes[i] ])
        for j = (i+1):length(nodes)
            m[ y[ nodes[i] ], y[ nodes[j] ] ] += Float64(G.weights[ nodes[i], nodes[j] ])
            m[ y[ nodes[j] ], y[ nodes[i] ] ] += Float64(G.weights[ nodes[j], nodes[i] ])
        end
    end
    return m
end

# Update array of sum of degrees in each cluster
function update_kappa(data, degree, y, nodes)
    kappa = zeros(Float64, data.k)
    [ kappa[ y[ nodes[i] ] ] += degree[ nodes[i] ] for i = 1:length(nodes) ]
    return kappa
end

# Compute the GMM log-likelihood from scratch
ll_GMM(data, sigma2) = data.n * data.d * (.5 + log(sqrt(sigma2)))

function ll_GMM(data, sigma2, counter)
    ll = 0.
    for r = 1:data.k
        if sigma2[r] > 0
            ll += counter[r] * data.d * (0.5 + log(sqrt(sigma2[r])))
        end
    end
    return ll
end

# Compute the GMM log-likelihood from scratch
function ll_GMM(data, mu, y, sigma2)
    ll = 0.
    [ ll += distance(data.X[i, :], mu[y[i], :]) for i = 1:data.n ]
    ll = ll/(2.0*sigma2) + data.n * data.d * log(sqrt(sigma2))
    return ll
end

# DCSBM likelihood -- L \omegas, one for each layer 
function ll_DCSBM(data, m, kappa)
    if data.input.ANNOTATION == 0
        return 0.0
    end
    log_ll = 0.
    [ log_ll += ll_DCSBM(data, m, kappa, l) for l = 1:data.L ]
    return log_ll
end

function ll_DCSBM(data, m, kappa, l)
    log_ll = 0.
    for r = 1:data.k
        for s = r:data.k
            ratio = 0
            if m[l, r, s] > 0
                ratio = m[l, r, s]/(kappa[l, r] * kappa[l, s])
            end
            contrib = m[l, r, s] * log(max(1e-10, ratio))
            log_ll -= contrib
            log_ll -= (r != s) * contrib
        end
    end
    log_ll = 0.5*log_ll - data.C[l]
    return log_ll
end

function ll_DCSBM_fixed(data, m, kappa, w)
    if data.input.ANNOTATION == 0
        return 0.0
    end
    log_ll = 0.
    [ log_ll += ll_DCSBM_fixed(data, m, kappa, w, l) for l = 1:data.L ]
    return log_ll
end

function ll_DCSBM_fixed(data, m, kappa, w, l)
    log_ll = 0.
    for r = 1:data.k
        for s = r:data.k
            contrib = m[l, r, s] * log(max(1e-10, w[l, r, s]))
            contrib -= w[l, r, s] * (kappa[l, r] * kappa[l, s])/(2 * data.nb_edges[l])
            log_ll -= contrib
            log_ll -= (r != s) * contrib
        end
    end
    return 0.5*log_ll
end

function ll_SBM_fixed(data, m, counter, w)
    if data.input.ANNOTATION == 0
        return 0.0
    end
    log_ll = 0.
    [ log_ll += ll_SBM_fixed(data, m, counter, w, l) for l = 1:data.L ]
    return log_ll
end

function ll_SBM_fixed(data, m, counter, w, l)
    log_ll = 0.
    for r = 1:data.k
        for s = r:data.k
            contrib = m[l, r, s] * log(max(1e-10, w[l, r, s])) - w[l, r, s] * counter[r] * counter[s]
            log_ll -= contrib
            log_ll -= (r != s) * contrib
        end
    end
    return 0.5*log_ll
end

function calc_beta(data, counter)
    beta = zeros(Float64, data.L, data.k, data.k)
    sumin = 0.0
    sumout = 0.0
    p = data.input.BETA

    if p == -1.0
        return beta
    end

    for r = 1:data.k
        sumin += counter[r]*(counter[r] + 1)
        for s = (r+1):data.k
            sumout += counter[r]*counter[s]
        end
    end
    sumin = 0.5*sumin

    beta_must = 1/((data.nb_edges[1]*p)/(p*sumin + (1-p)*sumout))
    beta_cann = 1/((data.nb_edges[2]*p)/((1-p)*sumin + p*sumout))

    for r = 1:data.k
        beta[1, r, r] = beta_must
        beta[2, r, r] = (p/(1-p))*beta_cann
        for s = (r+1):data.k
            beta[1, r, s] = (p/(1-p))*beta_must
            beta[2, r, s] = beta_cann
            beta[1, s, r] = beta[1, r, s]
            beta[2, s, r] = beta[2, r, s]
        end
    end
    return beta
end

function ll_SBM_fixed_prior(data, m, counter, w, beta)
    if data.input.ANNOTATION == 0
        return 0.0
    end

    log_ll = 0.
    for l = 1:data.L
        for r = 1:data.k
            for s = r:data.k
                # Likelihood
                contrib_ll = m[l, r, s] * log(max(1e-10, w[l, r, s])) - w[l, r, s] * counter[r] * counter[s]
                log_ll += contrib_ll
                log_ll += (r != s) * contrib_ll
                # Priors
                log_ll += 2.0 * (data.input.BETA > 0) * (log(max(1e-10, beta[l, r, s])) - beta[l, r, s] * w[l, r, s])
            end
        end
    end
    return -0.5*log_ll
end

function get_degrees_term(data, kappa)
    k = data.k
    T = zeros(Float64, 2, k, k)
    [ T[1, r, s] = (kappa[1, r]*kappa[1, s])/(2*data.nb_edges[1]) for r = 1:k for s = r:k ]
    [ T[2, r, s] = (kappa[2, r]*kappa[2, s])/(2*data.nb_edges[2]) for r = 1:k for s = r:k ]
    [ T[1, s, r] = T[1, r, s] for r = 1:k for s = (r+1):k ]
    [ T[2, s, r] = T[2, r, s] for r = 1:k for s = (r+1):k ]
    return T
end

function max_off(w)
    k = size(w)[1]
    x = maximum( (w[r, s]) for r = 1:k for s = (r+1):k )
    return x
end

function min_off(w)
    k = size(w)[1]
    x = minimum( (w[r, s]) for r = 1:k for s = (r+1):k )
    return x
end

function get_omega(data, m, counter)
    omega = zeros(Float64, data.L, data.k, data.k)
    for l = 1:data.L
        for r = 1:data.k
            for s = r:data.k
                if m[l, r, s] == 0 || counter[r] == 0 || counter[s] == 0
                    omega[l, r, s] = 0.0
                else
                    numer = 1.0*m[l, r, s]
                    denom = counter[r] * counter[s]
                    omega[l, r, s] = numer/denom
                end
                omega[l, s, r] = omega[l, r, s]
            end
        end
    end
    return omega
end

function get_omega_prior(data, m, counter, beta)
    omega = zeros(Float64, data.L, data.k, data.k)

    for l = 1:data.L
        for r = 1:data.k
            for s = r:data.k
                if m[l, r, s] == 0 || counter[r] == 0 || counter[s] == 0
                    omega[l, r, s] = 0.0
                else
                    numer = 1.0*m[l, r, s]
                    denom = counter[r] * counter[s]
                    denom += (1 + (r == s)) * beta[l, r, s]
                    omega[l, r, s] = numer/denom
                end
                omega[l, s, r] = omega[l, r, s]
            end
        end
    end
    return omega
end