const HASH_SIZE = 97

mutable struct Hash
    cost::Float64
    card::Array{Int}
end

function Hash(cost, card)
    return Hash(cost, card)
end

function getkey(x)
    seed = 0
    for i = 1:length(x.card)
        seed = seed + (i * x.card[i])
    end
    return Int(floor(abs(x.cost) + seed) % HASH_SIZE)
end

function isequal(x, y)
    tol = 0.0000001
    if x.cost < (y.cost - tol) || x.cost > (y.cost + tol)
        return false
    end
    for i = 1:length(x.card)
        if(x.card[i] != y.card[i])
            return false
        end
    end
    return true
end

function do_collide(x, bucket)
    for y in bucket
        if isequal(x, y)
             return true
        end
    end
    return false
end

function item_found(x, bucket)
    if length(bucket) > 0
        if do_collide(x, bucket)
            return true
        end
    end
    return false
end

function select_survivors(data, population, pi1, pi2)
    hash_table = Array{Array{Hash}}(undef, HASH_SIZE)
    [ hash_table[i] = Hash[] for i = 1:HASH_SIZE ]

    clones = Solution[]
    individuals = Solution[]

    for i = 1:length(population)
        card = sort(population[i].counter)
        item = Hash(population[i].ll, card)
        key = getkey(item) + 1

        if item_found(item, hash_table[key])
            # Insert item in the heap of clones
            push!(clones, population[i])
        else
            # Insert item in the hash table
            push!(hash_table[key], item)
            # Insert item in the heap of individuals
            push!(individuals, population[i])
        end
    end

    sort!(individuals, by = v -> v.ll)
    sort!(clones, by = v -> v.ll)

    if length(individuals) >= pi1
        return individuals[1:pi1]
    end
    return vcat(individuals, clones[1:(pi1 - length(individuals))])
end