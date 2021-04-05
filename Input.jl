struct Input
    SEED::Int
    MAX_IT::Int
    GAMMA::Float64
    ANNOTATION::Int
    BETA::Float64
end

function Input(seed, max_it, gamma, annotation, beta)
    return Input(seed, max_it, gamma, annotation, beta)
end