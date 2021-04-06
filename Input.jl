struct Input
    SEED::Int
    MAX_IT::Int
    ANNOTATION::Int
    BETA::Float64
end

function Input(seed, max_it, annotation, beta)
    return Input(seed, max_it, annotation, beta)
end