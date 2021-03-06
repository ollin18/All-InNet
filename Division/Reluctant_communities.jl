using LightGraphs
using RCall

R_kmeans = R"kmeans"

include("matrices.jl")
include("utils.jl")
include("read_graph.jl")

grad = grados(g)

NBM, edgeidmap = reluctant_matrix(g)

treshold = sqrt(mean(grad)/mean(grad.-1)/mean(grad))

R = mat_an(g,grad,NBM)

valores, vectores = eigens(R)

cuantos, index = num_com(valores)

matriz_embedded = real(vectores[:,index])

contraida = contraccion(g,index,matriz_embedded,edgeidmap)

grupos = membresia(length(index)+1,contraida)

comunidades = hcat(vertices(g),grupos)

println(comunidades)
