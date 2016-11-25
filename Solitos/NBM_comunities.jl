#!/usr/bin/env julia

using LightGraphs
using RCall
using ArgParse

parser = ArgParseSettings(description = "Find communities using Non Backtracking Matrix algorithm.")

@add_arg_table parser begin
    "--output", "-o"
    "file"
end

args = parse_args(s)

output_file = get(args, "output", STDOUT)

red = readdlm(args["file"])
g = Graph()

ultimovertice = Int64(maximum(red))
add_vertices!(g,ultimovertice)
for n in 1:Int64((length(red)/2))
    add_edge!(g,Int64(red[n,1]),Int64(red[n,2]))
end

function la_matriz(g)
    A = full(adjacency_matrix(g))
    ceros = zeros(A)
    D = zeros(A)
    menos = -1 * eye(A)
    for n in 1:nv(g)
        D[n,n] = degree(g,n)-1
    end
    sparse(hcat(vcat(ceros,menos),vcat(D,A)))
end

M = la_matriz(g)
valores, vectores = eigs(M, nev=20)
treshold = real(sqrt(valores[1]))
if real(last(valores)) > treshold
    valores, vectores = eigs(M, nev=30)
    if real(last(valores)) > treshold
        valores, vectores = eigs(M, nev=40)
        if real(last(valores)) > treshold
            valores, vectores = eigs(M, nev=50)
            if real(last(valores)) > treshold
                valores, vectores = eigs(M, nev=nv(g))
            end
        end
    end
end


cuantos = Array(Float64,0)
index = Array(Int64,0)
for i in 1:length(valores)
    if imag(valores[i]) == 0 && real(valores[i]) > treshold
        push!(cuantos,real(valores[i]))
        push!(index,i)
    end
end

print("There are $(length(cuantos)) communities in this Network \n")

matriz_embedded = real(vectores[:,index])
arriba = matriz_embedded[1:nv(g),:]
abajo = matriz_embedded[nv(g)+1:2*nv(g),:]
hola = arriba + abajo

R_kmeans = R"kmeans";

function membresia(n,v)
    membresia = Array(Int64,0)
    if n == 2
        for i in 1:length(v[:,1])
            if sign(v[i,2]) == sign(1)
                push!(membresia,1)
            else
                push!(membresia,2)
            end
        end
    else
        matriz_sirve = v
        cluster_p = R_kmeans(matriz_sirve[:,2:n-1],n,nstart=500)
        for i in 1:length(cluster_p[1])
            push!(membresia,cluster_p[1][i])
        end
    end
    membresia
end

grupos = membresia(length(index),hola)

comunidades = hcat(vertices(g),grupos)


#print("$(length(cuantos))")

writedlm("../Resultados/NBM_communities_$(la_red)",comunidades)
end

print("Your results will be saved as: NBM_communities_$(a[k]) \n")


writedlm(output_file, comunidades)
