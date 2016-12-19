#!/usr/bin/env julia

using LightGraphs
using RCall
using ArgParse

parser = ArgParseSettings(description = "Find communities using Flux matrix algorithm.")

@add_arg_table parser begin
    "--output", "-o"
    "file"
end

args = parse_args(parser)

output_file = get(args, "output", STDOUT)
if is(output_file, nothing)
  output_file = STDOUT
else
  output_file = args["output"]
end

function debug(msg)
  write(STDERR, msg)
end

debug("Output File: $(output_file)\n")

red = readdlm(args["file"])

g = Graph()

ultimovertice = Int64(maximum(red))
add_vertices!(g,ultimovertice)
for n in 1:Int64((length(red)/2))
    add_edge!(g,Int64(red[n,1]),Int64(red[n,2]))
end


function kron_δ(i,j)
  if i == j
    return 1
  else
    return 0
  end
end

"""
function la_NBM(g::SimpleGraph)
    edgeidmap = Dict{Edge, Int}()
    m = 0
    for e in edges(g)
        m += 1
        edgeidmap[e] = m
    end

    if !is_directed(g)
        for e in edges(g)
            m += 1
            edgeidmap[reverse(e)] = m
        end
    end


    B = spzeros(Float64, m, m)


    for (e,u) in edgeidmap
        i, j = src(e), dst(e)
        for (e2,v) in edgeidmap
          k, l = src(e2), dst(e2)
          B[v,u] = (kron_δ(j,k)*(1-kron_δ(i,l)))*(1/(degree(g,j)))
        end
    end
    return B, edgeidmap
end
"""
function la_NBM(g::SimpleGraph)
    edgeidmap = Dict{Edge, Int}()
    aristas = ne(g)
    m = 0
    for e in edges(g)
        m += 1
        edgeidmap[e] = m
	edgeidmap[reverse(e)] = m + aristas
    end

    B = spzeros(Float64, 2*aristas, 2*aristas)


    for (e,u) in edgeidmap
        i, j = src(e), dst(e)
	eles = neighbors(g,j)
	k = j
        for l in eles
          B[edgeidmap[k=>l],u] = (kron_δ(j,k)*(1-kron_δ(i,l)))*(1/(degree(g,j)))
        end
    end
    return B, edgeidmap
end

function contrae(g,v)
    y = zeros(Float64, nv(g))
    for i in 1:nv(g)
        for j in neighbors(g,i)
            u = NB2[Edge(j,i)]
            y[i] += v[u]
            y[j] += v[u]
        end
    end
    y
end

function nonBM_embedding(g)
  NB = la_NBM(g)[1]
  λ, eigve = eigs(NB, nev=20)
  los_reales = Array(Float64,0)
  los_indices_reales = Array(Int64,0)
  for i in 1:length(λ)
        if imag(λ[i]) == 0 && real(λ[i]) != 0 && i*real(λ[i]) >= sqrt(real(λ[1]))
          push!(los_reales,λ[i])
          push!(los_indices_reales,i)
      end
  end
  matriz_embedded = sub(real(eigve),:,los_indices_reales)
  length(los_reales), los_indices_reales
  ϕ = zeros(Float64, nv(g), length(los_reales))
  for n=1:(length(los_reales))
      v= matriz_embedded[:,n]
      ϕ[:,n] = contrae(g, v)
  end
  return length(los_reales), ϕ
end

grados = Array(Int64,0)
gradosdiv = Array(Int64,0)
for i in vertices(g)
    push!(grados,degree(g,i))
end

#treshold = sqrt(sum(grados/(grados-1))/nv(g))/(sum(grados)/nv(g))
treshold = sqrt(mean(grados)/mean(grados.-1)/mean(grados))
#treshold = sqrt(mean(grados./(grados.-1))/mean(grados))
NB1,NB2 = la_NBM(g)
unos = ones(2ne(g))/sqrt(2ne(g))
unnos = unos*unos'
R = NB1 - unnos

M = R
valores, vectores = eigs(M, nev=20)
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

print("There are $(length(cuantos)+1) communities in this Network \n")

matriz_embedded = real(vectores[:,index])

hola = zeros(Float64, nv(g), length(index))
for n in 1:length(index)
    hola[:,n] = contrae(g,matriz_embedded)
end

R_kmeans = R"kmeans"

function membresia(n,v)
    membresia = Array(Int64,0)
    if n == 2
        for i in 1:length(v[:,1])
            if sign(v[i,1]) == sign(1)
                push!(membresia,1)
            else
                push!(membresia,2)
            end
        end
    else
        matriz_sirve = v
        cluster_p = R_kmeans(matriz_sirve[:,1:n-1],n,nstart=500)
        for i in 1:length(cluster_p[1])
            push!(membresia,cluster_p[1][i])
        end
    end
    membresia
end

grupos = membresia(length(index)+1,hola)

if typeof(red[1,1]) != Int
  comunidades = hcat(Nodes,grupos)
else
  comunidades = hcat(vertices(g),grupos)
end

print("Your results will be saved as: Flux_communities_$(a[k]) \n")
writedlm("../Resultados/Flux_communities_$(a[k])",comunidades)

writedlm(output_file, comunidades)
