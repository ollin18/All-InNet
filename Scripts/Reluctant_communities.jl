include("read_graph.jl")

run(`cowsay Frankie says RELAX`)

using RCall

function kron_δ(i,j)
    if i == j
        return 1
    else
        return 0
    end
end

function inv_degree(g,i)
    1/degree(g,i)
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
          B[edgeidmap[k=>l],u] = (kron_δ(j,k)*(1-kron_δ(l,i)))+(kron_δ(j,k)*kron_δ(l,i)*inv_degree(g,j))
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

treshold = sqrt(sum(grados/(grados-1))/nv(g))/(sum(grados)/nv(g))
NB1,NB2 = la_NBM(g)
unos = ones(2ne(g))/sqrt(2ne(g))
unnos = unos*unos'
R = NB1 - unnos

mmm = eigvals(R)

cuantos = Array(Float64,0)
index = Array(Int64,0)
for i in 1:length(mmm)
    if imag(mmm[i]) == 0 && real(mmm[i]) > i*treshold
        push!(cuantos,real(mmm[i]))
        push!(index,i)
    end
end

print("There are $(length(cuantos)+1) communities in this Network \n")

matriz_embedded = eigvecs(R)[:,index]

hola = zeros(Float64, nv(g), length(index))
for n in 1:length(index)
    hola[:,n] = contrae(g,matriz_embedded)
end

R_kmeans = R"kmeans";

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

print("Your results will be saved as: Reluctant_communities_$(a[k])\n")
writedlm("../Resultados/Reluctant_communities_$(a[k])",comunidades)
