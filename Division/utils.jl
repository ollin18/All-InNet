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


function nonBM_embedding(g)
  NB = la_NBM(g)
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

function grados(g)
	grados = Array(Int64,0)
	gradosdiv = Array(Int64,0)
	for i in vertices(g)
		push!(grados,degree(g,i))
	end
	grados
end

function mat_an(g,v,NBM)
	treshold = sqrt(mean(v)/mean(v.-1)/mean(v))
	unos = ones(2ne(g))/sqrt(2ne(g))
	unnos = unos*unos'
	R = NBM - unnos
	R
end


function eigens(M)
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
	valores, vectores
end

function num_com(v)
	cuantos = Array(Float64,0)
	index = Array(Int64,0)
	for i in 1:length(v)
		if imag(v[i]) == 0 && real(v[i]) > treshold
		    push!(cuantos,real(v[i]))
		    push!(index,i)
		end
	end
	cuantos, index
end

function contrae(g,v,ei)
    y = zeros(Float64, nv(g))
    for i in 1:nv(g)
        for j in neighbors(g,i)
            u = ei[Edge(j,i)]
            y[i] += v[u]
            y[j] += v[u]
        end
    end
    y
end

function contraccion(g,vec,emb,ei)
	contraida = zeros(Float64, nv(g), length(vec))
	for n in 1:length(vec)
    	contraida[:,n] = contrae(g,emb,ei)
	end
	contraida
end

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

