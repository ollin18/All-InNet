print("\n")
print("This script takes tab separeted .dat files as argument")
print("\n")
exten = r".dat"

using LightGraphs

function ext(v::Regex)
  valid = Array(Any,0)
  for i in readdir("../Redes")
    if match(v,i) != nothing
      push!(valid,i)
    end
  end
  valid
end

a = ext(exten)

graph_files = Dict{String,Int64}(Dict(a[i]=>i for i in 1:length(a)))
print("\n")
print("Give the index of the file that contains the network to analize \n \n")

for i in sort(collect(keys(graph_files)))
  println("$(graph_files[i]) => $i ")
end

print("\n")

k = readline(STDIN)
k = parse(Int,k)
print("You selected $(a[k])")
print("\n \n")

ad_list = a[k]


red = readdlm("../Redes/$ad_list",'\t')
if typeof(red[1,1]) == Float64 && typeof(red[1,2]) == Float64
  red = round(Int64,red)
end
if size(red)[2] != 2
  print("This is a preview of the file, select the columns that correspond to the nodes")
  print("\n")
  for n in 1:4
    println(red[n,:])
  end
  print("\n")
  print("press return after every selection")
  print("\n")
  i = readline(STDIN)
  i = parse(Int,i)
  j = readline(STDIN)
  j = parse(Int,j)
else
  i = 1
  j = 2
end

if typeof(red[1,i]) != Int
  Nodes = red[:,[i,j]]
  Nodes = unique(Nodes)
  #weight = red[:,2]
  dic_nodes = Dict{String,Int64}(Dict(Nodes[i]=>i for i in 1:length(Nodes)))
  g = Graph()
  last_node = Int64(length(Nodes))
  add_vertices!(g,last_node)
  for n in 1:Int64((length(red)/j))
      add_edge!(g,dic_nodes[red[n,i]],dic_nodes[red[n,j]])
  end
  else
  g = Graph()
  last_node = Int64(maximum(red))
  add_vertices!(g,last_node)
  for n in 1:Int64((length(red)/2))
      add_edge!(g,Int64(red[n,i]),Int64(red[n,j]))
  end
end

print("$(a[k]) is a {vertices,edges} = ",g, "\n \n")#, "wait", "\n \n")
# run(`cowsay Frankie says RELAX`)
print("\n")
if typeof(red[1,i]) != Int
  return g, Nodes, dic_nodes
else
  return g
end
