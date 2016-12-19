using LightGraphs
using ArgParse

parser = ArgParseSettings(description = "Find communities using $(algorithm) algorithm.")

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
