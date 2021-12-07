using DelimitedFiles
using LinearAlgebra
using Random
using StatsBase
using Distributions

struct datos
    dist::Array{Int64}
    startPoint::Array{Int64}

    function datos(name::String)
        dist, startPoint = readData(name);
    end

    function start(num)
        return shuffle(MersenneTwister(1234), Vector(1:num))
    end

    function readData(name)
        # Saco arreglo de posiciones
        #   ini = time_ns()
        file = open(name);
        lines = parse(Int64,readline(file, keep=false))
        startPoint = start(lines)
        pos  = zeros(Int64, lines, 2)
        i = 1
        while eof(file) != true
            x,y = split(readline(file, keep=false), " ")[2:3]
            x = parse(Int64, x)
            y = parse(Int64, y)
            pos[i,1] = x
            pos[i,2] = y
            i = i+1
        end
        # Creo mtz de distancias
        dist = zeros(Int64, lines, lines)
        for j in 1:lines
            for jj in 1:lines
                dist[j,jj] = round(sqrt((pos[j,1]-pos[jj,1])^2 + (pos[j,2]-pos[jj,2])^2))
            end
        end
        # print("File loaded in ")
        # print((time_ns() - ini )/ 1e9)
        # print(" sec")
        return dist, pos
    end
end


function costoRuta(pnt)
    n = length(pnt)
    cost = 0;
    for p in 1:n-1
        cost += mtz[pnt[p], pnt[p+1]]
    end
    cost += mtz[pnt[1], pnt[n]]
end

# Construimos clase con archivo de datos 1
#data = datos("Data1.txt")
# Sacamos la mtz de costos y la ruta inicial
#matrizDistancias, rutaInicial = data
#global const mtz = matrizDistancas
# Verificamos costos de ruta
#costoRuta(data[2])
#costoRuta(rutaInicial)
