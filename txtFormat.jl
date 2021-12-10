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
        new(dist, startPoint)
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


function costoRuta(pnt, mtz)
    n = length(pnt)
    cost = 0;
    for p in 1:n-1
        cost += mtz[pnt[p], pnt[p+1]]
    end
    cost += mtz[pnt[1], pnt[n]]
end

function graficar(startPoint::Array{Int64}, solucion::Array{Int64}, califFinal::Int64)
    x = startPoint[1:floor(Int,length(startPoint)/2)]
    y = startPoint[floor(Int,length(startPoint)/2+1):length(startPoint)]

    plot(x, y, seriestype = :scatter, title = "Ciudades, costo = "*string( califFinal), size=(1000,1000))

    # Construimos ruta
    orden = Array{Tuple{Int, Int}}(undef, length(x))
    for num in 1:length(x)
        orden[num] = (startPoint[solucion[num],1] , startPoint[solucion[num],2])
    end

    # Ponemos Ruta
    beam = Shape(orden)
    plot!(beam, fillcolor = plot_color(:yellow, 0.3), fillalpha=0.0, alpha=0.2)
end
