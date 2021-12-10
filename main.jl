include("Geneticos.jl")
include("txtFormat.jl")

data=datos("Data1.txt")

mtz=data.dist

calif=x->costoRuta(x, mtz)

califFinal, genomaFinal = algoritmoGenetico(calif, "OX", 40, 1000; intStart=1, intEnd=floor(Int,length(data.startPoint)/2), random=5)

graficar(data.startPoint, genomaFinal, califFinal)

#califFinal, genomaFinal = algoritmoGenetico(calif, "PMX", 40, 10000; intStart=1, intEnd=131)
