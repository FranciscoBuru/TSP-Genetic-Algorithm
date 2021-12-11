include("Geneticos.jl")
include("txtFormat.jl")

data=datos("Data1.txt")

mtz=data.dist

calif=x->costoRuta(x, mtz)

genes = 10000
println(genes)
ini = time_ns()
califFinal, genomaFinal = algoritmoGenetico(calif, "OX", 100, genes; intStart=1, intEnd=floor(Int,length(data.startPoint)/2), random=5)
print((time_ns() - ini )/ 1e9)
graficar(data.startPoint, genomaFinal, califFinal)

#califFinal, genomaFinal = algoritmoGenetico(calif, "PMX", 40, 10000; intStart=1, intEnd=131)



# ini = time_ns()

# print("File loaded in ")
# print((time_ns() - ini )/ 1e9)
# print(" sec")

p, califFinal, genomaFinal = algoritmoGeneticoReporte(calif, "Aleatorio", 40, 10000; intStart=1, intEnd=floor(Int,length(data.startPoint)/2), random=0)
println(califFinal)
plot!()
gen=Vector(1:10)
graficar(data.startPoint, genomaFinal, califFinal)
