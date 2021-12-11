include("Geneticos.jl")
include("txtFormat.jl")

datas=datos("Data1.txt")

Pm_OX=0.1
pob_OX=5
Pm_CX=0.1
pob_CX=5
Pm_PMX=0.1
pob_PMX=5

mtz=data.dist

calif=x->costoRuta(x, mtz)

tipos=["OX"; "CX"; "PMX"]
pms=[Pm_OX; Pm_CX; Pm_PMX]
pobs=[pob_OX; pob_CX; pob_PMX]

#generaciones=[1000; 10000; 100000; 1000000]
#generaciones=[10]
graficas=Array{Plots.Plot}(undef, 3)
graficas2=Array{Plots.Plot}(undef, 3)

for i in 1:3
    for gen in generaciones
        ini = time_ns()
        val, califFinal, genomaFinal, p = algoritmoGeneticoReporte(calif, tipos[i], pobs[i], gen; intStart=1, intEnd=size(mtz,1 ), random=0, Pm=pms[i])
        println(tipo, ", Tiempo:",(time_ns() - ini )/ 1e9, ", generaciones:",gen, ", calif:", califFinal)
    end
    graficas[i]=graficar(data.startPoint, genomaFinal, califFinal)
    graficas2[i]=p
end
