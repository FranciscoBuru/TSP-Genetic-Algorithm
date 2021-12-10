tipo="OX"
numEval=40000
directorio= "Data1.txt"


data=datos(directorio)

mtz=data.dist

calif=x->costoRuta(x, mtz)
#Fijamos el numero de generaciones y optimizamos Pm
#Random lo dejamos en 0 y keepbest=true
n_exprimentos=101
Pm_bajo=0
Pm_alto=1

paso=(Pm_alto-Pm_bajo)/(n_exprimentos-1)
pob=41
random=0
gen=ceil(Int, numEval/(pob-1))


res=Array{Int64}(undef,gen, n_exprimentos)
i=1
best=1000000000
Pm_best=0
pms=Vector(Pm_bajo:paso:Pm_alto)
for pm in pms
    p, califFinal, genomaFinal = algoritmoGeneticoReporte(calif, tipo, pob, gen; intStart=1, intEnd=131, random=random, Pm=pm)
    res[:,i]=p
    if best>p[end]
        best=p[end]
        Pm_best=pm
        print("Pm:",pm, ", new best:",best,", ")
    end
    println("i:",i)
    i+=1

end
plot(1:gen, res, title="Variación en Pm en "*tipo, ylabel="Puntajes", xlabel="Generacón", labels="")

plot(pms, res[end,:], label="", title="Mejor resultado con variacion de pm", xlabel="Pm", ylabel="Resultado")
