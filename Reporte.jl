function reportePm(tipo, directorio, n_exprimentos, Pm_bajo, Pm_alto, numEval, pobsize, random)
    data=datos(directorio)
    mtz=data.dist
    calif=x->costoRuta(x, mtz)
    paso=(Pm_alto-Pm_bajo)/(n_exprimentos-1)
    n_gen=ceil(Int, numEval/(pobsize-1))
    res=Array{Int64}(undef,n_gen, n_exprimentos)
    i=1
    best=1000000000
    Pm_best=0
    pms=Vector(Pm_bajo:paso:Pm_alto)
    for pm in pms
        p, califFinal, genomaFinal = algoritmoGeneticoReporte(calif, tipo, pob, n_gen; intStart=1, intEnd=size(mtz, 1), random=random, Pm=pm)
        res[:,i]=p[:,1]
        if best>p[end,1]
            best=p[end,1]
            Pm_best=pm
            print("Pm:",pm, ", new best:",best,", ")
        end
        println("i:",i)
        i+=1
    end
    p1=plot(1:n_gen, res, title="Variación en Pm en "*tipo, ylabel="Puntajes", xlabel="Generacón", labels="")
    p2=plot(pms, res[end,:], label="", title="Mejor resultado con variacion de pm", xlabel="Pm", ylabel="Resultado")
    return res, p1, p2
end


tipo="OX"
numEval=100
directorio= "Data1.txt"
n_exprimentos=101
Pm_bajo=0
Pm_alto=1
pobsize=41
random=0


res=reportePm(tipo, directorio, n_exprimentos, Pm_bajo, Pm_alto, numEval, pobsize, random)
plot(res[2])
plot(res[3])

function reportePobsize(tipo, directorio, n_exprimentos, pob_bajo, pob_alto, numEval, Pm, random)
    data=datos(directorio)
    mtz=data.dist
    calif=x->costoRuta(x, mtz)
    paso=(pob_alto-pob_bajo)/(n_exprimentos-1)
    res=Array{Int64}(undef,n_exprimentos)
    i=1
    best=1000000000
    pob_best=0
    pobs=floor.(Int64, Vector(pob_bajo:paso:pob_alto))
    print(pobs)
    for pobsize in pobs
        n_gen=ceil(Int, numEval/(pobsize-1))
        p, califFinal, genomaFinal = algoritmoGeneticoReporte(calif, tipo, pobsize, n_gen; intStart=1, intEnd=size(mtz, 1), random=random, Pm=Pm)
        res[i]=p[end,1]
        if best>p[end,1]
            best=p[end,1]
            pob_best=pobsize
            print("pobsize:",pobsize, ", new best:",best,", ")
        end
        println("i:",i)
        i+=1
    end
    #p1=plot(1:n_gen, res, title="Variación en pobsize en "*tipo, ylabel="Puntajes", xlabel="Generacón", labels="")
    p=plot(pobs, res, label="", title="Mejor resultado con variacion de tamaño de población", xlabel="Tamaño de pblación", ylabel="Resultado")
    return res, p
end

tipo="OX"
numEval=40000
directorio= "Data1.txt"
n_exprimentos=99
pob_bajo=2
pob_alto=100
pobsize=41
random=0
Pm=0.1

res2 = reportePobsize(tipo, directorio, n_exprimentos, pob_bajo, pob_alto, numEval, Pm, random)
res2[2]
plot(res2[2])


function reporteHeatmap(tipo, directorio, n_exprimentos_pob, pob_bajo, pob_alto, numEval, Pm_bajo, Pm_alto, n_exprimentos_Pm, random)
    lk=ReentrantLock();
    data=datos(directorio)
    mtz=data.dist
    calif=x->costoRuta(x, mtz)
    res=Array{Int64}(undef,n_exprimentos_pob,n_exprimentos_Pm)
    i=1

    best=1000000000
    pob_best=0
    Pm_best=0

    paso_pob=(pob_alto-pob_bajo)/(n_exprimentos_pob-1)
    pobs=floor.(Int64, Vector(pob_bajo:paso_pob:pob_alto))

    paso_Pm=(Pm_alto-Pm_bajo)/(n_exprimentos_Pm-1)
    pms=Vector(Pm_bajo:paso_Pm:Pm_alto)
    k=1
    Threads.@threads for i in 1:n_exprimentos_pob
        pobsize=pobs[i]
        n_gen=ceil(Int, numEval/(pobsize-1))
        Threads.@threads for j in 1:n_exprimentos_Pm
            Pm=pms[j]
            p, califFinal, genomaFinal = algoritmoGeneticoReporte(calif, tipo, pobsize, n_gen; intStart=1, intEnd=size(mtz, 1), random=random, Pm=Pm)
            res[i,j]=p[end,1]
            lock(lk) do
                if best>res[i,j]
                    best=res[i,j]
                    print("pobsize:",pobsize, ", Pm:",Pm, ", new best:",best,", ")
                    if k%10!=0
                        println("i:",i,", j:",j,", k:",k)
                    end
                end
            end
            if k%10==0
                println("i:",i,", j:",j,", k:",k)
            end
            k+=1
            #
        end
    end
    p=heatmap(pobs, pms, res,
    c=cgrad([:blue, :white,:red, :yellow]),
        ylabel="Pm", xlabel="Tamaño de la poblacion",
        title="Modificacion de hiper parametros")
    return res, p
end


tipo="OX"
numEval=10000
directorio= "Data1.txt"
n_exprimentos_pob=20
pob_bajo=5
pob_alto=100
random=0
Pm_bajo=0.0
Pm_alto=0.1
n_exprimentos_Pm=11


res=reporteHeatmap(tipo, directorio, n_exprimentos_pob, pob_bajo, pob_alto, numEval, Pm_bajo, Pm_alto, n_exprimentos_Pm, random)
