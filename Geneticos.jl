using LinearAlgebra
using Random
using StatsBase
using Distributions
using Plots

abstract type Genoma end
abstract type GenomaEntero<:Genoma end
abstract type GenomaBinario<:Genoma end
abstract type GenomaEntero2Hijos<:GenomaEntero end
abstract type GenomaEntero1Hijo<:GenomaEntero end

function corregir!(x::Array{Int64}, y::Array{Int64})
    n=length(x)
    com=union(x,y)
    difx=setdiff(com, x)
    if isempty(difx)
        return
    end
    dify=setdiff(com, y)
    for i in 1:length(difx)
        x[findlast(isequal(dify[i]),x)]=difx[i]
        y[findlast(isequal(difx[i]),y)]=dify[i]
    end
end

function mutar!(x::Array{Int64}, Pm::Float64)
    n=length(x)
    n_mutaciones=rand(Binomial(n,Pm))
    for i in 1:n_mutaciones
        j=sample(Vector(1:n), 2; replace=false)
        x[j[1]],x[j[2]]=x[j[2]],x[j[1]]
    end
end

mutable struct PMX<:GenomaEntero2Hijos
    genoma::Array{Int64}
    calif::Int64

    function PMX(starter::Int64, ender::Int64, rng::MersenneTwister, funCalif::Function)
        genoma=shuffle(rng, Vector(starter:ender))
        calif=funCalif(genoma)
        new(genoma, calif)
    end
    function PMX(genoma::Array{Int64}, funCalif::Function)
        new(genoma, funCalif(genoma))
    end
end

function crossover(p1::PMX, p2::PMX, rng::MersenneTwister, funCalif::Function, Pm::Float64)
    n=length(p1.genoma)
    i=sample(Vector(1:n), 2; replace=false)
    x=copy(p1.genoma); y=copy(p2.genoma)
    if i[1]<i[2]
        i[1],i[2]=i[2],i[1]
    end
    y[1:i[1]],x[1:i[1]]=x[1:i[1]],y[1:i[1]]
    y[i[2]:n],x[i[2]:n]=x[i[2]:n],y[i[2]:n]
    #corregir las cadenas
    corregir!(x,y)
    mutar!(x, Pm); mutar!(y, Pm)
    return PMX(x,funCalif), PMX(y,funCalif)
end

mutable struct OX<:GenomaEntero2Hijos
    genoma::Array{Int64}
    calif::Int64
    function OX(starter::Int64, ender::Int64, rng::MersenneTwister, funCalif::Function)
        genoma=shuffle(rng, Vector(starter:ender))
        calif=funCalif(genoma)
        new(genoma, calif)
    end
    function OX(genoma::Array{Int64}, funCalif::Function)
        new(genoma, funCalif(genoma))
    end
end

function crossover(p1::OX, p2::OX, rng::MersenneTwister, funCalif::Function, Pm::Float64)
    n=length(p1.genoma)
    x=copy(p1.genoma); y=copy(p2.genoma)
    i=sample(Vector(1:n), 2; replace=false)
    if i[1]<i[2]
        j=[1:i[1];i[2]:n]
    else
        j=[1:i[2];i[1]:n]
    end

    #Extraemos los elementos que cambiamos
    xj=x[j]; yj=y[j]
    #Obtenemos su nuevo orden
    sort!(xj, by=i->findfirst(isequal(i),y))
    sort!(yj, by=i->findfirst(isequal(i),x))
    #Hacemos el cambio
    x[j]=xj
    y[j]=yj
    #Mutamos
    mutar!(x, Pm); mutar!(y, Pm)
    return OX(x,funCalif), OX(y,funCalif)
end


mutable struct Poblacion{T<:Genoma}
    pob::Array{T} #Poblacion en si
    rng::MersenneTwister #Seed
    funCalif::Function #Funcion calificadora, debe ser positiva entera y a maximizar
    keepbest::Bool #Numero de genomas que se guardan de la poblacion pasada
    random::Int64 #Numero de genomas que se inicializan al azar cada vuelta
    Pm::Float64 #Probabilidad de mutacion
    #Constructor entero
    function Poblacion{T}(n::Int64, starter::Int64, ending::Int64, funCalif::Function; keepbest::Bool=true, random::Int64=2, Pm::Float64=0.05, seed::Int64=1234) where T<:GenomaEntero
        rng=MersenneTwister(seed)
        pob=Array{T}(undef, n)
        for i in 1:n
            pob[i]=T(starter, ending, rng, funCalif)
        end
        new{T}(pob, rng, funCalif, keepbest, random, Pm)
    end
    #Constructor binario
    function Poblacion{T}(n::Int64, len::Int64, funCalif::Function; keepbest::Bool=true, random::Int64=2, Pm::Float64=0.05, seed::Int64=1234) where T<:GenomaBinario
        rng=MersenneTwister(seed)
        pob=Array{T}(undef, n)
        for i in 1:n
            pob[i]=T(len, funCalif)
        end
        new{T}(pob, rng, funCalif, keepbest, random, Pm)
    end
    #Constructor default
    function Poblacion{T}(pob::Array{T}, rng::MersenneTwister, funCalif::Function, keepbest::Bool, random::Int64, Pm::Float64) where T<:Genoma
        new(pob, rng, funCalif, keepbest, random, Pm)
    end
end

function rouletteselector(pob::Array{T}, n::Int64) where T<:Genoma
    f=x->x.calif; w=f.(pob)
    W=FrequencyWeights(-w.+(maximum(w)+1))
    sample(pob, W, n; replace=false)
end

function reproduce(p::Poblacion{T}) where T<:GenomaEntero2Hijos
    n=length(p.pob)
    pob=Array{T}(undef, n)
    m=n-p.keepbest-p.random
    for i in 1:ceil(Int,m/2)
        s=rouletteselector(p.pob, 2)
        c=crossover(s[1], s[2], p.rng, p.funCalif, p.Pm)
        #(p1::OX, p2::OX, rng::MersenneTwister, funCalif::Function, Pm::Float64)
        pob[2*i-1]=c[1]
        if 2*i<=m
            pob[2*i]=c[2]
        end
    end
    if p.random>0
        starter=minimum(p.pob[1].genoma)
        ender=maximum(p.pob[1].genoma)
        for i in 1:p.random
            pob[m+i]=T(starter, ender, p.rng, p.funCalif)
        end
    end
    if p.keepbest
        f=x->x.calif
        pob[n]=p.pob[findfirst(isequal(minimum(f,p.pob)),f.(p.pob))]
    end
    Poblacion{T}(pob, p.rng, p.funCalif, p.keepbest, p.random, p.Pm)
end
function getbest(pob::Poblacion)
    f=x->x.calif;
    return pob.pob[findfirst(isequal(minimum(f,pob.pob)),f.(pob.pob))]
end

function calif(pnt)
    n = length(pnt)
    #print(n)
    cost = 0;
    for p in 1:n-1
        cost += _mtz1[pnt[p], pnt[p+1]]
    end
    cost += _mtz1[pnt[1], pnt[n]]
end
#p=Poblacion{PMX}(10, 1, 6, calif; seed=1)


function algoritmoGenetico(funCalif::Function, tipo, pobsize, generations; Pm=0.01, intStart=0, intEnd=0, binlen=10, seed=1234, keepbest=true,random=2)
    if tipo=="PMX"
        pob=Poblacion{PMX}(pobsize, intStart, intEnd, funCalif)#; keepbest=keepbest, random=random, Pm=Pm, seed=seed)
        #Poblacion{T}(n::Int64, starter::Int64, ending::Int64, funCalif::Function; keepbest::Bool=true, random::Int64=2, Pm::Float64=0.05, seed::Int64=1234)
    elseif tipo=="OX"
        pob=Poblacion{OX}(pobsize, intStart, intEnd, funCalif; keepbest=keepbest, random=random, Pm=Pm, seed=seed)
    end

    for i in 2:generations
        println(i,getbest(pob))
        pob=reproduce(pob)
    end
    res=getbest(pob)
    return res.calif, res.genoma
end

## Lee Datos

data = datos("Data1.txt")

# AGUASSS!!!!!!!!!!!!!!!!!  Lee abajo

## Tienen que cambiar _mtz1 por _mtz2 en la funcion calif para poder correr el
# set de datos
# Julia no te deja redefinir globales entonces creamos una por dataset.

### ----------------------  Para data 1 ------------------------
global const _mtz1 = data[1]
### ----------------------  Para data 2 ------------------------
#global const _mtz2 = data[1]




califFinal, genomaFinal = algoritmoGenetico(calif, "OX", 40, 1000; intStart=1, intEnd=floor(Int,length(data[2])/2), random=5)
califFinal, genomaFinal = algoritmoGenetico(calif, "PMX", 40, 10000; intStart=1, intEnd=131)

## Graficamos problema y sol.

x = data[2][1:floor(Int,length(data[2])/2)]
y = data[2][floor(Int,length(data[2])/2+1):length(data[2])]

plot(x, y, seriestype = :scatter, title = "Ciudades, costo = "*string( califFinal), size=(1000,1000))

# Construimos ruta
orden = Array{Tuple{Int, Int}}(undef, length(x))
for num in 1:length(x)
    orden[num] = (data[2][genomaFinal[num],1] , data[2][genomaFinal[num],2])
end

# Ponemos Ruta
beam = Shape(orden)
plot!(beam, fillcolor = plot_color(:yellow, 0.3), fillalpha=0.0, alpha=0.2)
