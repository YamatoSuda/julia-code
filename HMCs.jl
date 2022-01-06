#HMCs.jl

#=
まず，配位を特徴づけるパラメータ群を配列に格納する，
もしくは，それぞれをセットして，もう1つ大きなパラメータベクトルに格納する．
具体的には，N=N+ + N- を全インスタントン数，Nfを系のクォークフレーバー数として次の自由度を持つ．

color orientation : SU(2)--> 3N
                  : SU(3)--> 7N
spacetime         : 4N
instanton size    : N
Quark Flavors     : Nf
=#

mutable struct ParamsofILM
    ρ::Array{Float64,1}     #インスタントンサイズ：１次元配列
    θ::Array{Float64,2}     #color orientationを決めるSU(N)のパラメータ：２次元配列
    z::Array{Float64,2}     #インスタントンの時空座標：２次元配列
    conjρ::Array{Float64,1}  #ρの共役変数：１次元配列
    conjθ::Array{Float64,2} #color orientationを決めるSU(N)のパラメータ：２次元配列
    conjz::Array{Float64,2} #インスタントンの時空座標：２次元配列
end

mutable struct ModelParams
    Nc::Int32  
    Nf::Int32  
    Λ::Float64 

    nI::Int64  
    nA::Int64  
    N::Int64   

    mag::Float64 
    A::Float64
    g::Float64

    gauge::String
end

mutable struct ColorWorkSpace
    u::Array{ComplexF64,1}
    sigma::Array{Any,1}
    lambda::Array{Any,1}
    tem::Array{ComplexF64,2}
    res::Array{ComplexF64,2}
end


function ChooseGaugeField!(params,nI,nA;gauge="SU(2)",flavor=0)
    if gauge == "SU(2)"
        params.Nc    = 3
        params.Nf    = flavor
        params.Λ     = 270

        params.nI    = nI
        params.nA    = nA
        params.N     = nI + nA

        params.mag   = 0.7 
        params.A     = 128
        params.g     = 1

        params.gauge = "SU(2)"
        println("----- Parameters Info. ---------------------------------------------")
        println("             Nc : 3")
        println("             Nf : $flavor")
        println("              Λ : 270 [MeV]")
        println("    Instanton # : $nI")
        println("AntiInstanton # : $nA")
        println("          Gauge : $gauge")
        println("--------------------------------------------------------------------")
    elseif gauge == "SU(3)"
        params.Nc    = 7
        params.Nf    = falvor
        params.Λ     = 270

        params.nI    = nI
        params.nA    = nA
        params.N     = nI + nA

        params.mag   = 0.7
        params.A     = 128
        params.g     = 1
        
        params.gauge = "SU(3)"
        println("----- Parameters Info. -----------------------------------")
        println("             Nc : 7")
        println("             Nf : $flavor")
        println("              Λ : 270 [MeV]")
        println("    Instanton # : $nI")
        println("AntiInstanton # : $nA")
        println("          Gauge : $gauge")
        println("----------------------------------------------------------")
    else 
        error("gauge is not valid!")
    end 
end


function Hamiltonian(ParamsofHMC,modelparams,colorworkspace)
    
    println("----- Hamiltonian Calc. --------------------------------------------")
    #運動量パート
    sumMomentum = 0.0

        # conjρ part
        sumMomentum += dot(ParamsofHMC.conjρ , ParamsofHMC.conjρ)
    
        #テスト(1) -----------------------------------------------------------

        if modelparams.gauge == "SU(2)" && modelparams.Nc == 3 
            println("test(1) passed!")
        elseif modelparams.gauge == "SU(3)" && modelparams.Nc == 7
            println("test(1) passed!")
        else
            error("couldn't pass test(1)!")
        end
        
        #--------------------------------------------------------------------

        # conjθ part (and conjz part)
        for j=1:modelparams.N
            sumMomentum += dot(ParamsofHMC.conjθ[j,1:modelparams.Nc],ParamsofHMC.conjθ[j,1:modelparams.Nc])
            # conjz part はまとめることもできるが，デバッグのために分離しておく．
            #sumMomentum += dot(ParamsofHMC.conjz[j,1:4],ParamsofHMC.conjz[j,1:4])
        end

        # conjz part

        for j=1:modelparams.N
            sumMomentum += dot(ParamsofHMC.conjz[j,1:4],ParamsofHMC.conjz[j,1:4])
        end

        sumMomentum *= 0.5
        @show sumMomentum
        #テスト(2) -----------------------------------------------------------
        println("test (2) passed!")
        #--------------------------------------------------------------------
    
    #運動量パート終

    #有効作用パート
    sumAction = 0.0
    #=
    有効作用の計算をする関数の引数を，全てここで作ったstructにする必要がある．
    そうでないと，引数を取り出してパラメータベクトルを作り，関数を呼び出すということになる．
    =#   
        # density part
        sumAction += -DensityPart(ParamsofHMC,modelparams)
        @show sumAction
        println("Density part passed!")
        # interacton part
        sumAction += InterAction(ParamsofHMC,modelparams,colorworkspace)
        println("Interaction part passed!")
        @show sumAction
        # determinant part



    #有効作用パート終わり


end

function LeapFog()

end

function MPtest()
    
end

function CalForce()

end
