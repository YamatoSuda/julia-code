#iDensity.jl

#=
# 2021-9-27 更新する前の旧バージョン（以下コメントアウト）

    function CNc(p)
        Nc, = p
        return 4.6*exp(-1.86*Nc)/(pi^2*factorial(Nc-1)*factorial(Nc-2))
    end

    function β1(x,p)
        Nc, Nf, Λ, = p 
        b = 11*Nc/3 - 2*Nf/3
        return -b*log(x*Λ)
    end

    function β2(x,p)
        Nc, Nf, Λ, = p 
        β = β1(x,p)
        b  = 11*Nc/3 - 2*Nf/3
        bb = 34*Nc^2/3 - 13*Nc*Nf/3 + Nf/Nc
        return β + bb*log(2*β/b)/(2*b)
    end

    function d_para(x,p)
        Nc, Nf, Λ, = p 
        b  = 11*Nc/3 - 2*Nf/3
        bb = 34*Nc^2/3 - 13*Nc*Nf/3 + Nf/Nc
        return CNc(p)*x^-5*β1(x,p)^(2*Nc)*exp(-β2(x,p)+(2*Nc-0.5*bb/b)*0.5*bb/b*log(β1(x,p))/β1(x,p))
    end

    #有効作用のインスタントン密度部分．作用内で用いるときは負号をつける．
    function density_part(ρ,p)
        return sum( log(d_para(ρ[i],p)) for i in eachindex(ρ) )
    end

    function SingleDensity(x,p)
        Nc, Nf, Λ, = p 
        b  = 11*Nc/3 - 2*Nf/3
        bb = 34*Nc^2/3 - 13*Nc*Nf/3 + Nf/Nc
        return CNc(p)*x^-5*β1(x,p)^(2*Nc)*exp(-β2(x,p)+(2*Nc-0.5*bb/b)*0.5*bb/b*log(β1(x,p))/β1(x,p))
    end
=#
    #新しい定義 ：2021-9-27 更新後の新バージョン（以下）

    function CNc(modelparams)
        Nc  = modelparams.Nc
        return 4.6exp(-1.86Nc)/(pi^2*factorial(Nc-1)*factorial(Nc-2))
    end

    function β1(x,modelparams)
        Nc = modelparams.Nc
        Nf = modelparams.Nf
        Λ  = modelparams.Λ
        b  = 11*Nc/3 - 2*Nf/3
        return -b * log( x * Λ )
    end

    function β2(x,modelparams)
        Nc = modelparams.Nc
        Nf = modelparams.Nf
        Λ  = modelparams.Λ
        β  = β1(x,modelparams)
        b  = 11*Nc/3 - 2*Nf/3
        bb = 34*Nc^2/3 - 13*Nc*Nf/3 + Nf/Nc
        return β + bb * log( 2 * β/b ) / ( 2 * b )
    end

    function SingleDensity(x,modelparams)
        Nc = modelparams.Nc
        Nf = modelparams.Nf
        Λ  = modelparams.Λ
        b  = 11*Nc/3 - 2*Nf/3
        bb = 34*Nc^2/3 - 13*Nc*Nf/3 + Nf/Nc
        return CNc(modelparams) * x^(-5) * β1(x,modelparams)^(2 * Nc) * exp( -β2(x,modelparams) + (2 * Nc - 0.5 * bb/b ) * 0.5 * (bb/b) * log( β1(x,modelparams)) / β1(x,modelparams))
    end

    function DensityPart(ParamsofHMC,modelparams)
        ρ = ParamsofHMC.ρ
        return sum( log(SingleDensity(ρ[i],modelparams)) for i in eachindex(ρ) )
    end


