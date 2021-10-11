#instanton_density.jl

    #新しい定義 ：2021-9-27 更新後の新バージョン（以下）

    function c_nc(Modelparams)
        Nc  = Modelparams.Nc
        return 4.66exp(-1.679Nc)/(pi^2*factorial(Nc-1)*factorial(Nc-2))
        #=
        cf. Instanton in QCD p.344, (93)
        =#
    end

    function β1(x,Modelparams)
        Nc = Modelparams.Nc
        Nf = Modelparams.Nf
        Λ  = Modelparams.Λ
        b  = 11*Nc/3 - 2*Nf/3
        return -b * log( x * Λ )
    end

    function β2(x,Modelparams)
        Nc = Modelparams.Nc
        Nf = Modelparams.Nf
        Λ  = Modelparams.Λ
        β  = β1(x,Modelparams)
        b  = 11*Nc/3 - 2*Nf/3
        bb = 34*Nc^2/3 - 13*Nc*Nf/3 + Nf/Nc
        return β + bb * log( 2 * β/b ) / ( 2 * b )
    end

    function single_density(x,Modelparams)
        Nc = Modelparams.Nc
        Nf = Modelparams.Nf
        Λ  = Modelparams.Λ
        b  = 11*Nc/3 - 2*Nf/3
        bb = 34*Nc^2/3 - 13*Nc*Nf/3 + Nf/Nc
        beta1 = β1(x,Modelparams)
        return c_nc(Modelparams) * x^(-5) * beta1^(2 * Nc) * exp( -β2(x,Modelparams) + (2 * Nc - 0.5 * bb/b ) * 0.5 * (bb/b) * log( beta1 ) / beta1)
    end

    function density_part(ParamsofHMC,Modelparams)
        ρ = ParamsofHMC.ρ
        return sum( log(single_density(ρ[i],Modelparams)) for i in eachindex(ρ) )
    end


