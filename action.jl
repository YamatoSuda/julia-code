#action.jl
    # conformal parameter : λ
#--------------------------------------------------------------------------#
    #新しく定義した関数     :2021-9-27 
    #関数名も変更          :2021-10-4

    function lambda(I,J,ParamsofHMC)
        ρ       = ParamsofHMC.ρ
        z       = ParamsofHMC.z
        RR      = sum( (z[I,i]-z[J,i])^2 for i=1:4 )
        #R_ratio = ( RR + ρ[I]^2 + ρ[J]^2 ) / ( ρ[I]*ρ[J] )
        R_ratio = (RR/(ρ[I]*ρ[J])) + (ρ[I]/ρ[J]) + (ρ[J]/ρ[I]) 
        return 0.5*R_ratio + 0.5*sqrt( R_ratio^2 - 4 )
    end

    function color_vector(thetaI,thetaJ,Colorworkspace)
        u      = Colorworkspace.u
        Sigma  = Colorworkspace.sigma
        Lambda = Colorworkspace.lambda
        tem    = Colorworkspace.tem
        res    = Colorworkspace.res
        mul!(tem,su2_matrix(thetaI,Sigma)',su2_matrix(thetaJ,Sigma)) 
        for i=1:4
            mul!(res,tem,tau_p(i))
            u[i] = 1/(2im)*tr(res)
        end
        return u
    end

    #2つのインスタントン間の中心を結ぶ単位ベクトル
    function unit_relative_vector(I,J,ParamsofHMC)
        z = ParamsofHMC.z
        r = zeros(Float64,4)
        R = sqrt(sum( (z[I,i]-z[J,i])^2 for i=1:4 ))
        if R == 0 
            error("@unit_relative_vector R = 0 !!")
        end
        for i=1:4
            r[i] = (z[I,i] - z[J,i])/R
        end
        return r
    end

    #2つのインスタントン間を結ぶ相対ベクトル
    function relative_vector(I,J,ParamsofHMC)
        z = ParamsofHMC.z
        R = zeros(Float64,4)
        for i=1:4
            R[i] = z[I,i] - z[J,i]
        end
        return R
    end

    # I,J でそれぞれインスタントンと反インスタントンの番号を指定する．
    function action_core(I,J,ParamsofHMC,Modelparams,Colorworkspace)
        θ            = ParamsofHMC.θ
        Nc           = Modelparams.Nc
        A            = Modelparams.A
        g            = Modelparams.g
        thetaI       = θ[I, 1:Nc]
        thetaJ       = θ[J, 1:Nc]
        u            = Colorworkspace.u
        
        λ    = 0.0
        uu   = 0.0
        # λの構成
        λ    = lambda(I,J,ParamsofHMC)
        # u²の構成
        u    = color_vector(thetaI,thetaJ,Colorworkspace)
        uu   = dot(u,u)

        #return (8*pi^2/g^2) * A*uu/λ^4 # S_core
        #@show A*uu/λ^4
        return A*uu/λ^4 # S_core/S_0
    end

    function action_stream(I,J,ParamsofHMC,Modelparams,Colorworkspace)
        θ            = ParamsofHMC.θ
        Nc           = Modelparams.Nc
        A            = Modelparams.A
        g            = Modelparams.g
        thetaI       = θ[I, 1:Nc]
        thetaJ       = θ[J, 1:Nc]
        u            = Colorworkspace.u
        
        λ     = 0.0
        s_int = 0.0
        term1 = 0.0
        term2 = 0.0
        term3 = 0.0
        R     = similar(u)

        # λの構成
        λ     = lambda(I,J,ParamsofHMC)
        # u,R_hatの構成
        u     = color_vector(thetaI,thetaJ,Colorworkspace)
        R_hat = unit_relative_vector(I,J,ParamsofHMC)
        # S_streamの計算
        #=
        出力が異常な値を返したら，sum(u.*R)の値の近くを確認する．R=0となり，unit_relative_vectorが
        定義されていない状況が考えられる．
        =#
        term1 = -4( 1 - λ^4 + 4λ^2 * log(λ) ) * ( dot(u,u) - 4*abs(sum(u.*R_hat))^2 )
        term2 = 2( 1 - λ^2 + (1 + λ^2) * log(λ) )
        term3 = (dot(u,u) - 4*abs(sum(u.*R_hat))^2)^2 + dot(u,u)^2 + 2*dot(u',u)*dot(u,u') 
        s_int = (term1 + term2*term3)*4/(λ^2-1)^3
        #return (8*pi^2/g^2) * s_int
        return s_int/4 # S_stream/S_0    4倍だけ違う
    end

    function action_int(ParamsofHMC,Modelparams,Colorworkspace)
        ret               = 0.0
        nI                = Modelparams.nI
        nA                = Modelparams.nA
        pair_numbers      = (1:nI+nA)
        pair_combinations = collect(combinations(pair_numbers,2))

        for i in eachindex(pair_combinations)
            I = pair_combinations[i][1]
            J = pair_combinations[i][2]
            if I <= nI && J <= nI  #II pair
                Sdum1  = action_core(I,J,ParamsofHMC,Modelparams,Colorworkspace)
                ret   += action_core(I,J,ParamsofHMC,Modelparams,Colorworkspace)
                #print("II! ")
                #println("($I,$J) ")
                #println("$Sdum1")
            elseif I > nI && J > nI #AA pair
                Sdum2  = action_core(I,J,ParamsofHMC,Modelparams,Colorworkspace)
                ret   += action_core(I,J,ParamsofHMC,Modelparams,Colorworkspace)
                #print("AA! ")
                #println("($I,$J) ")
                #println("$Sdum2")
            else  #IA pair
                Sdum3  = action_stream(I,J,ParamsofHMC,Modelparams,Colorworkspace)
                ret   += action_stream(I,J,ParamsofHMC,Modelparams,Colorworkspace)
                #print("IA! ")
                #println("($I,$J) ")
                #println("$Sdum3")
            end                             
        end
        return ret
    end
    
