#actions.jl
    # conformal parameter : λ
    # I,Jで1つ目，2つ目のインスタントンを指定する．
    # R=0,ρ1=ρ2のとき，λ=1
    function lambda(I,J,dynamical_p) 
        ρ, z,   = dynamical_p
        RR      = sum( (z[I,i]-z[J,i])^2 for i=1:4 )
        R_ratio = ( RR + ρ[I]^2 + ρ[J]^2 ) / ( ρ[I]*ρ[J] )
        return 0.5*R_ratio + 0.5*sqrt( R_ratio^2 - 4 )
    end

    function color_vector(su2I,su2J,dynamical_q)
        u, Sigma, tem, res = dynamical_q
        mul!(tem,SU2_matrix(su2I,Sigma)',SU2_matrix(su2J,Sigma)) 
        for i=1:4
            mul!(res,tem,tau_p(i))
            u[i] = 1/(2im)*tr(res)
        end
        return u
    end

    function relative_vector(I,J,z)
        r = zeros(Float64,4)
        R = sqrt(sum( (z[I,i]-z[J,i])^2 for i=1:4 ))
        for i=1:4
            r[i] = (z[I,i] - z[J,i])/R
        end
        return r
    end

    function space_ori_vector(I,J,z)
        R = zeros(Float64,4)
        for i=1:4
            R[i] = z[I,i] - z[J,i]
        end
        return R
    end

    # I,J でそれぞれインスタントンと反インスタントンの番号を指定する．
    # dynamical_p,dynamical_q で必要なパラメータを渡す
    function S_core(I,J,dynamical_p,dynamical_q)
        ρ, z, θ, A, g = dynamical_p
        u,            = dynamical_q
        su2I = θ[I, 1:3] # SU(2)の場合
        su2J = θ[J, 1:3] # SU(2)の場合
        
        λ    = 0.0
        uu   = 0.0
            
        # λの構成
        λ    = lambda(I,J,dynamical_p)
        
        # u²の構成
        u    = color_vector(su2I,su2J,dynamical_q)
        uu   = dot(u,u)
        
        #return (8*pi^2/g^2) * A*uu/λ^4 # S_core
        return A*uu/λ^4 # S_core/S_0
    end

    function S_stream(I,J,dynamical_p,dynamical_q)
        ρ, z, θ, A, g = dynamical_p
        u,            = dynamical_q
        su2I  = θ[I, 1:3] # SU(2)の場合
        su2J  = θ[J, 1:3] # SU(2)の場合
        
        λ     = 0.0
        s_int = 0.0
        term1 = 0.0
        term2 = 0.0
        term3 = 0.0
        R     = similar(u)
        
        # λの構成
        λ     = lambda(I,J,dynamical_p)
        
        # u,Rの構成
        u     = color_vector(su2I,su2J,dynamical_q)
        R     = relative_vector(I,J,z)
        
        # S_streamの計算
        term1 = -4( 1 - λ^4 + 4λ^2 * log(λ) ) * ( dot(u,u) - 4*abs(sum(u.*R))^2 )
        term2 = 2( 1 - λ^2 + (1 + λ^2) * log(λ) )
        term3 = (dot(u,u) - 4*abs(sum(u.*R))^2)^2 + dot(u,u)^2 + 2*dot(u',u)*dot(u,u') 
        s_int = (term1 + term2*term3)*4/(λ^2-1)^3
        #return (8*pi^2/g^2) * s_int
        return s_int/4 # S_stream/S_0    4倍だけ違う
    end

    #Sintを引数として，インスタントンと反インスタントンの数も含めてSintの計算を実行する．
    function S_interaction!(Sint,p,dynamical_p,dynamical_q)
        _, _, _, nI, nA = p
        IAnumbers = (1:N)
        IAcombinations = collect(combinations(IAnumbers,2))
        for i in eachindex(IAcombinations)
            I = IAcombinations[i][1]
            J = IAcombinations[i][2]
            if I <= nI && J <= nI  #II pair
                #Sdum1   = S_core(I,J,dynamical_p,dynamical_q)
                Sint += S_core(I,J,dynamical_p,dynamical_q)
                #print("II! ")
                #println("($I,$J) ")
                #println("$Sdum1")
            elseif I > nI && J > nI #AA pair
                #Sdum2   = S_core(I,J,dynamical_p,dynamical_q)
                Sint += S_core(I,J,dynamical_p,dynamical_q)
                #print("AA! ")
                #println("($I,$J) ")
                #println("$Sdum2")
            else  #IA pair
                #Sdum3   = S_stream(I,J,dynamical_p,dynamical_q)
                Sint += S_stream(I,J,dynamical_p,dynamical_q)
                #print("IA! ")
                #println("($I,$J) ")
                #println("$Sdum3")
            end                             
        end
        return Sint
    end

#--------------------------------------------------------------------------#
    #新しく定義した関数：2021-9-27 （以下コード）

    function lambda(I,J,ParamsofHMC)
        ρ       = ParamsofHMC.ρ
        z       = ParamsofHMC.z
        RR      = sum( (z[I,i]-z[J,i])^2 for i=1:4 )
        R_ratio = ( RR + ρ[I]^2 + ρ[J]^2 ) / ( ρ[I]*ρ[J] )
        return 0.5*R_ratio + 0.5*sqrt( R_ratio^2 - 4 )
    end

    function ColorVector(ParamofSUNI,ParamofSUNJ,colorworkspace)
        u      = colorworkspace.u
        Sigma  = colorworkspace.sigma
        Lambda = colorworkspace.lambda
        tem    = colorworkspace.tem
        res    = colorworkspace.res
        mul!(tem,SU2matrix(ParamofSUNI,Sigma)',SU2matrix(ParamofSUNJ,Sigma)) 
        for i=1:4
            mul!(res,tem,tau_p(i))
            u[i] = 1/(2im)*tr(res)
        end
        return u
    end

    function RelativeVector(I,J,ParamsofHMC)
        z = ParamsofHMC.z
        r = zeros(Float64,4)
        R = sqrt(sum( (z[I,i]-z[J,i])^2 for i=1:4 ))
        for i=1:4
            r[i] = (z[I,i] - z[J,i])/R
        end
        return r
    end

    # I,J でそれぞれインスタントンと反インスタントンの番号を指定する．
    function ActionCore(I,J,ParamsofHMC,modelparams,colorworkspace)
        θ            = ParamsofHMC.θ
        Nc           = modelparams.Nc
        A            = modelparams.A
        g            = modelparams.g
        ParamofSUNI  = θ[I, 1:Nc]
        ParamofSUNJ  = θ[J, 1:Nc]
        u            = colorworkspace.u
        
        λ    = 0.0
        uu   = 0.0
        # λの構成
        λ    = lambda(I,J,ParamsofHMC)
        # u²の構成
        u    = ColorVector(ParamofSUNI,ParamofSUNJ,colorworkspace)
        uu   = dot(u,u)

        #return (8*pi^2/g^2) * A*uu/λ^4 # S_core
        #@show A*uu/λ^4
        return A*uu/λ^4 # S_core/S_0
    end

    function ActionStream(I,J,ParamsofHMC,modelparams,colorworkspace)
        θ            = ParamsofHMC.θ
        Nc           = modelparams.Nc
        A            = modelparams.A
        g            = modelparams.g
        ParamofSUNI  = θ[I, 1:Nc]
        ParamofSUNJ  = θ[J, 1:Nc]
        u            = colorworkspace.u
        
        λ     = 0.0
        s_int = 0.0
        term1 = 0.0
        term2 = 0.0
        term3 = 0.0
        R     = similar(u)

        # λの構成
        λ     = lambda(I,J,ParamsofHMC)
        # u,Rの構成
        u     = ColorVector(ParamofSUNI,ParamofSUNJ,colorworkspace)
        R     = RelativeVector(I,J,ParamsofHMC)
        # S_streamの計算
        @show term1 = -4( 1 - λ^4 + 4λ^2 * log(λ) ) * ( dot(u,u) - 4*abs(sum(u.*R))^2 )
        @show term2 = 2( 1 - λ^2 + (1 + λ^2) * log(λ) )
        @show term3 = (dot(u,u) - 4*abs(sum(u.*R))^2)^2 + dot(u,u)^2 + 2*dot(u',u)*dot(u,u') 
        @show s_int = (term1 + term2*term3)*4/(λ^2-1)^3
        #return (8*pi^2/g^2) * s_int
        return s_int/4 # S_stream/S_0    4倍だけ違う
    end

    function InterAction(ParamsofHMC,modelparams,colorworkspace)
        ret            = 0.0
        nI             = modelparams.nI
        nA             = modelparams.nA
        IAnumbers      = (1:nI+nA)
        IAcombinations = collect(combinations(IAnumbers,2))

        for i in eachindex(IAcombinations)
            I = IAcombinations[i][1]
            J = IAcombinations[i][2]
            if I <= nI && J <= nI  #II pair
                Sdum1  = ActionCore(I,J,ParamsofHMC,modelparams,colorworkspace)
                ret   += ActionCore(I,J,ParamsofHMC,modelparams,colorworkspace)
                print("II! ")
                println("($I,$J) ")
                println("$Sdum1")
            elseif I > nI && J > nI #AA pair
                Sdum2  = ActionCore(I,J,ParamsofHMC,modelparams,colorworkspace)
                ret   += ActionCore(I,J,ParamsofHMC,modelparams,colorworkspace)
                print("AA! ")
                println("($I,$J) ")
                println("$Sdum2")
            else  #IA pair
                Sdum3  = ActionStream(I,J,ParamsofHMC,modelparams,colorworkspace)
                ret   += ActionStream(I,J,ParamsofHMC,modelparams,colorworkspace)
                print("IA! ")
                println("($I,$J) ")
                println("$Sdum3")
            end                             
        end
        return ret
    end
    
