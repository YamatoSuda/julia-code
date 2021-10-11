#hmc.jl

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

function choose_gauge_field!(params,nI,nA;gauge="SU(2)",flavor=0)
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
#=
        println("----- Parameters Info. ---------------------------------------------------------------------")
        println("             Nc : 3")
        if flavor == 0
        println("         Flavor : pure gauge")
        else
        println("             Nf : $flavor")
        end
        println("              Λ : 270 [MeV]")
        println("    Instanton # : $nI")
        println("AntiInstanton # : $nA")
        println("    Color Gauge : $gauge")
        println("--------------------------------------------------------------------------------------------")
 =#       
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
        println("----- Parameters Info. -----------------------------------------------------------------------")
        println("             Nc : 7")
        if flavor == 0
        println("         Flavor : pure gauge")
        else
        println("             Nf : $flavor")
        end
        println("              Λ : 270 [MeV]")
        println("    Instanton # : $nI")
        println("AntiInstanton # : $nA")
        println("    Color Gauge : $gauge")
        println("----------------------------------------------------------------------------------------------")
    else 
        error("gauge is not valid!")
    end 
end


function hamiltonian(ParamsofHMC,Modelparams,Colorworkspace)
    
    println("----- hamiltonian Calc. --------------------------------------------------------------------")
    #運動量パート
    sum_momentum = 0.0

        # conjρ part
        sum_momentum += dot(ParamsofHMC.conjρ , ParamsofHMC.conjρ)
    
        #テスト(1) -----------------------------------------------------------

        if Modelparams.gauge == "SU(2)" && Modelparams.Nc == 3 
            println("test(1) : [correct gauge?] passed!")
        elseif Modelparams.gauge == "SU(3)" && Modelparams.Nc == 7
            println("test(1) : [correct gauge?] passed!")
        else
            error("couldn't pass test(1)!")
        end
        
        #--------------------------------------------------------------------

        # conjθ part (and conjz part)
        for j=1:Modelparams.N
            sum_momentum += dot(ParamsofHMC.conjθ[j,1:Modelparams.Nc],ParamsofHMC.conjθ[j,1:Modelparams.Nc])
            # conjz part はまとめることもできるが，デバッグのために分離しておく．
            #sum_momentum += dot(ParamsofHMC.conjz[j,1:4],ParamsofHMC.conjz[j,1:4])
        end

        # conjz part

        for j=1:Modelparams.N
            sum_momentum += dot(ParamsofHMC.conjz[j,1:4],ParamsofHMC.conjz[j,1:4])
        end

        sum_momentum *= 0.5
        @show sum_momentum
        #テスト(2) -----------------------------------------------------------
        println("test(2) : [finish calculating momenta?] passed!")
        #--------------------------------------------------------------------
    #運動量パート終

    #有効作用パート
    sum_action = 0.0
    det_ret    = 0.0
    det_ret_vec = zeros(3)
    #=
    有効作用の計算をする関数の引数を，全てここで作ったstructにする必要がある．
    そうでないと，引数を取り出してパラメータベクトルを作り，関数を呼び出すということになる．
    =#   
        # density part
        sum_action += -density_part(ParamsofHMC,Modelparams)
        @show sum_action
        println("density part passed!")
        # interacton part
        sum_action += action_int(ParamsofHMC,Modelparams,Colorworkspace)
        println("action_int part passed!")
        @show sum_action
        # determinant part
        det_ret, det_ret_vec = det_stream_no2(1,2,ParamsofHMC,Modelparams,Colorworkspace)
        # sum_action += det_ret
        # sum_action += - Modelparams.Nf * tr

    #有効作用パート終わり
        #テスト(3) -----------------------------------------------------------
        println("test(3) : [finish calculating action?] passed!")
        #--------------------------------------------------------------------
        println("----- hamiltonian Calc. End ----------------------------------------------------------------")
    #

    #行列式パート
    sum_det = 0.0
    #=
    2021-9-29
    fermion行列式の計算は検討が必要なので，あとで加えることができる形でHMC法の実装を行う．
    =#
        #ここに行列式の実装をする

    #= 行列式パート終わり =#
    hamiltonian = sum_action + sum_momentum + sum_det
    println("hamiltonian : $hamiltonian")
    return hamiltonian
end

#=
分子軌道法におけるパラメータはtauの時間幅，tauの時間の分割数
これらを入力するパラメータとする．デフォルト値の設定もしておく必要がありそうか．

ρ，θ，zで分けて実装する．
ρはインスタントンを区別する添字 N で指定することができる．
θとzはインスタントンをしているする添字 N と何番目のパラメータかを指定する添字 k で指定できる．

これらを区別して実装する．区別するのは，dSdx = cal_force の実装を分けることを意味する．

=#

#N(=1...nI+nA)番目のrhoを１回更新する
function molecular_dynamics_rho!(N, ParamsofHMC, Modelparams, Colorworkspace, Nt=10, dt=0.01)
    println("----- Molecular Dynamics Calc. of rho -------------------------------------------------------------")
        println("Nt    : $Nt")
        println("dt    : $dt")
        println("Nt*dt : $(Nt*dt)")
        Λ         = Modelparams.Λ
        dum_force = 0.0
    
        #変数を個別に更新する場合
            #ParamsofHMC.ρの更新
                # j=0->1 ステップ
                P = ParamsofHMC.conjρ[N]
                X = ParamsofHMC.ρ[N] + P*0.5dt
                println("P = $P, X = $X")
                # j=1->N-1 ステップ
                for j=1:Nt-1
                    println("before X=$X")
                    dum_force = cal_force_rho(X, N, ParamsofHMC, Modelparams, Colorworkspace)*dt
                    @show dum_force
                    P = P - dum_force
                    X = X + P*dt
                    println("j=$j  P=$P, X=$X")
                    #=
                    @show j
                    @show P
                    @show X
                    =#
                end
                 # j=N-1->N ステップ
                P = P - cal_force_rho(X, N, ParamsofHMC, Modelparams, Colorworkspace)*dt
                X = X + P*0.5dt
                @show X
                
            ParamsofHMC.ρ[N] = X
        #変数を塊で更新する場合


        #変数を全て更新する場合


    println("----- Molecular Dynamics Calc. of rho End ---------------------------------------------------------")
end


function cal_force_rho(x, N, ParamsofHMC, Modelparams, Colorworkspace)

    densitypart = 0.0
    intpart     = 0.0
    detpart     = 0.0
    ret         = 0.0

    nI          = Modelparams.nI
    nA          = Modelparams.nA
    Nc          = Modelparams.Nc
    Nf          = Modelparams.Nf
    #=
    @show Nc
    @show Nf
    @show x
    @show Modelparams.Λ
    =#
    #Single Instanton Densityからの寄与------------------------------------------------------------------------------
        b             = 11*Nc/3 - 2*Nf/3
        bb            = 34*Nc^2/3 - 13*Nc*Nf/3 + Nf/Nc
        RatioB        = 0.5bb/b
        ColFlaFactor  = 2Nc - RatioB
        beta1         = β1(x,Modelparams)
        Single        = single_density(x,Modelparams)
        phase         = -β2(x,Modelparams) + ColFlaFactor * RatioB * log(beta1) / beta1
        factor        = 1 + RatioB * (beta1 + ColFlaFactor * (log(beta1)-1) ) / beta1^2
        DensityFactor = -5beta1 - 2Nc*b + beta1 * b * factor
        densitypart   = -c_nc(Modelparams) * x^(-6) * (beta1)^(2Nc-1) * exp(phase) * DensityFactor / Single
    #--------------------------------------------------------------------------------------------------------------


    #action_int からの寄与------------------------------------------------------------------------------------------

        #=  
        入力されたN,nI,nAの情報から，IIペアがどの数字の組か，IAペアがどの数字の組か，を判断すれば，
        ノートに書いてあるlambdaを用いた表式から計算することができる．

        例) N=2, nI=2, nA=2 の場合
        instanton,antiinstantonを合わせて，i=1,2,3,4と番号を振ることにする．
        このとき，考えるペアは (2,1),(2,3),(2,4)であり，それぞれ順に II,IA,IAペアとなる．
        この判断を N,nI,nAの情報からすることが目標．
        この部分が実装できれば，相互作用部分は計算することができる．
        =#

        #ペアを作る
        pairs = makepair(N, nI, nA)

        for ii = 1:(nI+nA-1)
            N1 = pairs[ii][1]
            N2 = pairs[ii][2]
            if (N1 <= nI) && (N2 <= nA)
                #@show pairs[ii]
                #println("II")
                #= core =#
                intpart += action_core_drho(N, N1, N2, ParamsofHMC, Modelparams, Colorworkspace)
            elseif (N1 <= nI) && (N2 > nI)
                #@show pairs[ii]
                #println("IA")
                #= core + stream =#
                intpart += ( 
                    action_core_drho(N, N1, N2, ParamsofHMC, Modelparams, Colorworkspace)
                  + action_stream_drho(N, N1, N2, ParamsofHMC, Modelparams, Colorworkspace)
                   )
            elseif (N1 > nI) && (N2 <= nI)
                #@show pairs[ii]
                #println("AI")
                #= core + stream =#
                intpart += (
                    action_core_drho(N, N1, N2, ParamsofHMC, Modelparams, Colorworkspace)
                  + action_stream_drho(N, N1, N2, ParamsofHMC, Modelparams, Colorworkspace)
                   )
            else 
                #@show pairs[ii]
                #println("AA")
                #= core =#
                intpart += action_core_drho(N, N1, N2, ParamsofHMC, Modelparams, Colorworkspace)
            end
        end
    
    #--------------------------------------------------------------------------------------------------------------

    #Fermion determinant からの寄与----------------------------------------------------------------------------------
    #= まだ実装できていない =#


    #--------------------------------------------------------------------------------------------------------------
    
    println("======================================")
    @show densitypart 
    @show intpart
    @show detpart
    println("======================================")
    
    #--------------------------------------------------------------------------------------------------------------
    if Nf == 0
        ret = densitypart + intpart
    else
        ret = densitypart + intpart + detpart
    end
    #--------------------------------------------------------------------------------------------------------------
    @show ret
    return real(ret)
end

function makepair(N, nI, nA; num=nI+nA)
    if N > num ; error("func. MakePair error."); end
    com = map(1:num) do i; (N, i) end
    num >= N && deleteat!(com,N)
    com
end

function action_core_drho(i, l, m, ParamsofHMC, Modelparams, Colorworkspace)
    
    θ            = ParamsofHMC.θ
    ρ            = ParamsofHMC.ρ
    z            = ParamsofHMC.z
    Nc           = Modelparams.Nc

    g            = Modelparams.g
    A            = Modelparams.A
    thetaI       = θ[l, 1:Nc]
    thetaJ       = θ[m, 1:Nc]
    u            = Colorworkspace.u 
    λ            = 0.0
    uu           = 0.0
    RR           = 0.0

    λ            = lambda(l,m,ParamsofHMC)
    u            = color_vector(thetaI,thetaJ,Colorworkspace)
    uu           = dot(u,u)

    RR           = sum( (z[l,ii]-z[m,ii])^2 for ii=1:4 )
    square_rho   = RR + ρ[l]*ρ[l] + ρ[m]*ρ[m]

    prefactor    = (8pi^2/g^2)*A*uu*(-4/λ^5)
    dλdρ         = 0.5 * ( 1 + ( ( square_rho/(ρ[l] * ρ[m]) )^2 - 4 )^(-0.5) * square_rho/( ρ[l] * ρ[m] ) ) * ( (2*ρ[l] - square_rho/ρ[l])*δ(i,l) + (2*ρ[m]-square_rho/ρ[m])*δ(i,m) )/(ρ[l]*ρ[m]) 

    return prefactor * dλdρ
end

function action_stream_drho(i, l, m, ParamsofHMC, Modelparams, Colorworkspace)
    
    θ            = ParamsofHMC.θ
    ρ            = ParamsofHMC.ρ
    z            = ParamsofHMC.z
    Nc           = Modelparams.Nc

    g            = Modelparams.g
    A            = Modelparams.A
    thetaI       = θ[l, 1:Nc]
    thetaJ       = θ[m, 1:Nc]
    u            = Colorworkspace.u 
    R            = similar(u)
    R_hat        = similar(u)
    λ            = 0.0
    uu           = 0.0
    RR           = 0.0


    RR           = sum( (z[l,ii]-z[m,ii])^2 for ii=1:4 )
    square_rho   = RR + ρ[l]*ρ[l] + ρ[m]*ρ[m]

    dλdρ         =  0.5 * ( 1 + ( ( square_rho/(ρ[l] * ρ[m]) )^2 - 4 )^(-0.5) * square_rho/( ρ[l] * ρ[m] ) ) *( (2*ρ[l] - square_rho/ρ[l])*δ(i,l) + (2*ρ[m]-square_rho/ρ[m])*δ(i,m) )/(ρ[l]*ρ[m]) 

    λ            = lambda(l,m,ParamsofHMC)
    u            = color_vector(thetaI,thetaJ,Colorworkspace)
    uu           = dot(u,u)
    R            = relative_vector(l,m,ParamsofHMC)
    R_hat        = unit_relative_vector(l,m,ParamsofHMC)


    term1 = - (24λ/(λ^2  - 1)^4) * ( 
        - 4( 1 - λ^4 + 4λ^2 * log(λ))*( uu-4*abs(sum(u.*R_hat))^2)
        + 2( 1 - λ^2 + (1 + λ^2) * log(λ) ) 
        * ( (uu - 4*abs(sum(u.*R_hat))^2)^2 + uu^2 + 2*dot(u',u)*dot(u,u') )
    )

    term2 = (4/(λ^2 - 1)^3) * (
          16( λ^3 - λ - 2λ*log(λ) )*( uu - 4*abs(sum(u.*R_hat))^2)
        + 2( -λ + 1/λ + 2λ*log(λ) )*( (uu - 4*abs(sum(u.*R_hat))^2)^2 + uu^2 + 2*dot(u',u)*dot(u,u') )
    )

    return dλdρ * (term1 + term2 )*(8pi^2/g^2)
end

#行列式のrho微分
function action_det_drho()

end



#N(=1...nI+nA)番目のインスタントンのk(=1...Nc)番目のθを１回更新する．
function molecular_dynamics_theta!(N, k, ParamsofHMC, Modelparams, Colorworkspace, Nt=10, dt=0.01)
    println("----- Molecular Dynamics Calc. of theta -------------------------------------------------------------")
        println("Nt    : $Nt")
        println("dt    : $dt")
        println("Nt*dt : $(Nt*dt)")
    
        #変数を個別に更新する場合
            #ParamsofHMC.θの更新
            # j=0->1 ステップ
            P = ParamsofHMC.conjθ[N, k]
            X = ParamsofHMC.θ[N, k] + P*0.5dt
            println("P = $P, X = $X")
            # j=1->N-1 ステップ
            for j=1:Nt-1
                println("before X=$X")
                dum_force = cal_force_theta(X, N, k, ParamsofHMC, Modelparams, Colorworkspace)*dt
                @show dum_force
                P = P - dum_force
                X = X + P*dt
                println("j=$j  P=$P, X=$X")
                #=
                @show j
                @show P
                @show X
                =#
            end
             # j=N-1->N ステップ
            P = P - cal_force_theta(X, N, k, ParamsofHMC, Modelparams, Colorworkspace)*dt
            X = X + P*0.5dt
            @show X
            
        ParamsofHMC.θ[N, k] = X
        #変数を塊で更新する場合


        #変数を全て更新する場合


    println("----- Molecular Dynamics Calc. of theta End ---------------------------------------------------------")
end

#i番絵mのインスタントンのk番目のθを時間発展させる力を計算する．

function dudθ(mu, l, m, a, i, ParamsofHMC, Modelparams, Colorworkspace)
    ret_dudtheta = 0.0
    sigma        = make_sigma()
    t_matrix     = sigma * 0.5    #SU(2)の生成子 t_a(a=1,2,3)
    theta_al     = ParamsofHMC.θ[l, a]   #l番目インスタントンのa番目θ
    #@show theta_al
    theta_am     = ParamsofHMC.θ[m, a]   #m番目インスタントンのa番目θ
    theta_l      = ParamsofHMC.θ[l]      #l番目インスタントンのθ
    theta_m      = ParamsofHMC.θ[m]      #m番目インスタントンのθ
    
    U_l          = su2_matrix(theta_l,sigma)
    U_m          = su2_matrix(theta_m,sigma)
    first_term   = - t_matrix[a] * exp( - im * theta_al * t_matrix[a] ) * U_m * tau_p(mu) * δ(i,l)
    second_term  =   U_l' * t_matrix[a] * exp(  im * theta_am * t_matrix[a] ) * tau_p(mu) * δ(i,m)

    ret_dudtheta = 0.5 * tr( first_term + second_term )

    return ret_dudtheta
end

#l,m:インスタントン,反インスタントンの添字,a:θのカラーの添字,i:インスタントンの添字.
function cal_force_theta(X, l, m, a, i, ParamsofHMC, Modelparams, Colorworkspace)
    A = Modelparams.A
    g = Modelparams.g
    pre_factor = 8pi^2/g^2
    λ            = 0.0

    λ            = lambda(l,m,ParamsofHMC)
end

function action_core_dtheta(l, m, a, i, ParamsofHMC, Modelparams, Colorworkspace)
    A            = Modelparams.A
    g            = Modelparams.g
    pre_factor   = 8pi^2/g^2
    main_factor  = 0.0
    λ            = 0.0

    # l,m に対する color vector : u
    theta_l      = ParamsofHMC.θ[l]      #l番目インスタントンのθ
    theta_m      = ParamsofHMC.θ[m]      #m番目インスタントンのθ
    u            = color_vector(theta_l,theta_m,Colorworkspace)

    # l,m に対する lambda
    λ            = lambda(l,m,ParamsofHMC)

    # main factor (= d(|u|^2)/dθ )
    main_factor  = sum( dudθ(mu, l, m, a, i, ParamsofHMC, Modelparams, Colorworkspace)' * u[mu] + u[mu]' * dudθ(mu, l, m, a, i, ParamsofHMC, Modelparams, Colorworkspace) for mu=1:4 )
    
    return pre_factor * (A/λ^4) * main_factor 
end

function action_stream_dtheta(l, m, a, i, ParamsofHMC, Modelparams, Colorworkspace)
    g            = Modelparams.g
    u            = Colorworkspace.u
    R_hat        = similar(u)

    pre_factor   = 8pi^2/g^2
    λ            = 0.0
    f1           = 0.0
    f2           = 0.0
    f3           = 0.0
    g1           = 0.0
    dg1dθ        = 0.0
    dg2dθ        = 0.0
    duudθ        = 0.0
    theta_l      = ParamsofHMC.θ[l]      #l番目インスタントンのθ
    theta_m      = ParamsofHMC.θ[m]      #m番目インスタントンのθ

    λ            = lambda(l,m,ParamsofHMC)
    R_hat        = unit_relative_vector(l,m,ParamsofHMC)

    u            = color_vector(theta_l,theta_m,Colorworkspace)
    duudθ        = sum( dudθ(mu, l, m, a, i, ParamsofHMC, Modelparams, Colorworkspace)' * u[mu] + u[mu]' * dudθ(mu, l, m, a, i, ParamsofHMC, Modelparams, Colorworkspace) for mu=1:4 )

    w            = sum( u[mu] * R_hat[mu] for mu=1:4 )
    dwdθ         = sum( dudθ(mu, l, m, a, i, ParamsofHMC, Modelparams, Colorworkspace)  * R_hat[mu] for mu=1:4 )

    alpha        = dot(u',u)
    dalphadθ     = 2 * sum( u[mu] * dudθ(mu, l, m, a, i, ParamsofHMC, Modelparams, Colorworkspace) for mu=1:4 )
    beta         = dot(u,u')
    dbetadθ      = 2 * sum( u[mu]' * dudθ(mu, l, m, a, i, ParamsofHMC, Modelparams, Colorworkspace)' for mu=1:4 )

    #各ファクターの計算
    f1    = func_f1(λ)
    f2    = func_f2(λ)
    f3    = func_f3(λ)

    g1    = dot(u,u) - 4*abs(sum(u.*R_hat))^2
    dg1dθ = duudθ - 4*( dwdθ' * w + w' * dwdθ )
    dg2dθ = 2*dot(u,u) * duudθ + 2 * ( dalphadθ * beta + alpha * dbetadθ )

    return pre_factor * f1 * ( dg1dθ * f2 + f3 * ( 2 * g1 * dg1dθ + dg2dθ ) )
end

func_f1(x)    = 4 / ( x^2 - 1 )^3
func_f2(x)    = -4 * ( 1 - x^4 + 4*x^2*log(x))
func_f3(x)    = 2 * ( 1 - x^2 + ( 1 + x^2 ) * log(x) )
func_df1dx(x) = -24 * x / ( x^2 - 1 )^4
func_df2dx(x) = 16 * ( x^3 - 2*x*log(x) - x )
func_df3dx(x) = 2 * ( -x + 1/x + 2*x*log(x) )
#func_g1(u) = dot(u,u) - 4*abs(sum(u.*R_hat))^2
#func_g2(u) = (dot(u,u))^2 + 2*dot(u',u)*dot(u,u')

function action_det_dtheta()
    
end

#N(=1...nI+nA)番目のインスタントンのk(=1,2,3,4)番目のzを１回更新する．
function molecular_dynamics_z!(N, k, ParamsofHMC, Modelparams, Colorworkspace, Nt=10, dt=0.01)
    println("----- Molecular Dynamics Calc. of z -------------------------------------------------------------")
        println("Nt    : $Nt")
        println("dt    : $dt")
        println("Nt*dt : $(Nt*dt)")
    
        #変数を個別に更新する場合
            #ParamsofHMC.zの更新
            # j=0->1 ステップ
            P = ParamsofHMC.conjz[N, k]
            X = ParamsofHMC.z[N, k] + P*0.5dt
            println("P = $P, X = $X")
            # j=1->N-1 ステップ
            for j=1:Nt-1
                println("before X=$X")
                dum_force = cal_force_z(X, N, k, ParamsofHMC, Modelparams, Colorworkspace)*dt
                @show dum_force
                P = P - dum_force
                X = X + P*dt
                println("j=$j  P=$P, X=$X")
                #=
                @show j
                @show P
                @show X
                =#
            end
             # j=N-1->N ステップ
            P = P - cal_force_z(X, N, k, ParamsofHMC, Modelparams, Colorworkspace)*dt
            X = X + P*0.5dt
            @show X
            
        ParamsofHMC.z[N] = X

        #変数を塊で更新する場合


        #変数を全て更新する場合


    println("----- Molecular Dynamics Calc. of z End ---------------------------------------------------------")
end

function cal_force_z(N, k, ParamsofHMC, Modelparams, Colorworkspace)

end

function action_core_dz(l, m, mu, i, ParamsofHMC, Modelparams, Colorworkspace)
    ρ            = ParamsofHMC.ρ
    g            = Modelparams.g
    A            = Modelparams.A
    u            = Colorworkspace.u
    R            = similar(u)
    λ            = 0.0
    main_factor  = 0.0
    ratio_factor = 0.0

    rho_l        = ρ[l]
    rho_m        = ρ[m]

    pre_factor1  = -32pi^2/g^2
    λ            = lambda(l,m,ParamsofHMC)
    R            = relative_vector(l,m,ParamsofHMC)
    pre_factor2  = ( A * dot(u,u) / λ^5 ) * ( δ(i,l) - δ(i,m) ) * R[mu] / ( rho_l * rho_m )
    ratio_factor = ( norm(R)^2 + rho_l^2 + rho_m^2 ) / ( rho_l * rho_m )
    main_factor  = 1 + ratio_factor / sqrt(ratio_factor^2 - 4) 
    
    @show λ
    @show pre_factor1
    @show pre_factor2
    @show main_factor

    return pre_factor1 * pre_factor2 * main_factor
end



function action_stream_dz(l, m, mu, i, ParamsofHMC, Modelparams, Colorworkspace)
    ρ            = ParamsofHMC.ρ
    z            = ParamsofHMC.z
    g            = Modelparams.g
    u            = Colorworkspace.u
    R            = similar(u)
    R_hat        = similar(u)

    pre_factor   = 8pi^2/g^2
    λ            = 0.0
    λ            = lambda(l,m,ParamsofHMC)
    rho_l        = ρ[l]
    rho_m        = ρ[m]
    z_mu_l       = z[l,mu] #添字の順序に注意．これが正しい順番．
    z_mu_m       = z[m,mu]

    R            = relative_vector(l,m,ParamsofHMC)
    R_hat        = unit_relative_vector(l,m,ParamsofHMC)
    norm_R       = norm(R)

    f1           = func_f1(λ)
    f2           = func_f2(λ)
    f3           = func_f3(λ)
    g1           = dot(u,u) - 4*abs(sum(u.*R_hat))^2
    g2           = (dot(u,u))^2 + 2*dot(u',u)*dot(u,u')


    ratio_factor = ( norm_R^2 + rho_l^2 + rho_m^2 ) / ( rho_l * rho_m )
    dλdz         = ( δ(i,l) - δ(i,m) ) * R[mu] / ( rho_l * rho_m ) * ( 1 + ratio_factor / sqrt(ratio_factor^2 - 4) )
    df1dλ        = func_df1dx(λ)
    df2dλ        = func_df2dx(λ)
    df3dλ        = func_df3dx(λ)
    df1dz        = dλdz * df1dλ
    df2dz        = dλdz * df2dλ
    df3dz        = dλdz * df3dλ

    w            = sum( u[mu] * R_hat[mu] for mu=1:4 )
    dwdz         = u[mu] * ( δ(i,l) - δ(i,m) ) * (1/norm_R) * ( 1 - ( ( z_mu_l - z_mu_m ) / norm_R )^2 )
    dg1dz        = -4 * (dwdz' * w + w' * dwdz )

    sub_factor1  = df1dz*f2*g1 + f1*df2dz*g1 + f1*f2*dg1dz
    sub_factor2  = ( df1dz*f3 + f1*df3dz )*( g1^2 + g2 ) + 2*f1*f2*g1*dg1dz

    @show λ
    @show sub_factor1
    @show sub_factor2
    @show pre_factor

    return pre_factor * ( sub_factor1 + sub_factor2 )
end


function action_det_dz()
    
end


function mp_check()
    
end
