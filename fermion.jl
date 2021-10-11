#fermion.jl

function det_integral_old(λ)
    det_integrand(x) = 6.0 * x^(3/2) / ((x + 1 / λ)^(3/2)*(x + λ)^(5/2))
    ret, err = quadgk(det_integrand,0.0,Inf)
end

function det_integral(f,λ)
    det_integrand(x) = f(x,λ)
    ret, err = quadgk(det_integrand,0.0,Inf)
end

#=
行列式を計算するときは，No.2を用いる．No.1は数値積分の結果を用いているが，No.2は解析的な積分を用いている．
比較した結果，論文の図を再現するのはNo.2の方だということがわかっている．

2021-10-9 追記
次元解析の結果，T要素の標識は (u・R) の形で含まれている必要がある．
これは，TがL^-1の次元を持つことに由来する．
これに伴い，コードの書き換えを行う．
=#

#structを用いた改訂版関数
function det_stream_no1(I,J,ParamsofHMC,Modelparams,Colorworkspace)
    ρ               = ParamsofHMC.ρ
    θ               = ParamsofHMC.θ
    Nc              = Modelparams.Nc
    u               = Colorworkspace.u
    Sigma           = Colorworkspace.sigma
    Lambda          = Colorworkspace.lambda
    R_hat           = similar(u)
    tau_times_R_hat = zeros(ComplexF64,2,2)
    prod_su2matrix  = similar(tau_times_R_hat)
    prod_matrix     = similar(tau_times_R_hat)
    prod_matrix2    = similar(tau_times_R_hat)

    ret_vector      = zeros(3)              #返り値を格納するベクトル

    thetaI          = θ[I, 1:Nc]
    thetaJ          = θ[J, 1:Nc]

    λ               = 0.0
    F               = 0.0
    error           = 0.0
    tr_matrix       = 0.0
    tr_matrix2      = 0.0
    uprodR          = Complex(0,0)

    # λの構成，F(λ)の計算
    λ               = lambda(I,J,ParamsofHMC)
    integ_func(x,a) = 6.0 * x^(3/2) / ((x + 1 / a)^(3/2)*(x + a)^(5/2))
    F, error        = det_integral(integ_func,λ)                              #F(λ) = 6.0 \int f(x,a) dx を返り値としている
    ret_vector[1]   = F
    ret_vector[2]   = error

    #行列トレース部分の計算

    R_hat           = unit_relative_vector(I,J,ParamsofHMC)
    tau_times_R_hat = tau_p(1)*R_hat[1] + tau_p(2)*R_hat[2] + tau_p(3)*R_hat[3] + tau_p(4)*R_hat[4]
    mul!(prod_su2matrix, su2_matrix(thetaI,Sigma)', su2_matrix(thetaJ,Sigma))
    mul!(prod_matrix, tau_times_R_hat, prod_su2matrix)
    tr_matrix       = tr(prod_matrix)
    geometric_mean  = sqrt(ρ[I]*ρ[J])
    #tr_matrix2  = tr(prod_su2matrix*tau_p(1)) * R_hat[1] + tr(prod_su2matrix*tau_p(2)) * R_hat[2] + tr(prod_su2matrix*tau_p(3)) * R_hat[3] + tr(prod_su2matrix*tau_p(4)) * R_hat[4] 

    #↓rhoに依存しない部分に確認
    ret_vector[3]    = abs(0.5 * tr_matrix)

    #それ以外の部分の確認
    ret_vector[2]    = F/geometric_mean

    val = - 0.5 * tr_matrix * F / geometric_mean
    #return abs(val), ret_vector #テスト用の返り値
    return ret_vector[2], ret_vector
    #return val,      ret_vector
end

function det_stream_no2(I,J,ParamsofHMC,Modelparams,Colorworkspace)
    ρ               = ParamsofHMC.ρ
    θ               = ParamsofHMC.θ
    Nc              = Modelparams.Nc
    u               = Colorworkspace.u
    Sigma           = Colorworkspace.sigma
    Lambda          = Colorworkspace.lambda
    tem             = Colorworkspace.tem
    res             = Colorworkspace.res
    R               = similar(u)
    R_hat           = similar(u)
    prod_su2matrix  = zeros(ComplexF64,2,2)
    prod_matrix     = similar(prod_su2matrix)
    prod_matrix2    = similar(prod_su2matrix)

    ret_vector       = zeros(3)              #返り値を格納するベクトル

    thetaI          = θ[I, 1:Nc]
    thetaJ          = θ[J, 1:Nc]

    λ               = 0.0
    F               = 0.0
    u_times_R       = 0.0
    u_times_R_hat   = 0.0

    #λの構成
    λ               = lambda(I,J,ParamsofHMC)

    #行列のトレース部分の計算
    R               = relative_vector(I,J,ParamsofHMC)
    R_hat           = unit_relative_vector(I,J,ParamsofHMC)
    u               = color_vector(thetaI,thetaJ,Colorworkspace)
    u_times_R       = u[1]*R[1]     + u[2]*R[2]     + u[3]*R[3]     + u[4]*R[4]
    u_times_R_hat   = u[1]*R_hat[1] + u[2]*R_hat[2] + u[3]*R_hat[3] + u[4]*R_hat[4]
    geometric_mean  = sqrt(ρ[I]*ρ[J])
    #@show u_times_R_hat
    #@show 2im * u_times_R_hat

    c1              = 3pi/8
    c2              = (3pi/32)^(4/3)
    F               = c1 * λ^(3/2) / (1+1.25*(λ^2-1) + c2*(λ^2-1)^2)^(3/4)
    ret_vector[1]   = F
    ret_vector[2]   = 0.0
    #@show F

    #↓rhoに依存しない部分の確認
    ret_vector[3]   = abs(im * u_times_R)/norm(R)
    #ret_vector[3]   = abs(im * u_times_R)

    #それ以外の部分の確認
    #ret_vector[2]   = norm(R)*F/geometric_mean^2
    ret_vector[2]   = F/geometric_mean
    
    val             = im * u_times_R * F / geometric_mean^2
    #@show abs(val)
    return abs(val),ret_vector
end



#この表式で計算した要素を行列式の要素として用いることにする． 2021-10-10
function det_stream_no3(I,J,ParamsofHMC,Modelparams,Colorworkspace)
    ρ               = ParamsofHMC.ρ
    θ               = ParamsofHMC.θ
    Nc              = Modelparams.Nc
    u               = Colorworkspace.u
    Sigma           = Colorworkspace.sigma
    Lambda          = Colorworkspace.lambda
    tem             = Colorworkspace.tem
    res             = Colorworkspace.res
    R               = similar(u)
    R_hat           = similar(u)
    prod_su2matrix  = zeros(ComplexF64,2,2)
    prod_matrix     = similar(prod_su2matrix)
    prod_matrix2    = similar(prod_su2matrix)

    ret_vector       = zeros(3)              #返り値を格納するベクトル

    thetaI          = θ[I, 1:Nc]
    thetaJ          = θ[J, 1:Nc]

    λ               = 0.0
    F               = 0.0
    u_times_R       = 0.0
    u_times_R_hat   = 0.0

    #λの構成
    λ               = lambda(I,J,ParamsofHMC)

    #行列のトレース部分の計算
    R               = relative_vector(I,J,ParamsofHMC)
    R_hat           = unit_relative_vector(I,J,ParamsofHMC)
    u               = color_vector(thetaI,thetaJ,Colorworkspace)
    u_times_R       = u[1]*R[1]     + u[2]*R[2]     + u[3]*R[3]     + u[4]*R[4]
    u_times_R_hat   = u[1]*R_hat[1] + u[2]*R_hat[2] + u[3]*R_hat[3] + u[4]*R_hat[4]
    geometric_mean  = sqrt(ρ[I]*ρ[J])
    #@show u_times_R_hat
    #@show 2im * u_times_R_hat

    c1              = 3pi/8
    c2              = (3pi/32)^(4/3)
    F               = c1 * λ^(3/2) / (1+1.25*(λ^2-1) + c2*(λ^2-1)^2)^(3/4)
    ret_vector[1]   = F
    ret_vector[2]   = 0.0
    #@show F

    #↓rhoに依存しない部分の確認
    ret_vector[3]   = abs(im * u_times_R)/norm(R)
    #ret_vector[3]   = abs(im * u_times_R)

    #それ以外の部分の確認
    #ret_vector[2]   = norm(R)*F/geometric_mean^2
    ret_vector[2]   = F/geometric_mean
    
    val             = im * u_times_R_hat * F / geometric_mean
    #@show abs(val)
    #return abs(val),ret_vector
    return abs(val)
    #return ret_vector[2],ret_vector
end

#無次元化済
function det_sum(I,J,ParamsofHMC,Modelparams,Colorworkspace)
    #setting
    ρ               = ParamsofHMC.ρ
    θ               = ParamsofHMC.θ
    Nc              = Modelparams.Nc
    u               = Colorworkspace.u
    R               = similar(u)
    thetaI          = θ[I, 1:Nc]
    thetaJ          = θ[J, 1:Nc]
    u_times_R       = 0.0

    #manipulation
    R               = relative_vector(I,J,ParamsofHMC)
    u               = color_vector(thetaI,thetaJ,Colorworkspace)
    u_times_R       = u[1]*R[1] + u[2]*R[2] + u[3]*R[3] + u[4]*R[4]
    norm_R          = norm(R)
    geometric_mean  = sqrt(ρ[I]*ρ[J])

    val = im * u_times_R * 4.0 / (2.0 + norm_R^2/geometric_mean^2)^2 / geometric_mean^2
    return abs(val)
end

#指定されたフレーバーの質量行列を返す．
function mass_matrix(f,Modelparams)
    #settings
    nI          = Modelparams.nI
    nA          = Modelparams.nA
    mass_diag   = zeros(nI+nA)
    mass_m      = zeros(nI+nA,nI+nA)
    #masses
    m_u         = 3 #Mev
    m_d         = 3 #Mev
    m_s         = 3 #Mev
    m_c         = 4 #Mev
    m_t         = 5 #Mev
    m_b         = 6 #Mev
    mass_vector = [m_u, m_d, m_s, m_c, m_t, m_b]
    
    fill!(mass_diag,mass_vector[f])
    mass_m = diagm(mass_diag)
    return mass_m
end

#指定されたフレーバーの行列式の値を返す．
function dirac_ope_matrix(ParamsofHMC,Modelparams,Colorworkspace)
    #settings
    nI                 = Modelparams.nI
    nA                 = Modelparams.nA
    dirac_matrix_upper = zeros(nI+nA,nI+nA)
    dirac_matrix_lower = zeros(nI+nA,nI+nA)

    #メインの行列計算
        #右上対角成分の計算
        for m = nI+1:nI+nA
            for l = 1:nI
                dirac_matrix_upper[l, m] = det_stream_no3(l,m,ParamsofHMC,Modelparams,Colorworkspace)
            end
        end
        #左下対角成分の計算
        for m = 1:nI
            for l = nI+1:nI+nA
                dirac_matrix_lower[l, m] = det_stream_no3(l,m,ParamsofHMC,Modelparams,Colorworkspace)
            end
        end
    #

    return dirac_matrix_upper + dirac_matrix_lower
end

#全てのフレーバーに対して，det(iD+mf) の積を計算
function det_total_product(total_flavor,ParamsofHMC,Modelparams,Colorworkspace)
    ret_det              = 1.0
    dirac_ope_matrix_tem = dirac_ope_matrix(ParamsofHMC,Modelparams,Colorworkspace) 

    for f = 1:total_flavor
        ret_det *= det( dirac_ope_matrix_tem + mass_matrix(f,Modelparams))
    end
    return ret_det
end

#全てのフレーバーに対して，trlog(iD + mf) の和を計算
function tr_log_fermion(total_flavor,ParamsofHMC,Modelparams,Colorworkspace)
    ret_trlog            = 0.0
    dirac_ope_matrix_tem = dirac_ope_matrix(ParamsofHMC,Modelparams,Colorworkspace) 

    for f = 1:total_flavor
        ret_trlog += tr( log(dirac_ope_matrix_tem + mass_matrix(f,Modelparams)) )
    end
    return abs(ret_trlog)
end

    
        
