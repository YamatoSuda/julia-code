#fermionDets.jl

function DetIntegralold(λ)
    DetIntegrand(x) = 6.0 * x^(3/2) / ((x + 1 / λ)^(3/2)*(x + λ)^(5/2))
    ret, err = quadgk(DetIntegrand,0.0,Inf)
end

function DetIntegral(f,λ)
    DetIntegrand(x) = f(x,λ)
    ret, err = quadgk(DetIntegrand,0.0,Inf)
end

function DetStreamlineNo1(I,J,dynamical_p,dynamical_q)
    ρ, z, θ, A, g       = dynamical_p
    u, Sigma, tem, res  = dynamical_q
    Rhat                = similar(u)
    TautimesR           = zeros(ComplexF64,2,2)
    ProdOfSU2matrix     = similar(TautimesR)
    ProdOfMatrix        = similar(TautimesR)
    ProdOfMatrix2       = similar(TautimesR)

    retVector           = zeros(3)
    
    ParamOfSU2I         = θ[I, 1:3]
    ParamOfSU2J         = θ[J, 1:3]

    λ               = 0.0
    F               = 0.0
    error           = 0.0
    TraceOfMatrix   = 0.0
    TraceOfMatrix2  = 0.0
    uprodR          = Complex(0,0)
            
    # λの構成，F(λ)の計算
    λ               = lambda(I,J,dynamical_p)
    IntegFunc(x,a)  = 6.0 * x^(3/2) / ((x + 1 / a)^(3/2)*(x + a)^(5/2))
    F, error        = DetIntegral(IntegFunc,λ)                              #F(λ) = 6.0 \int f(x,a) dx を返り値としている
    retVector[1]    = F
    retVector[2]    = error
    #@show λ
    #@show F
    #F, error        = DetIntegralold(λ)
    
    # 行列トレース部分の計算
    u               = color_vector(ParamOfSU2I,ParamOfSU2J,dynamical_q)
    R_hat           = relative_vector(I,J,z)
    TautimesR       = tau_p(1)*R_hat[1] + tau_p(2)*R_hat[2] + tau_p(3)*R_hat[3] + tau_p(4)*R_hat[4]
    mul!(ProdOfSU2matrix, SU2_matrix(ParamOfSU2I,Sigma)', SU2_matrix(ParamOfSU2J,Sigma))
    mul!(ProdOfMatrix, TautimesR, ProdOfSU2matrix)
    TraceOfMatrix   = tr(ProdOfMatrix)
    geoMean         = sqrt(ρ[I]*ρ[J])
    #TraceOfMatrix2  = tr(ProdOfSU2matrix*tau_p(1)) * R_hat[1] + tr(ProdOfSU2matrix*tau_p(2)) * R_hat[2] + tr(ProdOfSU2matrix*tau_p(3)) * R_hat[3] + tr(ProdOfSU2matrix*tau_p(4)) * R_hat[4] 

    retVector[3]    = abs(0.5 * TraceOfMatrix)

    val = 0.5 * TraceOfMatrix * F / geoMean^2
    #val = 0.5 * TraceOfMatrix * F 
    #@show abs(val)
    return abs(val),retVector
    #return val
end

function DetStreamlineNo2(I,J,dynamical_p,dynamical_q)
    ρ, z, θ, A, g       = dynamical_p
    u, Sigma, tem, res  = dynamical_q
    R                   = similar(u)
    R_hat               = similar(u)
    
    ParamOfSU2I         = θ[I, 1:3]
    ParamOfSU2J         = θ[J, 1:3]

    retVector           = zeros(3)

    λ               = 0.0
    F               = 0.0
    utimesR         = 0.0
    utimesR_hat     = 0.0
    
    # λの構成
    λ               = lambda(I,J,dynamical_p)
    
    # 行列トレース部分の計算
    R               = space_ori_vector(I,J,z)
    R_hat           = relative_vector(I,J,z)
    u               = color_vector(ParamOfSU2I,ParamOfSU2J,dynamical_q)
    utimesR         = u[1]*R[1]     + u[2]*R[2]     + u[3]*R[3]     + u[4]*R[4]
    utimesR_hat     = u[1]*R_hat[1] + u[2]*R_hat[2] + u[3]*R_hat[3] + u[4]*R_hat[4]
    geoMean         = sqrt(ρ[I]*ρ[J])
    #@show utimesR_hat
    #@show 2im * utimesR_hat

    c1              = 3pi/8
    c2              = (3pi/32)^(4/3)
    F               = c1 * λ^(3/2) / (1+1.25*(λ^2-1) + c2*(λ^2-1)^2)^(3/4)
    retVector[1]    = F
    retVector[2]    = 0.0
    #@show F
    retVector[3]    =  abs(im * utimesR_hat)
    val = im * utimesR_hat * F / geoMean^2
    #val = im * utimesR_hat * F 
    #val = im * utimesR * F / (ρ[I]*ρ[J])
    #@show abs(val)
    return abs(val),retVector
    #return val
end

function DetSumAnsatz(I,J,dynamical_p,dynamical_q)
    ρ, z, θ, A, g       = dynamical_p
    u, Sigma, tem, res  = dynamical_q
    R                   = similar(u)
    R_hat               = similar(u)
    
    ParamOfSU2I         = θ[I, 1:3]
    ParamOfSU2J         = θ[J, 1:3]

    utimesR         = 0.0
    utimesR_hat     = 0.0
    
    R               = space_ori_vector(I,J,z)
    R_hat           = relative_vector(I,J,z)
    u               = color_vector(ParamOfSU2I,ParamOfSU2J,dynamical_q)
    utimesR         = u[1]*R[1] + u[2]*R[2] + u[3]*R[3] + u[4]*R[4]
    utimesR_hat     = u[1]*R_hat[1] + u[2]*R_hat[2] + u[3]*R_hat[3] + u[4]*R_hat[4]
    NormR           = norm(R)
    geoMean         = sqrt(ρ[I]*ρ[J])
    #@show utimesR_hat
    #@show 2im * utimesR_hat

    #val = im * utimesR_hat * 4.0 / ( geoMean^2 * (2.0 + NormR^2/geoMean^2)^2 )
    val = im * utimesR * 4.0 / (2.0 + NormR^2/geoMean^2)^2 / geoMean^2
    #val = im * utimesR * 4.0 / (2.0 + NormR^2/geoMean^2)^2 
    #@show val
    #@show abs(val)
    return abs(val)
    #return val
end

#fermionの行列式を計算するには，クォークのフレーばごとに質量行列を考える必要がある．
function DetDiracOpe()

end



    
        
