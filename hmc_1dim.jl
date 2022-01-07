using Plots
using Random
using Distributions
using ForwardDiff
using QuadGK

# Generating Gaussian distribution (1-dim) by using HMC method
function S(v,x,n_size)
    sum(v(x[i]) for i=1:n_size)
end

function dSdx(v,x)
    return ForwardDiff.derivative(v,x)
end

function H(x,p,n_size)
    K = sum(p[i]^2 for i=1:n_size)
    return S(v,x,n_size) + 0.5*K
end

#１変数更新のHMC法
function HMC_1dim(x,parameter)
        n_size,n_tau,τ = parameter
        #初期化
        p_ini = rand(Normal(0,1),n_size)
        p_fin = zeros(n_size)
        x_ini = zeros(n_size)
        x_fin = zeros(n_size)
    
        #x -> x_ini
        x_ini = x[:]
    
        #cal H_ini -> leapfrog -> H_fin
        H_ini = H(x_ini,p_ini,n_size)
        leapfrog!(x,p_ini,parameter)

        x_fin = x[:]
        p_fin = p_ini[:]
        H_fin = H(x_fin,p_fin,n_size)
    
        #Metropolis check
        r=rand()
        if r<exp(H_ini - H_fin)
            x = x_fin[:]
        else
            x = x_ini[:]
        end
end

#x,pをリープフロッグ法で更新する．
function leapfrog!(x,p,parameter)
    n_size,n_tau,τ = parameter
    Δτ  = τ/n_tau 
    p0  = 0.0
    p1  = 0.0
    x05 = 0.0
    x15 = 0.0
    for k=1:n_size
        #initial setting : τ0
        p0  = p[k]
        x05 = x[k] + p0 * 0.5*Δτ
    
        #(n-1)step : τ0 -> τn-1 
        for i=1:n_tau-1
            p1  = p0 - dSdx(v,x05)*Δτ
            x15 = x05 + p1*Δτ
            p0  = p1
            x05 = x15
        end
    
        #final step : τn-1 -> τn
        p_end = p1 - dSdx(v,x15)*Δτ
        x_end = x15 + p_end*0.5*Δτ
        
        x[k]  = x_end
        p[k]  = p_end
    end
end

function gen_configs(Ω,n_con,parameter)
    n_size,n_tau,τ = parameter
    x = zeros(n_size)
    for j=1:n_con
        HMC_1dim(x,parameter)
        Ω[j] = x[1]
    end
end

square(x) = x^2

function x_ev(Ω)
    ret = sum(Ω)
    return ret/length(Ω)
end

function xx_ev(Ω)
    ret = sum(square,Ω)
    return ret/length(Ω)
end

function main()
    #パラメータベクトル
    n_config = 1000000          # # of configurations
    n_size   = 1
    n_tau    = 100
    τ        = 5
    q        = [n_size,n_tau,τ]
    Ω        = zeros(n_config)

    gen_configs(Ω,n_config,q)
    x_value  = x_ev(Ω)
    xx_value = xx_ev(Ω)

    ###plot###
    #gauss(x) = exp(-0.5*x^2)/sqrt(2*pi)
    pdf_func(x) = exp(-v(x))
    Z0,error = quadgk(pdf_func,-Inf,Inf)
    normed_func(x) = pdf_func(x)/Z0

    xr = range(-5,5,length=10001)

    println("<x>  = ",x_value)
    println("<xx> = ",xx_value)
    plot(Ω,st=:histogram,nbins=50,norm=:pdf,alpha=0.3,label="HMC's sample",legend=:topleft)
    plot!(xr,normed_func.(xr),label="pdf_ana")
end

#１次元の関数は適用可能
v(x) = 0.5*x^2
#v(x) = 2*x^2*(x^2-3)
@time main()