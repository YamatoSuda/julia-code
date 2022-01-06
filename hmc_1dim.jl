using Plots
using Random
using Distributions
using ForwardDiff

# Generating Gaussian distribution (1-dim) by using HMC method
v(x) = 0.5*x^2

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
    
        #初期配位を保存
        for i=1:n_size x_ini[i] = x[i] end
    
        #初期配位のハミルトニアンを計算
        H_ini = H(x_ini,p_ini,n_size)
    
        #リープフロッグで時間発展
        leapfrog!(x,p_ini,parameter)
    
        #時間発展後の配位を保存
        for i=1:n_size x_fin[i] = x[i]; p_fin[i] = p_ini[i] end
    
        #時間発展後のハミルトニアンを計算
        H_fin = H(x_fin,p_fin,n_size)
    
        #Metropolis check
        r=rand()
        if r<exp(H_ini - H_fin)
            for i=1:n_size x[i] = x_fin[i] end
        else
            for i=1:n_size x[i] = x_ini[i] end
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

function main()
    #パラメータベクトル
    n_config = 1_000_000          # # of configurations
    n_size   = 1
    n_tau    = 100
    τ        = 5
    q        = [n_size,n_tau,τ]
    Ω        = zeros(n_config)

    #result = @benchmark gen_configs(n_config,q)
    @time gen_configs(Ω,n_config,q)

    ##########################################################

    ###plot###
    gauss(x) = exp(-0.5*x^2)/sqrt(2*pi)
    #pdf_v(x) = exp(-v(x))
    #Z0,error = quadgk(pdf_v,-Inf,Inf)
    #normed_v(x) = pdf_v(x)/Z0

    xr = range(-5,5,length=10001)

    plot(Ω,st=:histogram,nbins=50,norm=:pdf,alpha=0.3,label="HMC's sample",legend=:topleft)
    plot!(xr,gauss.(xr),label="pdf_ana")
end
main()