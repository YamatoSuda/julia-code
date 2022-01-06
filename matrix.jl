#matrix.jl

    #'tHooft symbol
    macro swap!(x,y)
        quote
        local tmp = $(esc(x))
        $(esc(x)) = $(esc(y))
        $(esc(y)) = tmp
        end
    end
    
    # τ+ 行列
    function tau_p(μ)
        τ = zeros(ComplexF64,2,2)
        if μ==0
            println("mu of tau is 0!")
            return
        elseif μ > 4
            println("mu of tau is lager than 4!")
            return
        end
        
        if μ == 4
            τ[1,1] = -im
            τ[2,2] = -im
            return τ
        else
            return sigma(μ)
        end
    end
    
    
    # τ- 行列
    function tau_n(μ)
        τ = zeros(ComplexF64,2,2)
        if μ==0
            println("mu of tau is 0!")
            return
        elseif μ > 4
            println("mu of tau is lager than 4!")
            return
        end
        
        if μ == 4
            τ[1,1] = im
            τ[2,2] = im
            return τ
        else
            return sigma(μ)
        end
    end
    
    # Edinton epsilon symbol
    function ϵ(a,b,c)
        if a > 3 || b > 3 || c > 3
            println("a or b or c is bigger than 3!") 
            return 
        end
        
        if a == b || a == c || b == c
            return 0
        end
        
        #sorting 
        count = 0
        dum   = 0
        n     = 2 # (引数の個数)-1
        for i = 1:n
            if a > b 
                @swap!(a,b)
                count += 1
            elseif b > c
                @swap!(b,c)
                count += 1
            end
        end
        return (-1)^count
        
    end
    
    # Kronecker delta symbpl
    function δ(a,b)
        a == b ? 1 : 0
    end
    
    # 'tHooft eta symbol
    function η(a,μ,ν)
        if ν == 4
            return δ(a,μ)
        elseif μ == 4
            return -δ(a,ν)
        else
            return ϵ(a,μ,ν)
        end
    end
    
    # 'tHooft eta_bar symbol
    function ηbar(a,μ,ν)
        if ν == 4
            return -δ(a,μ)
        elseif μ == 4
            return δ(a,ν)
        else
            return ϵ(a,μ,ν)
        end
    end
    
    # Pauli matrix
    function sigma(i)
        if i < 1 || i > 3 return println("In Pauli matrix,i is not defined!") end
        sig = zeros(ComplexF64,2,2)
        if i == 1
            sig[1,2] = 1   
            sig[2,1] = 1
            return sig 
        elseif i == 2
            sig[1,2] = -im 
            sig[2,1] =  im
            return sig
        else 
            sig[1,1] = 1   
            sig[2,2] = -1
            return sig
        end
    end
    
    # Gell-Mann matrices
    function lambda(a)
        if a < 1 || a > 8 return println("In Gell-Mann matrix,a is not defined!") end
        lam = zeros(ComplexF64,3,3)
        C8  = 1/sqrt(3)
        if a == 1
            lam[1,2] = 1
            lam[2,1] = 1
            return lam
        elseif a == 2
            lam[1,2] = -im
            lam[2,1] =  im
            return lam
        elseif a == 3
            lam[1,1] =  1
            lam[2,2] = -1
            return lam
        elseif a == 4
            lam[1,3] =  1
            lam[3,1] = -1
            return lam
        elseif a == 5
            lam[1,3] = -im
            lam[3,1] =  im
            return lam
        elseif a == 6 
            lam[2,3] = 1
            lam[3,2] = 1
            return lam
        elseif a == 7
            lam[2,3] = -im
            lam[3,2] =  im
            return lam
        else 
            lam[1,1] =  1
            lam[2,2] =  1
            lam[3,3] = -2
            return lam * C8
        end
    end

    function make_sigma()
        any = []
        for i=1:3
            push!(any,sigma(i))
        end
        return any
    end

    function make_lambda()
        any = []
        for a=1:8
            push!(any,lambda(a))
        end
        return any
    end

    # SU(2) color orientation matrix (random)
    function su2_matrix(su2,sig)
        return exp(im * 0.5 * sum(su2 .* sig))
    end

    # SU(3) color orientation matrix (random)
    function su3_matrix(su3,lam)
        return exp(im * 0.5 * sum(su3 .* lam))
    end

