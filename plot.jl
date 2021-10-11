#plot用の関数を作成するファイル

function plot_func(lengthx,lwidth; func="sin", color="red")
    if func == "sin"
        xr = range(-lengthx,lengthx,length=1000)
    plot(xr,sin.(xr))
    elseif func == "cos"
        xr = range(-lengthx,lengthx,length=1000)
    plot(xr,cos.(xr))
    end
end

function plot_overlap_emlements(Rstart, Rend, rho1, rho2, ParamsofHMC, Modelparams, Colorworkspace; type="sl", re=0, sample=1000)
    plot_R_scaled   = range(Rstart, Rend, length=sample) #無次元化した R' = R/ρ に対応
    plot_Tia        = similar(plot_R_scaled)
    plotname        = ""

    #プロットする領域だけ変数を設定
    ParamsofHMC.z[1,1:4] .= 0.0
    ParamsofHMC.z[2,1:4] .= 0.0
    ParamsofHMC.θ[1,1:3] .= 1.0
    ParamsofHMC.θ[2,1:3] .= 1.0
    ParamsofHMC.ρ[1]      = rho1
    ParamsofHMC.ρ[2]      = rho2

    #scaling factor
    geometric_mean        = sqrt(rho1 * rho2)
    @show geometric_mean

    #scalingを行う
    ParamsofHMC.ρ[1]      = ParamsofHMC.ρ[1] / geometric_mean
    ParamsofHMC.ρ[2]      = ParamsofHMC.ρ[2] / geometric_mean

    #scalingされた変数
    @show ParamsofHMC.ρ[1] 
    @show ParamsofHMC.ρ[2]
    
    if type == "sl1"
        plotname = "streamlineNo1"
        for i in eachindex(plot_R_scaled)
            ParamsofHMC.z[2,4]     = plot_R_scaled[i]
            plot_Tia[i],ret_vector = det_stream_no1(1,2, ParamsofHMC, Modelparams, Colorworkspace)
            #@show ret_vector[1]
        end
    elseif type == "sl2"
        plotname = "streamlineNo2"
        for i in eachindex(plot_R_scaled)
            ParamsofHMC.z[2,4]     = plot_R_scaled[i]
            #@show ret_vector[1]
            plot_Tia[i],ret_vector = det_stream_no2(1,2, ParamsofHMC, Modelparams, Colorworkspace)
        end
    elseif type == "sl3"
        plotname = "streamlineNo3"
        for i in eachindex(plot_R_scaled)
            ParamsofHMC.z[2,4]     = plot_R_scaled[i]
            #@show ret_vector[1]
            plot_Tia[i],ret_vector = det_stream_no3(1,2, ParamsofHMC, Modelparams, Colorworkspace)
        end
    elseif type == "sum"
        plotname = "sum ansatz     "
        for i in eachindex(plot_R_scaled)
            ParamsofHMC.z[2,4] = plot_R_scaled[i] # R'=R/ρ の値で 1~5 の範囲を計算するので，plot_Rを用いてよい．
            plot_Tia[i]        = det_sum(1,2, ParamsofHMC, Modelparams, Colorworkspace)
        end
    end

    if re == 0
        plot(plot_R_scaled, plot_Tia,
        xlim=(Rstart,Rend), ylim=(0,1.3),
        lw=2, lc=:blue, alpha=0.7,
        xlabel="R/ρ", ylabel="T(1/ρ)",
        label="$plotname [ρ1=$rho1, ρ2=$rho2, ρ=$geometric_mean]")
    elseif re == 1
        plot!(plot_R_scaled, plot_Tia,
        lw=2, lc=:red, alpha=0.7,
        label="$plotname [ρ1=$rho1, ρ2=$rho2, ρ=$geometric_mean]")
    elseif re == 2
        plot!(plot_R_scaled, plot_Tia,
        lw=2, lc=:magenta, alpha=0.7,
        label="$plotname [ρ1=$rho1, ρ2=$rho2, ρ=$geometric_mean]")
    elseif re == 3
        plot!(plot_R_scaled, plot_Tia,
        lw=2, lc=:cyan, alpha=0.7,
        label="$plotname [ρ1=$rho1, ρ2=$rho2, ρ=$geometric_mean]")
    end
end


#============================================================================================#
#============================================================================================#

function plot_overlap_integral(x_min, x_max; type="F", re = "false", sample=1000)

    plot_r      = range(x_min,x_max, length=sample)
    plot_func   = similar(plot_r)
    label       = ""

    if type == "F"
        integ_func(x,a) = 6.0 * x^(3/2) / ((x + 1 / a)^(3/2)*(x + a)^(5/2))
        for i in eachindex(plot_r)
            F, error     = det_integral(integ_func,plot_r[i])   
            #@show plot_r[i]    
            plot_func[i] = F
        end
        label = "F(λ)"
    else
        c1  = 3pi/8
        c2  = (3pi/32)^(4/3)
        for i in eachindex(plot_r)
            λ            = plot_r[i]
            #@show λ
            plot_func[i] = c1 * λ^(3/2) / (1+1.25*(λ^2-1) + c2*(λ^2-1)^2)^(3/4)
        end
        label = "G(λ)"
    end

    if re == "false"
        plot(plot_r, plot_func, 
            xlabel="λ", ylabel="$label",label="$label",title="Overlap integral",
            lw=2, lc=:blue
        )
    else 
        plot!(plot_r, plot_func, 
            label="$label",
            lw=2, lc=:red
        )
    end
end

#============================================================================================#
#============================================================================================#