using Plots

print("Hellow world")

#=
x -> x^2
のような記法は，xをx^2に対応させるというニュアンスで用いることができる．
それ以降には，プロットのオプションを列挙している．
=#
plot(
    x->x^3+x^2+x, 
    xlabel="x",             
    ylabel="y",
    label="y=x^3+x^2+x",
    xlim=(-3,3),
    ylim=(-10,10),
    xticks=-3:0.5:3,
    title="Test plot",
    linecolor=:blue,
    linewidth=3,
    linestyle=:dot,
)

#パラメータを含むような関数を書くこともできる．
a=1.0
plot(
    x->(sin(a*x)/(a*x)-cos(a*x))^2,
    xlim=(-10,10)
)


#===============================================================#
x = y = range(-5,5, length = 40)
anime =  @animate for i=-100:250
κ = -i*0.1
λ = 1.0
f(x,y) = κ*(x^2 + y^2) + λ*((x^2 + y^2)^2)

plt = surface(x,y,f,
        zlim = (-100,30),
        transparency = true,
        alpha = 0.5,
        colorbar = false,
        camera = (15,80),
        resolution = (1000,1000),
        c = :thermal
) 
    plot!(plt[1])
end

gif(anime, "./sym_break.gif");