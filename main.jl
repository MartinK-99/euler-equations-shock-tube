include("numerical_methods.jl")

function main()
    # number of cells in discretization
    N = 1000
    
    # duration of simulation
    tn = 1.8

    @time x,t,uLF = laxFriedrichs(N,tn)
    @time x,t,uGod = godunov(N,tn)

    #@time plottingSolution(x,t,uGod,"Godunov Scheme")
    @time plottingTwoSolutions(x,t,uLF,uGod,["Lax-Friedrichs' Scheme","Godunov's Scheme"])
end

function plottingSolution(x,t,u,labelstring)
    skip = maximum([1, Int(floor(size(u,3)/200))])
    println("Start animating")
    res = maximum([1,round(Int,size(u,2)/500)])
    titel = "Shu Osher Test"
    anim = @animate for i ∈ append!(ones(Int,25), 1:skip:(size(t,1)),(size(t,1))*ones(Int,50))
        plot(x[1:res:end],u[1,1:res:end,i],ylim=(-0.5,5),label = labelstring,legend=:bottomleft)
        title!(titel * ", N = " * string(size(u,2)) * ", t = " * rpad(string(floor(t[i],digits=2)),4,"0"))
        xlabel!("x")
        ylabel!("ρ")
    end
    gif(anim, "animation.gif", fps = 25)
end

function plottingTwoSolutions(x,t,u1,u2,labels)
    skip = maximum([1, Int(floor(size(u1,3)/200))])
    println("Start animating for comparison")
    res = maximum([1,round(Int,size(u1,2)/500)])
    titel = "Shu Osher Test"

    anim = @animate for i ∈ append!(ones(Int,25), 1:skip:(size(t,1)),(size(t,1))*ones(Int,50))
        plot(x[1:res:end],u1[1,1:res:end,i],ylim=(-0.5,5),label = labels[1],legend=:bottomleft)
        plot!(x[1:res:end],u2[1,1:res:end,i],ylim=(-0.5,5),label = labels[2],legend=:bottomleft)
        title!(titel * ", N = " * string(size(u1,2)) * ", t = " * rpad(string(floor(t[i],digits=2)),4,"0"))
        xlabel!("x")
        ylabel!("ρ")
    end
    gif(anim, "animationComparison.gif", fps = 25)
end

main()