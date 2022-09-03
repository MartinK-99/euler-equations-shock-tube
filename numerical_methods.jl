# flux for Euler equations
function f(u)
    ρ = u[1]
    v = u[2] / ρ
    E = u[3]

    γ = 7/5

    p = (γ-1)*(E - 1/2 * ρ*v^2)

    return [u[2];
            u[2]*v + p;
            v * (E + p)]
end

# initial conditions for shu osher test
function initialCondition(N)
    Δx = 10/(N-1)

    x = 0:Δx:10
    ε = 0.2
    γ = 7/5
    ρ₀(x) = x < 1 ? 3.857153 : 1 + ε * sin(5*x)
    v₀(x) = x < 1 ? 2.629 : 0
    p₀(x) = x < 1 ? 10.333 : 1

    return Δx, x, hcat(ρ₀.(x), ρ₀.(x) .* v₀.(x), p₀.(x) / (γ - 1) + 1/2 * ρ₀.(x) .* v₀.(x) .^ 2)'
end

function laxFriedrichs(N,tn)
    Δx, x, u0 = initialCondition(N)

    Δt = Δx * 0.2

    t = 0:Δt:(tn+Δt) # dependent on the CFL number of the method used
    n = size(t,1)

    u = zeros(Float64,3,N,n+1)
    u[:,:,1] = u0

    for i = 1:n
        # constant boundary condition
        u[:,1,i+1] = u[:,1,i]
        u[:,end,i+1] = u[:,end,i]
        
        for k = 2:N-1
            # Lax-Friedrichs scheme
            u[:,k,i+1] = (u[:,k+1,i] + u[:,k-1,i]) / 2 - Δt / (2 * Δx) * (f(u[:,k+1,i]) - f(u[:,k-1,i]))
        end
    end

    return x, t, u
end

# Godunov flux for local Riemann problems
function fG(ul,ur)
    flux = zeros(Float64,3,1)
    fl = f(ul)
    fr = f(ur)
    f0 = f(zeros(Float64,3))
    for i = 1:size(ul,1)
        if ul[i] > ur[i]
            if ul[i,1] + ur[i,1] >= 0
                flux[i,1] = fl[i,1]
            else
                flux[i,1] = fr[i,1]
            end
        elseif ul[i,1] <= ur[i,1]
            if ul[i,1] >= 0
                flux[i,1] = fl[i,1]
            elseif ul[i,1] <= 0 <= ur[i,1]
                flux[i,1] = f0[i,1]
            else
                flux[i,1] = fr[i,1]
            end
        end
    end
    return flux
end


function godunov(N,tn)
    Δx, x, u0 = initialCondition(N)

    Δt = Δx * 0.2 # dependent on the CFL number of the method used

    t = 0:Δt:(tn+Δt)
    n = size(t,1)
    
    u = zeros(Float64,3,N,n+1)
    u[:,:,1] = u0

    for i = 1:n
        # constant boundary conditions
        u[:,1,i+1] = u[:,1,i]
        u[:,end,i+1] = u[:,end,i]
        
        for k = 2:N-1
            # Godunov scheme
            u[:,k,i+1] = u[:,k,i] + Δt / Δx * (fG(u[:,k-1,i],u[:,k,i]) - fG(u[:,k,i],u[:,k+1,i]))
        end
    end

    return x, t, u
end