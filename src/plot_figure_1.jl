# Plots the trace of the optimizer in 3d parameter space

# Initialize the plot
using Plots
plt1 = plot3d(1, title = "Optimizer trace", marker = 2)

for λ in Optim.x_trace(results)
    push!(plt1, λ...)
end

display(plt1)