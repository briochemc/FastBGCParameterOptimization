#
##===================================================
#Load AIBECS metadata
#===================================================#
#include("../load_packages.jl")
#
#Circulation = OCIM1 # from AIBECS
#
#include("../setup.jl")
#
#path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-2]...)
#
#jld_file = joinpath(path_to_package_root, "data", "Optim_trace.jld2")
#
#@load jld_file results x₀ λ₀ xopt λopt popt μDIPobs xs
#
#
#
#ENV["MPLBACKEND"]="qt5agg"
#using PyPlot, PyCall
#ccrs = pyimport("cartopy.crs")

fig = figure("mismatch", figsize=(12,4))
clf()


#===================================================
Maps
===================================================#
# Projection
central_lon = 200.0

# Add black land
cfeature = pyimport("cartopy.feature")

# lat and lon
lon = vec(grd["xt"])
lat = vec(grd["yt"])

# Contour levels
DIPlevels = 0:0.25:4

# Make the data cyclic
lon_cyc = [lon; 360+lon[1]] # making it cyclic for Cartopy

iz = 10
#i_depths = [1, 15]
nrows = length(i_depths)
#xplots = [μDIPobs, xs[iDIP,1], xs[iDIP,2], xs[iDIP,3], xopt[iDIP]]
xplots = [μDIPobs, xopt[iDIP]]
depth = Int(round(grd["zt"][iz]))
titles = ["Observed DIP at $(depth)m", "Simulated DIP (optimal parameters)"]

ncols = length(xplots)
p = [] # empty array for storing handles of plots

gridspec = PyPlot.matplotlib.gridspec

gs = gridspec.GridSpec(2, 3,
                       width_ratios=[2.3,2.3,1],
                       height_ratios=[10,1]
                       )

for (icol, x) in enumerate(xplots)
    ax = subplot(gs[icol], projection=ccrs.EqualEarth(central_longitude=central_lon))
    ax.add_feature(cfeature.COASTLINE, edgecolor="#000000") # black coast lines
    ax.add_feature(cfeature.LAND, facecolor="#AAAAAA")      # gray land
    ax.gridlines()
    # Convert vector to 3D to 2D
    x3D = NaN * wet3d
    x3D[iwet] = x[iDIP] * ustrip(1.0u"mol/m^3" .|> u"mmol/m^3")
    # make it cyclic for Cartopy
    map_2D = x3D[:,:,iz]
    map_cyc = hcat(map_2D, map_2D[:,1])
    push!(p, PyPlot.contourf(lon_cyc, lat, map_cyc, levels=DIPlevels, transform=ccrs.PlateCarree(), zorder=-1, extend="both"))
    for c in p[end].collections # Removes white line artifact at filled contour edges in PDF
        c.set_edgecolor("face")
        c.set_linewidth(0.1)
    end
    title(titles[icol])
end

#Colorbar
cbar_ax1 = subplot(gs[4])
cbar1 = colorbar(p[1], cax=cbar_ax1, extend="both", orientation="horizontal")
cbar1.set_label("mmol m⁻³")
cbar1.solids.set_edgecolor("face")
cbar1.solids.set_linewidth(0.1)

cbar_ax1bis = subplot(gs[5])
cbar1bis = colorbar(p[2], cax=cbar_ax1bis, extend="both", orientation="horizontal")
cbar1bis.set_label("mmol m⁻³")
cbar1bis.solids.set_edgecolor("face")
cbar1bis.solids.set_linewidth(0.1)

#===================================================
joint PDFs
===================================================#
DIPobs = μDIPobs
using StatsBase, KernelDensity
weights = Weights(v ./ σ²DIPobs / sum(v ./ σ²DIPobs))
lims = [-0.1, 3.6]

cmap = ColorMap("PuBu")

icol = 3
x = xopt[iDIP]
DIPmod = x[iDIP]
bandwidth = (0.00001, 0.00001)
npoints = (2^10, 2^10)
k = kde((DIPobs, DIPmod), weights=weights, bandwidth = bandwidth, npoints = npoints)
I = sortperm(vec(k.density))
dx = k.x[2]-k.x[1]
dy = k.y[2]-k.y[1]
q = k.density[I] * dx * dy
D = zeros(size(k.density))
D[I] .= 100cumsum(q)
levels = 5:5:95
ax = subplot(gs[3])
cmap = ColorMap("cividis")
cmap.set_under([1.0,1.0,1.0])
push!(p, contourf(1e3k.x, 1e3k.y, map(x -> x < 0.05 ? NaN : x, permutedims(D, [2, 1])), levels=levels, cmap=cmap, extend="both"))
for c in p[end].collections # Removes white line artifact at filled contour edges in PDF
    c.set_edgecolor("face")
    c.set_linewidth(0.1)
end
xlim(lims)
ylim(lims)
#ax.axis("equal")
ax.set_aspect("equal", "box")
plot(lims, lims, "r", linewidth=0.5)
title("DIP joint PDF")
xlabel("Observed [mmol m⁻³]")
ylabel("Simulated [mmol m⁻³]")
ax.grid()

cbar_ax2 = subplot(gs[6])
cbar2 = colorbar(p[end], cax=cbar_ax2, extend="both", orientation="horizontal")
cbar2.set_label("Percentile")
cbar2.solids.set_edgecolor("face")
cbar2.solids.set_linewidth(0.1)

eps_file = joinpath(path_to_package_root, "fig", "mismatch.eps")
savefig(eps_file)
png_file = joinpath(path_to_package_root, "fig", "mismatch.png")
savefig(png_file)
svg_file = joinpath(path_to_package_root, "fig", "mismatch.svg")
savefig(svg_file)
pdf_file = joinpath(path_to_package_root, "fig", "mismatch.pdf")
savefig(pdf_file)

