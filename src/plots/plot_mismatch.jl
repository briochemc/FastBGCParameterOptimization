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

fig = figure("mismatch", figsize=(2.7*4,4))
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
δDIPlevels = -25:2.5:25

# Make the data cyclic
lon_cyc = [lon; 360+lon[1]] # making it cyclic for Cartopy

iz = 12
depth = Int(round(grd["zt"][iz]))
xplots = [μDIPobs, xopt[iDIP]]

p = [] # empty array for storing handles of plots

gridspec = PyPlot.matplotlib.gridspec

gs = gridspec.GridSpec(1, 2,
                       width_ratios=[2,1],
                       )

ax = subplot(gs[1], projection=ccrs.PlateCarree(central_longitude=central_lon))
ax.add_feature(cfeature.COASTLINE, edgecolor="#000000") # black coast lines
ax.add_feature(cfeature.LAND, facecolor="#EEEEEE")      # gray land
#gl = ax.gridlines(draw_labels=false)
#gl.xlabels_top = false
#gl.ylabels_right = false
#mticker = PyPlot.matplotlib.ticker
#gl.xlocator = mticker.FixedLocator(-180:60:180)
ax.set_xticks(-180:60:180, crs=ccrs.PlateCarree())
ax.set_yticks(-80:20:80, crs=ccrs.PlateCarree())
cticker = pyimport("cartopy.mpl.ticker")
LongitudeFormatter = cticker.LongitudeFormatter
LatitudeFormatter = cticker.LatitudeFormatter
lon_formatter = LongitudeFormatter(
                                   degree_symbol="",
                                   dateline_direction_label=true)
lat_formatter = LatitudeFormatter(
                                      degree_symbol="")
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
# Convert vector to 3D to 2D
x3D = NaN * wet3d
x3D[iwet] = 100(xopt[iDIP] - μDIPobs) ./ μDIPobs
# make it cyclic for Cartopy
map_2D = x3D[:,:,iz]
map_cyc = hcat(map_2D, map_2D[:,1])
push!(p, PyPlot.contourf(lon_cyc, lat, map_cyc, cmap="PiYG_r", levels=δDIPlevels, transform=ccrs.PlateCarree(), zorder=-1, extend="both"))
for c in p[end].collections # Removes white line artifact at filled contour edges in PDF
    c.set_edgecolor("face")
    c.set_linewidth(0.1)
end
title("(a) DIP mismatch at $(depth)m depth")

#Colorbar
cbar1 = colorbar(p[1], ax=ax, extend="both", orientation="vertical", ticks=δDIPlevels[1:2:end], shrink=0.8)
cbar1.set_label("δDIP / DIP [%]")
cbar1.solids.set_edgecolor("face")
cbar1.solids.set_linewidth(0.1)

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
bandwidth = (0.000005, 0.000005)
npoints = (2^10, 2^10)
k = kde((DIPobs, DIPmod), weights=weights, bandwidth = bandwidth, npoints = npoints)
I = sortperm(vec(k.density))
dx = k.x[2]-k.x[1]
dy = k.y[2]-k.y[1]
q = k.density[I] * dx * dy
D = zeros(size(k.density))
D[I] .= 100cumsum(q)
levels = 5:5:95
ax = subplot(gs[2])
cmap = ColorMap("magma_r")
#cmap.set_under([1.0,1.0,1.0]) # White color for values below minimum
push!(p, contourf(1e3k.x, 1e3k.y, map(x -> x ≤ 0.1 ? NaN : x, permutedims(D, [2, 1])), levels=levels, cmap=cmap, extend="both"))
for c in p[end].collections # Removes white line artifact at filled contour edges in PDF
    c.set_edgecolor("face")
    c.set_linewidth(0.1)
end
xlim(lims)
ylim(lims)
xticks(0:1:3)
yticks(0:1:3)
#ax.axis("equal")
ax.set_aspect("equal", "box")
plot(lims, lims, "k:", linewidth=0.5)
title("(b) DIP joint PDF")
xlabel("Observed DIP [mmol m⁻³]")
ylabel("Simulated DIP [mmol m⁻³]")
#ax.grid(linestyle=":")

cbar2 = colorbar(p[end], ax=ax, extend="both", orientation="vertical", ticks=levels[1:2:end], shrink=0.8)
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
savefig(pdf_file, bbox_inches="tight")

