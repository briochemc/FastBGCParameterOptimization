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
clf()

fig = figure("mismatch", figsize=(10,15))




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

i_depths = [1, 5, 10, 15]
#i_depths = [1, 15]
nrows = length(i_depths) + 1
#xplots = [μDIPobs, xs[iDIP,1], xs[iDIP,2], xs[iDIP,3], xopt[iDIP]]
xplots = [μDIPobs, xs[iDIP,1], xopt[iDIP]]
ncols = length(xplots)
p = [] # empty array for storing handles of plots

for (irow, iz) in enumerate(i_depths), (icol, x) in enumerate(xplots)
    ax = subplot2grid((6,6), (irow-1, 2*(icol-1)), colspan=2, projection=ccrs.EqualEarth(central_longitude=central_lon))
    ax.add_feature(cfeature.COASTLINE, edgecolor="#000000") # black coast lines
    ax.add_feature(cfeature.LAND, facecolor="#AAAAAA")      # gray land
    # Convert vector to 3D to 2D
    x3D = NaN * wet3d
    x3D[iwet] = x[iDIP] * ustrip(1.0u"mol/m^3" .|> u"mmol/m^3")
    # make it cyclic for Cartopy
    map_2D = x3D[:,:,iz]
    map_cyc = hcat(map_2D, map_2D[:,1])
    push!(p, PyPlot.contourf(lon_cyc, lat, map_cyc, levels=DIPlevels, transform=ccrs.PlateCarree(), zorder=-1, extend="both"))
    for c in p[end].collections # Removes white line artifact at filled contour edges in PDF
        c.set_edgecolor("face")
        c.set_edgewidth(0.1)
    end
    title("$irow, $icol")
end

#Colorbar
fig = gcf()
fig.subplots_adjust(right=0.85, top=0.85)
cbar_ax1 = fig.add_axes([0.9, 0.3, 0.025, 0.65])
cbar1 = fig.colorbar(p[1], cax=cbar_ax1, extend="both")
cbar1.set_label("mmol m⁻³")
cbar1.solids.set_edgecolor("face")


#===================================================
joint PDFs
===================================================#
DIPobs = μDIPobs
using StatsBase, KernelDensity
weights = Weights(v ./ σ²DIPobs / sum(v ./ σ²DIPobs))
lims = [-0.1, 5]

cmap = ColorMap("PuBu")

for (icol, x) in enumerate(xplots)
    icol == 1 && continue
    DIPmod = x[iDIP]
    bandwidth = (0.000005, 0.000005)
    npoints = (2^11, 2^11)
    k = kde((DIPobs, DIPmod), weights=weights, bandwidth = bandwidth, npoints = npoints)
    I = sortperm(vec(k.density))
    dx = k.x[2]-k.x[1]
    dy = k.y[2]-k.y[1]
    q = k.density[I] * dx * dy
    D = zeros(size(k.density))
    D[I] .= cumsum(q)
    levels = 0.1:0.1:0.9
    ax = subplot2grid((6,6), (4, 2*(icol-1)), colspan=2, rowspan=2)
    push!(p, contourf(1e3k.x, 1e3k.y, map(x -> x < 0.01 ? NaN : x, permutedims(D, [2, 1])), levels=levels, cmap="cividis", extend="both"))
    for c in p[end].collections # Removes white line artifact at filled contour edges in PDF
        c.set_edgecolor("face")
        c.set_edgewidth(0.1)
    end
    xlim(lims)
    ylim(lims)
    plot(lims, lims, "r", linewidth=1)
    title("joint PDF")
end

cbar_ax2 = fig.add_axes([0.9, 0.025, 0.025, 0.2])
cbar2 = fig.colorbar(p[end], cax=cbar_ax2, extend="both")
cbar2.set_label("percentile")
cbar2.solids.set_edgecolor("face")

eps_file = joinpath(path_to_package_root, "fig", "mismatch.eps")
savefig(eps_file)
png_file = joinpath(path_to_package_root, "fig", "mismatch.png")
savefig(png_file)
svg_file = joinpath(path_to_package_root, "fig", "mismatch.svg")
savefig(svg_file)
pdf_file = joinpath(path_to_package_root, "fig", "mismatch.pdf")
savefig(pdf_file)

