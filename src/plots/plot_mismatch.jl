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
#
#using MAT
#OCIM_vars = MAT.matread("/Users/benoitpasquier/.julia/datadeps/OCIM1_MATLAB/CTL.mat")
#OCIM_masks = OCIM_vars["output"]["MSKS"]

fig = figure("mismatch", figsize=(2.7*4,1.5*4))
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
DIPlevels = 0:0.25:3.5
δDIPlevels = -25:2.5:25

# Make the data cyclic
lon_cyc = [lon; 360+lon[1]] # making it cyclic for Cartopy

iz = 12
depth = Int(round(grd["zt"][iz]))
xplots = [μDIPobs, xopt[iDIP]]
#xplots = [μDIPobs, xs[iDIP,1], xs[iDIP,2], xs[iDIP,3], xopt[iDIP]]

p = [] # empty array for storing handles of plots

gridspec = PyPlot.matplotlib.gridspec

gs = gridspec.GridSpec(2, 2,
                       width_ratios=[2,1],
                       height_ratios=[1,1],
                       )

#===================================================
Map 1
===================================================#
ax = subplot(gs[1], projection=ccrs.PlateCarree(central_longitude=central_lon))
ax.add_feature(cfeature.COASTLINE, edgecolor="#000000") # black coast lines
ax.add_feature(cfeature.LAND, facecolor="#EEEEEE")      # gray land
#gl = ax.gridlines(draw_labels=false)
#gl.xlabels_top = false
#gl.ylabels_right = false
#mticker = PyPlot.matplotlib.ticker
#gl.xlocator = mticker.FixedLocator(-180:60:180)
ax.set_xticks(-180:60:179, crs=ccrs.PlateCarree())
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
x3D[iwet] = xopt[iDIP] * ustrip(1.0u"mol/m^3" |> u"mmol/m^3")
# make it cyclic for Cartopy
map_2D = x3D[:,:,iz]
map_cyc = hcat(map_2D, map_2D[:,1])
push!(p, PyPlot.contourf(lon_cyc, lat, map_cyc, cmap="viridis", levels=DIPlevels, transform=ccrs.PlateCarree(), zorder=-1, extend="max"))
for c in p[end].collections # Removes white line artifact at filled contour edges in PDF
    c.set_edgecolor("face")
    c.set_linewidth(0.1)
end
title("(a) Modeled DIP at $(depth)m depth")

#Colorbar
cbar1 = colorbar(p[1], ax=ax, extend="max", orientation="vertical", ticks=DIPlevels[1:2:end], shrink=0.8)
cbar1.set_label("DIP [mmol m⁻³]")
cbar1.solids.set_edgecolor("face")
cbar1.solids.set_linewidth(0.1)

#===================================================
Map 2
===================================================#
ax = subplot(gs[3], projection=ccrs.PlateCarree(central_longitude=central_lon))
ax.add_feature(cfeature.COASTLINE, edgecolor="#000000") # black coast lines
ax.add_feature(cfeature.LAND, facecolor="#EEEEEE")      # gray land
#gl = ax.gridlines(draw_labels=false)
#gl.xlabels_top = false
#gl.ylabels_right = false
#mticker = PyPlot.matplotlib.ticker
#gl.xlocator = mticker.FixedLocator(-180:60:180)
ax.set_xticks(-180:60:179, crs=ccrs.PlateCarree())
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
title("(b) DIP mismatch at $(depth)m depth")

#Colorbar
cbar2 = colorbar(p[end], ax=ax, extend="both", orientation="vertical", ticks=δDIPlevels[1:2:end], shrink=0.8)
cbar2.set_label("δDIP / DIP [%]")
cbar2.solids.set_edgecolor("face")
cbar2.solids.set_linewidth(0.1)


#===================================================
Depth profiles
===================================================#


gspro = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[2])

# basins from OCIM MAT file
basins = ["ATL", "PAC", "IND"]

for i in 1:3
    ax = subplot(gspro[i])
    mask3d = OCIM_masks[basins[i]]
    x3D = 0 * wet3d

    x3D[iwet] = xopt[iDIP] * ustrip(1.0u"mol/m^3" |> u"mmol/m^3")
    profile_mod = [mean(x3D[:,:,iz], StatsBase.weights(array_of_volumes(grd)[:,:,iz] .* mask3d[:,:,iz])) for iz in 1:24]

    #x3D[iwet] = xs[iDIP,1] * ustrip(1.0u"mol/m^3" |> u"mmol/m^3")
    #profile_x0 = [mean(x3D[:,:,iz], StatsBase.weights(array_of_volumes(grd)[:,:,iz] .* mask3d[:,:,iz])) for iz in 1:24]

    x3D[iwet] = μDIPobs * ustrip(1.0u"mol/m^3" |> u"mmol/m^3")
    profile_obs = [mean(x3D[:,:,iz], StatsBase.weights(array_of_volumes(grd)[:,:,iz] .* mask3d[:,:,iz])) for iz in 1:24]

    z = vec(grd["zt"])

    #push!(p, plot(profile_x0, z))
    push!(p, plot(profile_obs, z))
    push!(p, plot(profile_mod, z))
    ax.invert_yaxis()
    xticks(0:1:2.5)
    xlim([0, 2.8])
    if i == 1
        ylabel("depth [m]")
    else
        ax.set_yticks([])
    end

    if i == 2
        title("(c) Ocean-basin mean DIP depth profiles")
        xlabel("DIP [mmol m⁻³]")
        #legend(("x0", "opt", "obs"), handlelength=0.8)
        legend(("mod", "obs"), handlelength=0.8)
    end

    ax_c = ax.twiny()
    ax_c.set_xlabel(basins[i])
    ax_c.set_xticks([])
end

#title(gs[2],"(c) DIP depth profile")
#ax.axis("equal")

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
ax = subplot(gs[4])
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
title("(d) Modeled–observed DIP joint PDF")
xlabel("observed DIP [mmol m⁻³]")
ylabel("modeled DIP [mmol m⁻³]")
#ax.grid(linestyle=":")

cbar3 = colorbar(p[end], ax=ax, extend="both", orientation="vertical", ticks=levels[1:2:end], shrink=0.8)
cbar3.set_label("percentile")
cbar3.solids.set_edgecolor("face")
cbar3.solids.set_linewidth(0.1)

PyPlot.tight_layout()

#eps_file = joinpath(path_to_package_root, "fig", "mismatch.eps")
#savefig(eps_file)
#png_file = joinpath(path_to_package_root, "fig", "mismatch.png")
#savefig(png_file)
#svg_file = joinpath(path_to_package_root, "fig", "mismatch.svg")
#savefig(svg_file)
pdf_file = joinpath(path_to_package_root, "fig", "mismatch.pdf")
savefig(pdf_file, bbox_inches="tight")

