
#include("../load_packages.jl")
#
#Circulation = OCIM1 # from AIBECS
#
#include("../src/setup.jl")
#
#path_to_package_root = joinpath(splitpath(@__DIR__)[1:end-2]...)
#
#jld_file = joinpath(path_to_package_root, "data", "Optim_trace.jld2")
#
#@load jld_file results x₀ λ₀ xopt λopt popt μDIPobs
#
ENV["MPLBACKEND"]="qt5agg"
using PyPlot, PyCall
ccrs = pyimport("cartopy.crs")
clf()



#===================================================
Observations map
===================================================#
# Projection
central_lon = 200.0
ax1 = subplot(211, projection=ccrs.EqualEarth(central_longitude=central_lon))

# Add black land
cfeature = pyimport("cartopy.feature")
ax1.add_feature(cfeature.COASTLINE, edgecolor="#000000") # black coast lines
ax1.add_feature(cfeature.LAND, facecolor="#AAAAAA")      # gray land

# lat and lon
lon = vec(grd["xt"])
lat = vec(grd["yt"])

# Convert vector to 3D to 2D
DIPobs3D = NaN * wet3d
DIPobs3D[iwet] = μDIPobs * ustrip(1.0u"mol/m^3" .|> u"mmol/m^3")
iz = 7
depth = round(grd["zt"][iz])
map_2D1 = DIPobs3D[:,:,iz]

# Contour levels
DIPlevels = 0:0.2:4

# Make the data cyclic
lon_cyc = [lon; 360+lon[1]] # making it cyclic for Cartopy
map_cyc1 = hcat(map_2D1, map_2D1[:,1])

p1 = PyPlot.contourf(lon_cyc, lat, map_cyc1, levels=DIPlevels, transform=ccrs.PlateCarree(), zorder=-1, extend="max")

title("Observed DIP (at $(depth)m)")


#===================================================
Optimal solution map
===================================================#
ax2 = subplot(212, projection=ccrs.EqualEarth(central_longitude=central_lon))
# Add black land
cfeature = pyimport("cartopy.feature")
ax2.add_feature(cfeature.COASTLINE, edgecolor="#000000") # black coast lines
ax2.add_feature(cfeature.LAND, facecolor="#AAAAAA")      # gray land
# Convert vector to 3D to 2D
DIPopt3D = NaN * wet3d
DIPopt3D[iwet] = xopt[iDIP] * ustrip(1.0u"mol/m^3" .|> u"mmol/m^3")
map_2D2 = DIPopt3D[:,:,iz]

map_cyc2 = hcat(map_2D2, map_2D2[:,1])

p2 = PyPlot.contourf(lon_cyc, lat, map_cyc2, levels=DIPlevels, transform=ccrs.PlateCarree(), zorder=-1, extend="max")
title("Simulated DIP (after optimization)")



#===================================================
Colorbar
===================================================#
fig = gcf()
fig.subplots_adjust(right=0.86)
cbar_ax = fig.add_axes([0.81, 0.15, 0.025, 0.7])
cbar = fig.colorbar(p1, cax=cbar_ax, extend="max")
cbar.set_label("mmol m⁻³")


#===================================================
joint PDF?
===================================================#

#ax3 = subplot(212, projection=ccrs.EqualEarth(central_longitude=central_lon))
## Add black land
#cfeature = pyimport("cartopy.feature")
#ax2.add_feature(cfeature.COASTLINE, edgecolor="#000000") # black coast lines
#ax2.add_feature(cfeature.LAND, facecolor="#AAAAAA")      # gray land
## Convert vector to 3D to 2D
#DIPopt3D = NaN * wet3d
#DIPopt3D[iwet] = xopt[iDIP] * ustrip(1.0u"mol/m^3" .|> u"mmol/m^3")
#map_2D2 = DIPopt3D[:,:,iz]
#
#map_cyc2 = hcat(map_2D2, map_2D2[:,1])
#
#p2 = PyPlot.contourf(lon_cyc, lat, map_cyc2, levels=DIPlevels, transform=ccrs.PlateCarree(), zorder=-1, extend="max")
#title("Optimal simulated DIP")

DIPopt = xopt[iDIP]
k = kde((μDIPobs, DIPopt), bandwidth = (0.00001, 0.00001))
I = sortperm(vec(k.density))
q = k.density[I] * dx * dy
D = zeros(size(k.density))
D[I] .= cumsum(q)
levels = 0.05:0.05:0.95
contourf(1e3k.x, 1e3k.y, permutedims(D, [2, 1]), levels=levels)
colorbar
