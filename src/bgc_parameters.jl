# Build the parameters type and p₀
t = empty_parameter_table()    # initialize table of parameters
add_parameter!(t, :xgeo, 2.17u"mmol/m^3",
    optimizable = true,
    description = "Geological mean P concentration",
    LaTeX = "\\state^\\mathrm{geo}")
add_parameter!(t, :τg, 1.0u"Myr",
    description = "Geological restoring timescale",
    LaTeX = "\\tau_\\mathrm{geo}")
add_parameter!(t, :ku, 10.0u"μmol/m^3",
    optimizable = true,
    description = "Half-saturation constant (Michaelis-Menten)",
    LaTeX = "k_\\vec{u}")
add_parameter!(t, :z₀, 80.0u"m",
    description = "Depth of the euphotic layer base",
    LaTeX = "z_0")
add_parameter!(t, :w₀, 1.0u"m/d",
    optimizable = true,
    description = "Sinking velocity at surface",
    LaTeX = "w_0")
add_parameter!(t, :w′, 1/4.4625u"d",
    optimizable = true,
    description = "Vertical gradient of sinking velocity",
    LaTeX = "w'")
add_parameter!(t, :κ, 1/5.25u"d",
    optimizable = true,
    description = "Dissolution rate constant (POP to DIP)",
    LaTeX = "\\kappa")
add_parameter!(t, :τu, 30.0u"d",
    optimizable = true,
    description = "Maximum uptake rate timescale",
    LaTeX = "\\tau_\\vec{u}")
initialize_Parameters_type(t)   # Generate the parameter type

