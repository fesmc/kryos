# A miniature yelmox checkout for the converter tests: the files runme and the
# converter read (.runme/info.json, the defaults and library files of input/,
# one par file), small enough to read and complete enough for the whole path.

const YELMO_GROUPS = ["ytopo", "ycalv", "ydyn", "ytill", "yhyd", "ymat", "ytrc", "ytherm",
                      "yelmo_masks", "yelmo_init_topo", "yelmo_data"]
const YELMO_NML = Dict("ytopo" => "nml_ytopo", "ycalv" => "nml_ycalv", "ydyn" => "nml_ydyn", "ytill" => "nml_ytill",
                       "yhyd" => "nml_yhyd", "ymat" => "nml_ymat", "ytrc" => "nml_ytrc", "ytherm" => "nml_ytherm",
                       "yelmo_masks" => "nml_masks", "yelmo_init_topo" => "nml_init_topo", "yelmo_data" => "nml_data")

function yelmo_defaults_text()
    io = IOBuffer()
    println(io, "&yelmo")
    println(io, "    domain = \"None\"\n    grid_name = \"None\"\n    phys_const = \"Earth\"\n    restart = \"None\"\n    dt_min = 0.1")
    for g in YELMO_GROUPS
        println(io, "    ", YELMO_NML[g], " = \"", g, "\"")
    end
    println(io, "/")
    for g in YELMO_GROUPS
        println(io, "&", g, "\n    a = 1.0\n    b = \"x\"\n/")
    end
    return String(take!(io))
end

const YELMOX_DEFAULTS = """
&forcing
    reference = "pd"
    coefs     = ""
    terms     = ""
    var_defs  = "input/clim_db.nml" "input/forcing_db.nml"
/
&coef
    method      = "const"
    series_file = ""
    f_min       = 0.0
    units       = "1"
/
&term
    coef    = ""
    pattern = "snapshots"
    snap0   = ""
    snap1   = ""
    fields  = ""
/
&snapshots
    grid_global = "GLOBAL-1DEG"
    fields      = "tas" "pr" "to" "so" "zs"
    lapse       = 0.0080 0.0065
    f_p         = 0.05
/
&snapshot
    exp        = ""
    monthly    = True
    clim_sigma = 0.0
    ocn_sigma  = 0.0
    tas_season = ""
    tas        = ""
    pr         = ""
    zs         = ""
    to         = ""
    so         = ""
/
&taper
    active = False
    weight = ""
/
&bias
    active = False
    dtas   = ""
/
&ice_bnd
    active = False
    method = "file"
/
&orbit
    method = "fixed_time"
    time   = 0.0
    ecc    = 0.016702
    obl    = 23.4393
    per    = 102.918
    fldr   = "input"
/
"""

const CHION_DEFAULTS = """
&chion
    model     = "bessi"
    nml_bessi = "bessi"
    nml_pdd   = "pdd"
    nml_itm   = "itm"
    nml_const = "chion_const"
    restart   = "None"
/
&chion_const
    alpha_ice = 0.4
/
&bessi
    layers = 5
/
&pdd
    ddf_snow = 3.0
/
&itm
    alb_ice = 0.4
    itm_c   = -45.0
/
"""

const CLIM_DB = """
&tas_pd
    filename = "fesmdata/{domain}/{grid_name}/{grid_name}_TAS.nc"
    name     = "t2m"
    units    = "K" "K"
    scaling  = 1.0 0.0
    time     = 0 0 0 0 1
/
&tas_lgm
    filename = "fesmdata/Global/{grid_global}/lgm.nc"
    name     = "tas"
/
"""

const FORCING_DB = """
&coef_dT
    method = "ramp-slope"
    units  = "K"
/
&term_dT
    coef    = "dT"
    pattern = "uniform"
    fields  = "tas" "pr"
/
"""

const PHYS_CONST = """
&Earth
    sec_year = 31556926.0
    rho_ice  = 910.0
/
"""

# the par file; `edits` replace text (pairs old => new) to make variants
function par_text(; edits=Pair{String,String}[])
    s = """
&ctrl
    run_step        = "ctrl"
    calendar        = False
    calendar_ref    = 1950.0
    tstep_method    = "const"
    tstep_const     = 0.0
    time_init       = 0.0
    time_end        = 15e3
    dtt             = 10.0
    restart         = "None"
    write_lossless  = False
/
&comps
    with_ice_sheet = True
    with_isostasy  = True
    with_climate   = True
    with_surface   = True
    with_shelf     = True
    climate        = "snapshots"
    surface_method = "chion"
    dt_clim        = 10.0
    dt_isos        = 20.0
/
&sim
    init_method       = "equil"
    init_equil_time   = 10.0
    init_marine_H     = False
    init_kill_shelves = False
    init_time_thrm    = 0.0
    relax             = False
    opt               = False
    use_negis         = False
    scale_glacial_smb = False
    lim_pd_ice        = False
/
&output
    write_yelmo   = True
    write_isos    = False
    write_shelf   = True
    write_surface = True
    write_clim    = True
    write_htopo   = True
    write_cmip    = False
    dt_cmip       = 1.0
/
&domain
    name         = "Testland"
    phys_const   = "Earth"
    grid_hub     = "TST-16KM"
    grid_ice     = "TST-32KM"
    grid_isos    = "TST-32KM"
    grid_clim    = "TST-32KM"
    grid_surface = "TST-32KM"
    grid_shelf   = "TST-16KM"
    basins       = "Zwally2012"
/
&topo
    path = "fesmdata/{domain}/{grid_name}/{grid_name}_TOPO.nc"
    vars = "z_bed" "z_srf" "H_ice"
/
&regions
    path_regions = "fesmdata/{domain}/{grid_name}/{grid_name}_REGIONS.nc"
/
&tm_1D
    method = "const"
    dt     = 1.0
    file   = "input/timeout.txt"
    times  = 0, 1, 2
/
&tm_2D
    method = "const"
    dt     = 200.0
    file   = "input/timeout.txt"
    times  = 0
/
&tm_2Dsm
    method = "const"
    dt     = 50.0
    file   = "input/timeout.txt"
    times  = 0
/
&tm_rst
    method = "const"
    dt     = 20e3
    file   = "input/timeout.txt"
    times  = 0
/
&forcing
    reference = "pd"
    coefs     = "alpha" "dT"
    terms     = "dT"
/
&coef_alpha
    method      = "series"
    series_file = "input/alpha.dat"
/
&snapshot_pd
    tas = "tas_pd"
/
&orbit
    method = "fixed_time"
    time   = 0.0
/
&relax
    topo_rel = 4
    tau1     = 100.0
    tau2     = 2000.0
    time1    = 2e3
    time2    = 6e3
    m        = 2.0
/
&opt
    opt_cf       = "L21"
    cf_time_init = 0.0
    cf_time_end  = 15e3
    cf_init      = 1e-2
    tau_c        = 500.0
    H0           = 100.0
    sigma_err    = 50e3
    sigma_vel    = 100.0
    fill_method  = "cf_min"
    opt_tf       = "none"
    tf_time_init = 0.0
    tf_time_end  = 15e3
    H_grnd_lim   = 500.0
    tf_sigma     = 50e3
    tau_m        = 100.0
    m_temp       = 10.0
    tf_min       = -1.0
    tf_max       = 1.0
    tf_basins    = -1
    basin_fill   = False
/
&isos
    method  = 2
    restart = "none"
/
&barysealevel
    method  = "fastiso"
    restart = "None"
/
&sed
    method   = 0
    obs_path = "fesmdata/{domain}/{grid_name}/{grid_name}_Sediment.nc"
    obs_name = "z_sed"
/
&ghf
    use_obs      = True
    obs_path     = "fesmdata/{domain}/{grid_name}/{grid_name}_GHF.nc"
    obs_name     = "ghf"
    obs_err_name = "ghf_sd"
    f_stdev      = 0.0
    ghf_const    = 50.0
/
&yelmo
    phys_const = "Earth"
    restart    = "None"
    dt_min     = 0.01
/
&ytopo
    a = 2.0
/
&snapshots
    f_p = 0.06
/
&smbpal
    abl_method = "itm"
/
&surface_chion
    sf_a = 0.273
/
&chion
    model   = "itm"
    nml_itm = "itm_chion"
/
&itm_chion
    alb_ice = 0.31
/
&marine_shelf
    bmb_method    = "lin"
    restart       = "none"
    basins        = "Zwally2012"
    corr_method   = "tf"
    corr_names    = "b1" "b2"
    corr_values   = 0.3 -1.5
    corr_b1       = "Zwally2012:1"
    corr_b2       = "Zwally2012:2"
    tf_correction = False
    tf_path       = "ice_data/{domain}/{grid_name}/OCEAN.nc"
    tf_name       = "dT_nl"
/
"""
    for (a, b) in edits
        occursin(a, s) || error("mini par: no text $(repr(a))")
        s = replace(s, a => b; count=1)
    end
    return s
end

"Write the miniature yelmox checkout under `root` with the par file `yelmox/yelmox_Test.nml`; returns root."
function mini_yelmox(root::AbstractString; edits=Pair{String,String}[], par="yelmox_Test.nml")
    mkpath(joinpath(root, ".runme")); mkpath(joinpath(root, "input")); mkpath(joinpath(root, "yelmox"))
    write(joinpath(root, ".runme", "info.json"),
          """{"par_defaults": ["input/yelmo_defaults.nml", "input/yelmox_defaults.nml"], "par_library": ["input/forcing_db.nml"], "grp_aliases": {"bsl": "barysealevel", "ghfx": "ghf"}}""")
    write(joinpath(root, "input", "yelmo_defaults.nml"), yelmo_defaults_text())
    write(joinpath(root, "input", "yelmox_defaults.nml"), YELMOX_DEFAULTS)
    write(joinpath(root, "input", "chion_defaults.nml"), CHION_DEFAULTS)
    write(joinpath(root, "input", "clim_db.nml"), CLIM_DB)
    write(joinpath(root, "input", "forcing_db.nml"), FORCING_DB)
    write(joinpath(root, "input", "yelmo_phys_const.nml"), PHYS_CONST)
    write(joinpath(root, "input", "alpha.dat"), "0 0\n1 1\n")
    write(joinpath(root, "input", "timeout.txt"), "0\n")
    write(joinpath(root, "yelmox", par), par_text(; edits=edits))
    return root
end
