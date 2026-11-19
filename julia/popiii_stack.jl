# The one supported Pop III configuration (docs/POPIII_STACK.md, docs/POPIII_SWITCH_AUDIT2.md): the settings of the main run
# z300mm (16 Nov 2026) as the defaults of popiii_m2c.jl. Included first by popiii_m2c.jl: every key that the environment does not set
# is set here, before any module reads ENV, so a run needs only its own settings (mesh, cloud, output) and the start. An explicit
# setting still overrides (tests, other meshes). Not used by any other script (the S4/S5, quasi-star and harness scripts keep the
# code defaults).
#
# Groups: backend; IC; material (network everywhere outside the dense interior: G2 SCVH95+x_e table, HM88, network continuum
# opacity, three-level H); chemistry (coupled implicit, reduced CK network, HD, He); H₂ lines (hybrid ALI/EP, at the iterate);
# star light at the iterate; radiation (Eddington factors at the iterate, sub-zone front closure, flux tie in the thin outer
# zones); energy (total-energy row); hydro (zone-width viscosity); mesh (log mesh from a regular centre, GRIDEQ moving mesh);
# step control (TR-BDF2 error estimate); run control.
const POPIII_STACK = [
    # backend
    "TITAN_NOLIB" => "1", "TITAN_MGR" => "600",
    # initial condition: ABN02 z = 19 cloud, the reduced-x_e variant (7137 M☉), T_CMB at z = 19
    "IC" => "abn02red", "MCLOUD" => "7137", "TCMB" => "54.5", "RAD" => "1",
    # material
    "CHEM" => "netmodernli", "HOT" => "1", "COLDMLT" => "1e4", "ALPHAMLT" => "1.5",
    "SCVHXE" => "1", "HM88" => "1", "G2OPAC" => "2", "NLOHBFDEEP" => "1", "TL3DEEP" => "1",
    "PSWITCH" => "1e30", "PSWTAU" => "0", "PSWSTAY" => "0", "PSWTEX" => "5000", "PSWNH" => "1e19", "PSWHYS" => "0.3",
    "PSWTMIN" => "4500", "PSWRAMP" => "10",
    "NLOPAC" => "1", "NLOPACJ" => "1",
    "LOCALJ" => "1", "CLOSJ" => "1",     # CLOSJ: band Jacobian carries the Eddington/tie closure dependence (TitanDriver.CSENS)
    # chemistry
    "COUPLED" => "1",   # (read by TitanFortran: meqn ≥ 9 for the two species rows; popiii_m2c always couples)
    "CHEMCK" => "1", "CHEMCKHSS" => "100", "CHEMSUB" => "1.3,3", "CHEMFBNET" => "16",
    "HDCHEM" => "1", "HECHEM" => "1", "CHEMTHREADS" => "1", "MATTHREADS" => "8",
    # H₂ lines
    "NLTE" => "hyb", "LINEINC" => "1", "LINECOLL" => "1", "LINESRC" => "zone", "LINEQUAD" => "p2",
    "LINECORE" => "8", "LINEFAST" => "1", "EPTAB" => "1", "ALIWMIN" => "1e-3", "ALISHRINK" => "0.5", "ALIGROW" => "0.1",
    # star light
    "STARLIGHT" => "1", "SL3" => "1", "SLBLEND" => "1", "SLEB" => "1", "SLWFLUX" => "1", "SLSHIELD" => "sobolev", "LYABB" => "1",
    "SLIMP" => "2",
    # radiation
    "VEFIMP" => "1", "EDDGMU0" => "1", "EDDSCAT" => "1", "SUBZONE" => "1,2", "FLUXTIE" => "0.1,1", "ERTAUMIN" => "3000",
    # energy
    "ETOTAL" => "1", "ETASW" => "0.6", "ETASWTH" => "1e-3", "ETASWTIE" => "0.5", "EXCAP" => "1e8",
    # hydro
    "AVZONE" => "2", "AVTENSION" => "0",
    # mesh: 290 zones, Δln r 0.0794 from a regular centre (one central zone to 0.1 R☉); GRIDEQ moving mesh (iterate-only monitor)
    "ENVEUL" => "0,0.0793995,0.0793995,1e30", "ENVCORE" => "0.1,0.0793995,1e30", "ENVLIN" => "1", "ENVCONS" => "1",
    "GRIDEQ" => "1", "GEQRAMP" => "1e-3,300,10,20", "GEQC" => "100,1.5,30,0.5,0.5,0.1,1e4,0.3,0.3", "GEQW" => "1,1,1,1,1,0,0,3",
    # step control: TR-BDF2 embedded error estimate (TITAN's rule below is then diagnostic only)
    "DTERR" => "0.1", "CONV" => "1e-3", "NITER" => "25", "NTRY" => "60", "CDYN" => "3e5", "STOL" => "0.3", "FAILTFAC" => "0.5",
    # run control
    "TSTOP" => "2e8", "NSTOP" => "1e27", "DTSTOP" => "1e-9", "CKEVERY" => "100", "SNAPEVERY" => "50", "PRINTEVERY" => "100",
    "MSG" => "1", "FLUSHOUT" => "1",
]
for (k, v) in POPIII_STACK
    haskey(ENV, k) || (ENV[k] = v)
end
