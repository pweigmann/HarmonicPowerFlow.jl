var documenterSearchIndex = {"docs":
[{"location":"functions/#Detailed-API","page":"Functions","title":"Detailed API","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"HarmonicPowerFlow","category":"page"},{"location":"functions/#HarmonicPowerFlow","page":"Functions","title":"HarmonicPowerFlow","text":"A simple and fast module to solve the harmonic power flow problem for given distribution grids.\n\n\n\n\n\n","category":"module"},{"location":"functions/","page":"Functions","title":"Functions","text":"Modules = [HarmonicPowerFlow]\nOrder   = [:constant, :type, :function, :macro]","category":"page"},{"location":"functions/#HarmonicPowerFlow.PowerGrid","page":"Functions","title":"HarmonicPowerFlow.PowerGrid","text":"Struct containing information about nodes and lines.\n\n\n\n\n\n","category":"type"},{"location":"functions/#HarmonicPowerFlow.Settings","page":"Functions","title":"HarmonicPowerFlow.Settings","text":"Struct containing simulation settings.\n\n\n\n\n\n","category":"type"},{"location":"functions/#HarmonicPowerFlow.THD-Tuple{Any}","page":"Functions","title":"HarmonicPowerFlow.THD","text":"THD(u)\n\nCalculate the Total Harmonic Distortion at all nodes.\n\nTHD is a measure for the amount of voltage distortion present at a node.  This function calculates two alternative definitions \n\nTHD_F: relative to fundamental, more commonly used and \nTHD_R: RMS, relative to all frequencies\nTHD_40: relative to fundamental, only until h = 40 (or 39) as used by   standard EN 61000.\n\n\n\n\n\n","category":"method"},{"location":"functions/#HarmonicPowerFlow.admittance_matrices-Tuple{Any, Any}","page":"Functions","title":"HarmonicPowerFlow.admittance_matrices","text":"admittance_matrices(nodes, lines, harmonics)\n\nBuild the nodal admittance matrices (admittance laplacian) for all harmonics. \n\nAdmittance scales with frequency: Xh = Xf * h. Return a dictionary of DataFrames.\n\n\n\n\n\n","category":"method"},{"location":"functions/#HarmonicPowerFlow.barplot_THD-Tuple{Any, Any}","page":"Functions","title":"HarmonicPowerFlow.barplot_THD","text":"barplot_THD(u, nodeID; h_min=1, h_max=0, yticks=0:0.1:1)\n\nBarplot of harmonic voltage magnitudes at one node.\n\n\n\n\n\n","category":"method"},{"location":"functions/#HarmonicPowerFlow.build_harmonic_jacobian-NTuple{5, Any}","page":"Functions","title":"HarmonicPowerFlow.build_harmonic_jacobian","text":"build_harmonic_jacobian(net, settings, u, LY, NE)\n\nConstruct harmonic Jacobian J containing all derivatives of power and current mismatches wrt voltages from harmonic state vector.\n\n\n\n\n\n","category":"method"},{"location":"functions/#HarmonicPowerFlow.create_nodes_manually-Tuple{}","page":"Functions","title":"HarmonicPowerFlow.create_nodes_manually","text":"createnodesmanually()\n\nManually create a nodes DataFrame. \n\nNote that linear nodes are added first, then nonlinear ones.  Within the linear ones, first add PV, then PQ nodes. First node contains the slack.\n\n\n\n\n\n","category":"method"},{"location":"functions/#HarmonicPowerFlow.current_balance-NTuple{5, Any}","page":"Functions","title":"HarmonicPowerFlow.current_balance","text":"current_balance(net, settings, u, LY, NE)\n\nEvaluate current balance at all nodes and frequencies, except linear nodes at fundamental frequency.\n\n\n\n\n\n","category":"method"},{"location":"functions/#HarmonicPowerFlow.current_injections-NTuple{5, Any}","page":"Functions","title":"HarmonicPowerFlow.current_injections","text":"current_injections(nodes, nodeID, u, NE, harmonics)\n\nCalculate the harmonic current injections at one node.\n\nUses coupled or uncoupled method based on dimensions of NE.\n\n\n\n\n\n","category":"method"},{"location":"functions/#HarmonicPowerFlow.fund_jacobian-Tuple{Any, Any, Any}","page":"Functions","title":"HarmonicPowerFlow.fund_jacobian","text":"fund_jacobian(net, u, LY)\n\nConstruct fundamental Jacobian, containing powers derived wrt voltages.\n\n\n\n\n\n","category":"method"},{"location":"functions/#HarmonicPowerFlow.fund_mismatch-Tuple{Any, Any, Any}","page":"Functions","title":"HarmonicPowerFlow.fund_mismatch","text":"fund_mismatch(net, u, LY)\n\nEvaluate active and reactive power mismatch at fundamental frequency.\n\n\n\n\n\n","category":"method"},{"location":"functions/#HarmonicPowerFlow.fund_state_vec-Tuple{Any, Any}","page":"Functions","title":"HarmonicPowerFlow.fund_state_vec","text":"fund_state_vec(net, u)\n\nTake the voltage DataFrame u and return a vector of voltage phases and magnitudes.\n\n\n\n\n\n","category":"method"},{"location":"functions/#HarmonicPowerFlow.harmonic_mismatch-NTuple{5, Any}","page":"Functions","title":"HarmonicPowerFlow.harmonic_mismatch","text":"harmonic_mismatch(net, settings, u, LY, NE)\n\nCombine current and power mismatch to construct harmonic mismatch vector f.\n\n\n\n\n\n","category":"method"},{"location":"functions/#HarmonicPowerFlow.harmonic_state_vec-Tuple{Any, Any, Any}","page":"Functions","title":"HarmonicPowerFlow.harmonic_state_vec","text":"harmonic_state_vec(net, u, harmonics)\n\nConstruct harmonic state vector x from fundamental and harmonic voltages.\n\n\n\n\n\n","category":"method"},{"location":"functions/#HarmonicPowerFlow.heatmap_NE-Tuple{Any, Any}","page":"Functions","title":"HarmonicPowerFlow.heatmap_NE","text":"heatmap_NE(nodes, settings; save=png)\n\nExport heatmap indicating magnitude of Norton equivalent parameters tp .png or .pdf.\n\n\n\n\n\n","category":"method"},{"location":"functions/#HarmonicPowerFlow.hpf-Tuple{Any, Any}","page":"Functions","title":"HarmonicPowerFlow.hpf","text":"hpf(net, settings)\n\nSolve the harmonic power flow problem for a given PowerGrid and Settings.\n\n\n\n\n\n","category":"method"},{"location":"functions/#HarmonicPowerFlow.import_Norton_Equivalents","page":"Functions","title":"HarmonicPowerFlow.import_Norton_Equivalents","text":"import_Norton_Equivalents(nodes, settings, folder_path=\"devices\\\")\n\nImport NE parameters for all nonlinear devices in nodes from folder devices.\n\n\n\n\n\n","category":"function"},{"location":"functions/#HarmonicPowerFlow.init_power_grid-Tuple{DataFrames.DataFrame, DataFrames.DataFrame, Any}","page":"Functions","title":"HarmonicPowerFlow.init_power_grid","text":"init_power_grid(nodes, lines, settings)\n\nCreate a PowerGrid struct from nodes and lines DataFrames and converts all quantities to the p.u. system.\n\nIndexing (c, m, n) could be done better by being based on number of nodes of certain type. It could make sense to also add LY to PowerGrid struct.\n\n\n\n\n\n","category":"method"},{"location":"functions/#HarmonicPowerFlow.init_settings-Tuple{Any, Any}","page":"Functions","title":"HarmonicPowerFlow.init_settings","text":"init_settings(coupled, harmonics; ...)\n\nInitialize settings struct. \n\nDefault values:\n\nbase_frequency=50\nbase_voltage=230 \nbase_power=1000 \nthresh_f = 1e-6 \nmaxiterf = 30 \nthresh_h=1e-4\nmaxiterh=50\nv_f=1\nϕ_f = 0 \nv_h=0.1 \nϕ_h=0\n\n\n\n\n\n","category":"method"},{"location":"functions/#HarmonicPowerFlow.init_voltages-Tuple{Any, Any}","page":"Functions","title":"HarmonicPowerFlow.init_voltages","text":"init_voltages(nodes, settings)\n\nInitialize voltages for all nodes and harmonics using guesses from on Settings.\n\n\n\n\n\n","category":"method"},{"location":"functions/#HarmonicPowerFlow.limits-Tuple{Any, Any}","page":"Functions","title":"HarmonicPowerFlow.limits","text":"limits(u, nodeID)\n\nCheck if harmonic limits according to EN 50160 and EN 61000 are met at node.\n\n\n\n\n\n","category":"method"},{"location":"functions/#HarmonicPowerFlow.pf","page":"Functions","title":"HarmonicPowerFlow.pf","text":"pf(net, settings, LY, plt_convergence = false)\n\nPerform fundamental power flow, iterate until mismatch is smaller than threshold.\n\n\n\n\n\n","category":"function"},{"location":"functions/#HarmonicPowerFlow.update_fund_voltages!-Tuple{Any, Any, Any}","page":"Functions","title":"HarmonicPowerFlow.update_fund_voltages!","text":"update_fund_voltages!(net, u, x)\n\nWrite new voltages from x to u.\n\n\n\n\n\n","category":"method"},{"location":"functions/#HarmonicPowerFlow.update_harmonic_voltages!-NTuple{4, Any}","page":"Functions","title":"HarmonicPowerFlow.update_harmonic_voltages!","text":"update_fund_voltages!(net, u, x, settings)\n\nWrite new voltages from x to u for all harmonics.\n\n\n\n\n\n","category":"method"},{"location":"functions/#HarmonicPowerFlow.update_state_vec!-Tuple{Any, Any, Any}","page":"Functions","title":"HarmonicPowerFlow.update_state_vec!","text":"update_state_vec!(J, x, f)\n\nPerform Newton-Raphson iteration and update the fundamental state vector.\n\n\n\n\n\n","category":"method"},{"location":"conventions/#Conventions","page":"Conventions","title":"Conventions","text":"","category":"section"},{"location":"conventions/","page":"Conventions","title":"Conventions","text":"In this project \"node\" and \"bus\" are used to describe the same object and fully interchangeable terms.","category":"page"},{"location":"conventions/","page":"Conventions","title":"Conventions","text":"A linear/nonlinear bus is a network bus with a linear/nonlinear device connected to it.","category":"page"},{"location":"conventions/#Abbreviations","page":"Conventions","title":"Abbreviations","text":"","category":"section"},{"location":"conventions/","page":"Conventions","title":"Conventions","text":"Abbr. Meaning\nPF Power Flow\nHPF Harmonic Power Flow\nHCNE Harmonic Coupled Norton Equivalent (method)\nNL Nonlinear\nNE Norton Equivalent(s)\nSMPS Switched-mode power supply\nNR Newton-Raphson (algorithm)\nTHD Total Harmonic Distortion\nRMS Root Mean Square","category":"page"},{"location":"conventions/#Nomenclatur","page":"Conventions","title":"Nomenclatur","text":"","category":"section"},{"location":"conventions/","page":"Conventions","title":"Conventions","text":"Symbol (Code) Symbol (Math) Object/Parameter\nLY_h Y^h (Nodal) Admittance matrix at harmonic h\nu_h u^h complex voltage at harmonic h","category":"page"},{"location":"conventions/#Indexing","page":"Conventions","title":"Indexing","text":"","category":"section"},{"location":"conventions/","page":"Conventions","title":"Conventions","text":"The indexing can look a bit confusing at first glance. It was chosen this way to make it comparable to some of the literature sources and for easier implementation when using the indexes for slicing DataFrames.","category":"page"},{"location":"conventions/#Buses","page":"Conventions","title":"Buses","text":"","category":"section"},{"location":"conventions/","page":"Conventions","title":"Conventions","text":"Variable Meaning\ni bus index, i = 1 for slack bus\nn total number of buses\nc-1 number of PV buses (excluding slack), i = 2  c\nd-1 number of PQ buses, i = c+1  m-1\nm-1 total number of linear buses, i = 1  m-1\nn-m+1 number of nonlinear buses, i = m  n","category":"page"},{"location":"conventions/","page":"Conventions","title":"Conventions","text":"This implies:","category":"page"},{"location":"conventions/","page":"Conventions","title":"Conventions","text":"total buses (n) = slack bus (1) + PV buses (c-1) + PQ buses (d-1) + nonlinear buses (n-m+1)","category":"page"},{"location":"conventions/","page":"Conventions","title":"Conventions","text":"linear buses (m-1) = slack bus (1) + PV buses (c-1) + PQ buses (d-1)","category":"page"},{"location":"conventions/","page":"Conventions","title":"Conventions","text":"or","category":"page"},{"location":"conventions/","page":"Conventions","title":"Conventions","text":"n   = 1 + (c-1) + (d-1) + (n-m+1)","category":"page"},{"location":"conventions/","page":"Conventions","title":"Conventions","text":"(m-1) = 1 + (c-1) + (d-1) Rightarrow m = c + d","category":"page"},{"location":"conventions/#Harmonics","page":"Conventions","title":"Harmonics","text":"","category":"section"},{"location":"conventions/","page":"Conventions","title":"Conventions","text":"Variable Meaning\nh harmonic index, h = 1 (or sometimes h = f) for fundamental frequency\nK number of considered harmonics (excluding fundamental)\nK_1 number of considered harmonics (including fundamental)\nL highest harmonic considered","category":"page"},{"location":"conventions/#Other","page":"Conventions","title":"Other","text":"","category":"section"},{"location":"conventions/","page":"Conventions","title":"Conventions","text":"Variable Meaning\nY_N I_N Norton Equivalent parameter (admittance and current)\nI_inj nonlinear device current injections","category":"page"},{"location":"conventions/#Putting-indices-together","page":"Conventions","title":"Putting indices together","text":"","category":"section"},{"location":"conventions/","page":"Conventions","title":"Conventions","text":"Generally harmonics will be upper indices, bus numbers and other indices lower ones.","category":"page"},{"location":"conventions/","page":"Conventions","title":"Conventions","text":"Indices k and p are used as to denote bus and harmonic respectively when used to form derivatives. So fracpartial I_inj i^hpartial theta_k^p is the derivative of the current injection at bus i and harmonic h with respect to the voltage angle theta at bus k and harmonic p.","category":"page"},{"location":"conventions/","page":"Conventions","title":"Conventions","text":"When it comes to admittances, this means Y_ij^h is the network admittance at harmonic h  between bus i and j, while Y_N i^hq is the Norton parameter which describes the coupling of harmonics h and q at bus i.","category":"page"},{"location":"conventions/#Per-Unit-System","page":"Conventions","title":"Per-Unit System","text":"","category":"section"},{"location":"conventions/","page":"Conventions","title":"Conventions","text":"'base_voltage' is also grid voltage.","category":"page"},{"location":"import/#Import","page":"Import","title":"Import","text":"","category":"section"},{"location":"import/","page":"Import","title":"Import","text":"PowerGrids can only be imported from .csv files at the moment. ","category":"page"},{"location":"import/","page":"Import","title":"Import","text":"The files should be named <grid_name>_buses.csv and <grid_name>_lines.csv and contain the columns","category":"page"},{"location":"import/","page":"Import","title":"Import","text":"ID ; type ; component ; S ; X_shunt ; P ; Q","category":"page"},{"location":"import/","page":"Import","title":"Import","text":"and","category":"page"},{"location":"import/","page":"Import","title":"Import","text":"ID ; fromID ; toID ; R ; X","category":"page"},{"location":"import/","page":"Import","title":"Import","text":"respectively. Use SI Units and enter nothing if you want to leave the field component empty.","category":"page"},{"location":"import/","page":"Import","title":"Import","text":"See this example for reference.","category":"page"},{"location":"#HarmonicPowerFlow.jl","page":"General","title":"HarmonicPowerFlow.jl","text":"","category":"section"},{"location":"","page":"General","title":"General","text":"(Image: codecov) (Image: CI) (Image: Code on Github.)","category":"page"},{"location":"#Module-index","page":"General","title":"Module index","text":"","category":"section"},{"location":"","page":"General","title":"General","text":"Modules = [HarmonicPowerFlow]\nOrder   = [:constant, :type, :function, :macro]","category":"page"}]
}
