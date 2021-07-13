# Structs
struct Settings
    coupled::Bool  # choose coupled or uncoupled Norton Equivalent model
    harmonics::Array  # only uneven harmonics, must contain "1", example: [1,3,5,7]
    base_frequency::Int  # in Hz, default: 50
    base_voltage::Number  # in V, default: 230
    base_power::Number  # in W, default: 1000

    # algorithm parameters
    thresh_f
    thresh_h
    max_iter_f::Int
    max_iter_h::Int

    # voltage start values [pu]
    v_f 
    ϕ_f 
    v_h 
    ϕ_h

    # calculated parameters
    K::Int  # number of harmonics (without fundamental)
    K1::Int # number of harmonics (including fundamental)
    base_current::Number
    base_admittance::Number
    base_impedance::Number
end


# Functions
function init_settings(coupled, harmonics; 
    base_frequency=50, base_voltage=230, base_power=1000, 
    thresh_f = 1e-6, max_iter_f = 30, thresh_h=1e-4, max_iter_h=50,
    v_f=1, ϕ_f = 0, v_h=0.1, ϕ_h=0)

    K = length(harmonics) - 1
    K1 = length(harmonics)

    # derived pu base values
    base_current = base_power/base_voltage
    base_admittance = base_current/base_voltage
    base_impedance = 1/base_admittance
    
    Settings(
        coupled, harmonics, base_frequency, base_voltage, base_power, 
        thresh_f, thresh_h, max_iter_f, max_iter_h, 
        v_f, ϕ_f, v_h, ϕ_h, 
        K, K1, base_current, base_admittance, base_impedance)
end