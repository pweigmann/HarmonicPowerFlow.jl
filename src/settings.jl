# Structs
struct Settings
    coupled::Bool  # choose coupled or uncoupled Norton Equivalent model
    harmonics::Array  # only uneven harmonics, must contain "1", example: [1,3,5,7]
    base_frequency::Int  # in Hz, default: 50
    base_voltage::Number  # in V, default: 230
    base_power::Number  # in kVA, default: 1000

    # algorithm parameters
    thresh_f
    thresh_h
    max_iter_f::Int
    max_iter_h::Int

    # voltage start values [pu]
    v1 
    ϕ1 
    vh 
    ϕh

    # calculated parameters
    K::Int  # number of harmonics (without fundamental)
    base_current::Number
    base_admittance::Number
end


# Functions
function init_settings(coupled, harmonics; base_frequency=50, base_voltage=230,
    base_power=1000, thresh_f = 1e-6, max_iter_f = 30, thresh_h=1e-4, max_iter_h=50,
    v1=1, ϕ1 = 0, vh=0.1, ϕh=0)

    K = length(harmonics) - 1

    # derived pu base values
    base_current = 1000*base_power/base_voltage
    base_admittance = base_current/base_voltage
    
    Settings(
        coupled, harmonics, base_frequency, base_voltage, base_power, 
        thresh_f, thresh_h, max_iter_f, max_iter_h, 
        v1, ϕ1, vh, ϕh, 
        K, base_current, base_admittance)
end