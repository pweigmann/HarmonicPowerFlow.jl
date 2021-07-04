# Conventions

In this project "node" and "bus" are used to describe the same object and fully interchangeable terms.

A linear/nonlinear bus is a network bus with a linear/nonlinear device connected to it.



## Abbreviations
Abbr. | Meaning
--- | :---
PF | Power Flow
HPF | Harmonic Power Flow
HCNE | Harmonic Coupled Norton Equivalent (method)
NL | Nonlinear
NE | Norton Equivalent(s)
SMPS | Switched-mode power supply
NR | Newton-Raphson (algorithm)
THD | Total Harmonic Distortion
RMS | Root Mean Square

## Nomenclatur

Symbol (Code) | Symbol (Math) | Object/Parameter
:---: | :---: | :---
LY_h | $Y^h$ | (Nodal) Admittance matrix at harmonic h
u_h | $u^h$ | complex voltage at harmonic h


## Indexing
The indexing can look a bit confusing at first glance. It was chosen this way to make it comparable to some of the literature sources and for easier implementation when using the indexes for slicing DataFrames.

### Buses

Variable | Meaning
--- | :---
$i$ | bus index, $i = 1$ for slack bus
$n$ | total number of buses
$c-1$ | number of PV buses (excluding slack), $i = 2, ..., c$
$d-1$ | number of PQ buses, $i = c+1, ..., m-1$
$m-1$ | total number of linear buses, $i = 1, ..., m-1$
$n-m+1$ | number of nonlinear buses, $i = m, ..., n$

This implies:

total buses (n) = slack bus (1) + PV buses (c-1) + PQ buses (d-1) + nonlinear buses (n-m+1)

linear buses (m-1) = slack bus (1) + PV buses (c-1) + PQ buses (d-1)

or

$n   = 1 + (c-1) + (d-1) + (n-m+1)$

$(m-1) = 1 + (c-1) + (d-1) \Rightarrow m = c + d$

### Harmonics

Variable | Meaning
--- | :---
$h$ | harmonic index, $h = 1$ (or sometimes $h = f$) for fundamental frequency
$K$ | number of considered harmonics (excluding fundamental)
$K_1$ | number of considered harmonics (including fundamental)
$L$ | highest harmonic considered

### Other

Variable | Meaning
--- | :---
$Y_N, I_N$ | Norton Equivalent parameter (admittance and current)
$I_{inj}$ | nonlinear device current injections


### Putting indices together

Generally harmonics will be upper indices, bus numbers and other indices lower ones.

Indices $k$ and $p$ are used as to denote bus and harmonic respectively when used to form derivatives. So $\frac{\partial I_{inj, i}^h}{\partial \theta_k^p}$ is the derivative of the current injection at bus $i$ and harmonic $h$ with respect to the voltage angle $\theta$ at bus $k$ and harmonic $p$.

When it comes to admittances, this means $Y_{ij}^h$ is the network admittance at harmonic $h$  between bus $i$ and $j$, while $Y_{N, i}^{hq}$ is the Norton parameter which describes the coupling of harmonics $h$ and $q$ at bus $i$.

## Per-Unit System

'base_voltage' is also grid voltage.