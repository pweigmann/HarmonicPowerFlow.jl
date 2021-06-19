# Import

PowerGrids can only be imported from `.csv` files at the moment. 

The files should be named `<grid_name>_buses.csv` and `<grid_name>_lines.csv` and contain the columns
```
ID ; type ; component ; S ; X_shunt ; P ; Q
```
and
```
ID ; fromID ; toID ; R ; X
```
respectively. Use SI Units and enter `nothing` if you want to leave the field `component` empty.

See this [example](https://github.com/pweigmann/HarmonicPowerFlow.jl/blob/41fecd02f4d184b9d1dff49137af24e8678c13d1/examples/simple_coupled_hpf/) for reference.

