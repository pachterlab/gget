> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python. The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
## gget setup 🔧

Function to install/download third-party dependencies for a specified gget module.

**Positional argument**  
`module`  
gget module for which dependencies should be installed (currently only "alphafold").  

### Example
```bash
gget setup alphafold
```
```python
# Python
gget.setup("alphafold")
```
&rarr; Installs all (modified) third-party dependencies and downloads model parameters (~4GB) required to run `gget alphafold`. 
