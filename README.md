# saq_analysis
analysis of data from Q-Pix SAQ

So far:
1. `make_templates_root.py`  
   Simulate reset data (aka "templates"). Run many instances in parallel (on noether).  
   Produces several `.root` files containing parameter values and simulated resets (per area).  
2. `templates_merge.py`  
   Produce a single Pandas DataFrame containing *all* templates. Store in a pickle file
3. `template_compute_chisq_pandas.py`  
   Compute and save the optimal scaling factor and sum of squares of residuals between data and model
4. TBD  
   For each set of nuissance parameters, find the best set of diffusion parameters for the data.

Helper functions
* `offset_simulation.py`
  A collection of useful helper functions. (`import offset_simulation as off`)
* `pressure_scan_data_to_root.py`
  Convert a collection of `.csv` files containing reset data into a single `.root` file
* `pressure_scan_data_read_root.py`
  Read in ROOT files containing pressure scan data

To do:
* could likely avoid ROOT altogether and just use Pandas (e.g. for reading/saving pressure scan data)
