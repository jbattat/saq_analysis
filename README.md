# saq_analysis
analysis of data from Q-Pix SAQ

So far:
1. `make_templates_root.py`  
   Simulate reset data (aka "templates"). Run many instances in parallel (on noether).  
   Produces several `.root` files containing parameter values and simulated resets (per area).  
2. `templates_merge.py`  
   Produce a single Pandas DataFrame containing *all* templates. Store in a pickle file
3. `template_compute_chisq.py`  
   Compute and save the optimal scaling factor and sum of squares of residuals between data and model. Store in a pickle file.
4. `template_compute_optimal_parameter_sets.py`  
   For each set of nuissance parameters, find the set of diffusion parameters that minimizes the overall chi-squared. Store in a pickle file.
5. `template_analysis.py`  
   Plot the data and best fit models (not yet mature).

Helper functions
* `offset_simulation.py`
  A collection of useful helper functions. (`import offset_simulation as off`)
* `pressure_scan_data_to_root.py`
  Convert a collection of `.csv` files containing reset data into a single `.root` file
* `pressure_scan_data_read_root.py`
  Read in ROOT files containing pressure scan data

To do:
* could likely avoid ROOT altogether and just use Pandas (e.g. for reading/saving pressure scan data)
