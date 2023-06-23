import sys
import numpy as np
import matplotlib.pyplot as plt
import lmfit as lm

# see: https://lmfit.github.io/lmfit-py/faq.html#how-can-i-fit-multiple-data-sets

def simulate_data(x):
    nsets = 2
    y = np.zeros( (nsets, len(x)), dtype=float )

    slope = 5.0 # shared slope
    yint0 = 2.0
    yint1 = -3.0
    
    y[0,:] = linear(x, yint0, slope) + np.random.normal(size=x.size, scale=1.0)
    y[1,:] = linear(x, yint1, slope) + np.random.normal(size=x.size, scale=1.0)
    
    return y

def linear(x, yint, slope):
    return yint + slope*x

def linear_model(params, ii, x):
    slope = params['slope']
    yint = params[f'yint_{ii}']
    return linear(x, yint, slope)

def objective(params, x, data):
    ndata, _ = data.shape
    resid = 0.0*data[:]
    
    models = np.zeros_like(data)
    for ii in range(ndata):
        models[ii,:] = params['slope']*x + params[f'yint_{ii}']

    # make residuals per data set
    for ii in range(ndata):
        resid[ii,:] = data[ii,:] - models[ii,:]

    return resid.flatten()

        
if __name__ == "__main__":
    
    x = np.linspace(-3, 3, num=100)
    y = simulate_data(x)
    ndata = len(y)
    
    params = lm.Parameters()
    params.add(f"slope", value=0)
    for ii in range(ndata):
        params.add(f"yint_{ii}", value=0)
    params.pretty_print()


    out = lm.minimize(objective, params, args=(x,y))
    lm.report_fit(out.params)

    for ii in range(ndata):
        yfit = linear_model(out.params, ii, x)
        plt.plot(x, y[ii,:], '.', x, yfit, '-')
    plt.show()
