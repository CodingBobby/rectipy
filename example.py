"""
rectipy module by
@author: Bob Walter

This file provides an example of usage.
To design a rectification process, this is the simplest way to do so.
"""

import rectipy as rp

# first, we create the components, our mixture will be made of
# the argument is a string of the lowercase name of the component
# for available components and their exact spelling, check the README
water = rp.Component('water')
ethanol = rp.Component('ethanol')

# then, we mix them together
# it is important in which order you pass the components,
# as the first one will be reference to concentration
# values in the following procedure
mixture = rp.Mixture(ethanol, water)

# and now initialize the rectification process
# the required arguments are:

# mixture: the created Mixture object
# xA: concentration of the swamp
# xF: concentration of the feed
# xD: concentration of the product
rectification = rp.Rectification(mixture, 0.12, 0.4, 0.78)

# with this, we get the McCabe-Thile plot of our mixture
# with operation lines and rectification steps drawn in
rectification.plot()
