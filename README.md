#Semiconductor

This is a place I used to get together a bunch of different analytical models and
tabulated data sets for semiconductor properties. Its main focus is Silicon, as
that is what I work with. 

All the models are implemented in a similar way. A class is built to allows switching
of model written by different authors through the author command. The model will 
then result in the values from the authors implementation. That is, this is module
is meant to allow replication of model.

## Example 

Here is an example of how to use this module. We will look at two different band gap narrowing models. The default model is that from Yan in 2014.
model.

```python
    from semiconductor.matterial.bandgap_narrowing import BandGapNarrowing as BGN
    import numpy as np

    # initialise the class
    BGN_class = BGN(matterial='Si')
    # define the number of dopants
    Na = 0.
    Nd = np.logspace(16, 20)
    # Set the excess carriers to zero
    min_car_den = 0
    bgn_yan = BGN_class.update_BGN(Na, Nd, min_car_den)
```

If a different band gap narrowing model is desired, 
pick from the available ones. The available ones can be found
using the available_models() function.

```python
    print BGN_class.available_models()
```

For the band gap narrowing function it returns.

```python
    ['DelAlamo1985', 'Cuevas1996', 'Yan2013bz', 'Yan2014bz', 'Schenk1988fer', 'Schenk1988_reparamitisation_Yan2013', 'Yan2013fer', 'Yan2014fer']
```

Changing to a model by a different author is done using the author input in either
the  initalisation of the class, or through the "update" function. Lets 
choose Schenk's from 1988 and set it through the "update" function. All classes
have a similar update function. If we look at the models inputs, we see it also needs
an input for temperature. This is just passed to the update function, which passes
it to the appropriate places.  
```python
    temp = 300
    bgn_sch = BGN_class.update_BGN(Na, Nd, min_car_den, temp=300, author='Schenk1988fer')
```

Finally we can plot, and compare the differences in the models.

``` python
    plt.plot(Nd, bgn_yan, label = 'Yan')
    plt.plot(Nd, bgn_sch, label = 'Schenk')
    plt.legend(loc=0, title='Author')
    plt.xlabel('Doping')
    plt.ylabel('Band Gap Narrowing (eV)')
    plt.semilogx()
    plt.show()
```

![Comparison of Yan's and Schenk's band gap narrowing models](comparison.png)