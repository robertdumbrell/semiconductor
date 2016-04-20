#!/usr/local/bin/python
# UTF-8

import matplotlib.pylab as plt
import numpy as np

"""
To do:

    Need to fix the check_doping function
"""

def change_model( Models, author=None):

    author = author or Models.get('default', 'model')

    model = Models.get(author, 'model')

    vals = dict(Models.items(author))

    del vals['model']

    for k, v in vals.iteritems():

        # checks if float or list
        try:
            vals[k] = float(v)

        except:
            try:
                vals[k] = [float(i) for i in v.split(';')]
                # print self.vals[k]
            except:
                pass

    return vals, model


class HelperFunctions():

    def change_model(self, author, Models=None):

        Models = Models or self.Models

        self.vals, self.model = change_model(Models, author)

    def plot_all_models(self, update_function, xvalues=None, **kwargs):
        '''
        cycles through all the models and plots the result
        inputs:
            update_function: str
                 the name of the specific function used to update the author
                i.e 'update_Eg'
            **kwargs:
                variables to be passed to the update function. 
        '''
        fig, ax = plt.subplots(1)
        # ax = plt.add_subplot(111)
        for model in self.available_models():
            print update_function, model
        # ax.plot(np.inf,np.inf,'k-',label = 'Auger')
            self.change_model(model)
            result = getattr(self, update_function)(**kwargs)
            if xvalues is None:
                ax.plot(result, label=model)
            else:
                ax.plot(xvalues, result, label=model)
        plt.legend(loc=0)

    def plotting_colours(self, n_colours, fig, ax, repeats=None):
        '''
        a function to help with plotting
        Lets you get a unique range of colours, 
        and have repeats of colours
        '''

        colours = []

        # Define the colours for the figures, using CMRmap to get the colours
        for i in np.linspace(.1, .7, n_colours):

            # if don't want any repeat of colours
            if repeats is None:
                colours.append(plt.cm.CMRmap(i))
            else:
                for j in range(repeats):
                    colours.append(plt.cm.CMRmap(i))

        # incase its a 2D axis
        try:
            for axes in ax.flatten():
                axes.set_color_cycle(colours)
        except:
            ax.set_color_cycle(colours)

    def available_models(self, Filter=None, Filter_value=None):
        '''
        Returns a list of all authors for a given parameter this a class of.
        The returned list can be filtered in field "Filter" by the value "Filter_value"

        inputs: (all optional)
            Filter: (str)
                The model field that is to be checked        
            Filter_value:
                The value to be checked for
        returns:
            list of authors
        '''

        author_list = self.Models.sections()
        author_list.remove('default')

        # remove modles that are not implimented
        for author in list(author_list):
            # print  dict(self.Models.items(author))['model']
            # print 'not_implimented' in dict(self.Models.items(author))['model']
            # print author
            if 'not_implimented' in dict(self.Models.items(author))['model']:
                # print author, 'removed'
                author_list.remove(author)

        # does the filtering
        if Filter is not None:
            for c, author in enumerate(author_list):
                if self.Models.get(author, Filter) not in Filter_value:
                    author_list.remove(author)

        # prints no models available
        if not author_list:
            print 'No authors for this models available'

        return author_list

    def check_doping(self, Na, Nd):
        '''
        checks that Na*Nd is at least bigger than ni^2
        then increases it

        This is rubish
        '''

        # if np.any(Na * Nd) < self.ni**2:
            # if Na > Nd:
                # Nd = self.ni**2 / Na
            # else:
                # Nd = self.ni**2 / Na
        return Na, Nd

    def print_model_notes(self, model=None):
        '''
         prints the notes about the modells
         inputs:
            model: str 
                prints notes on model, if not model is seltecte prints all model notes
        '''

        if model is None:
            models = self.available_models()
        else:
            models = [model]

        for mdl in models:
            print '{0}:\t'.format(mdl),
            try:
                print
                print dict(self.Models.items(mdl))['notes']
            except:
                print 'No notes'
