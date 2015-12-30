
import matplotlib.pylab as plt
import json
import sys


"""
To do:

    Need to fix the check_doping function
"""

class HelperFunctions():

    def change_model(self, model_author=None):

        if model_author is None:
            self.model_author = self.Models.get('default', 'model')
        else:
            # Need a check to make sure craNhcan't be passed
            self.model_author = model_author
        # print self.model_author
        self.model = self.Models.get(self.model_author, 'model')

        self.vals = dict(self.Models.items(self.model_author))
        # List.remove('model')
        del self.vals['model']

        for k, v in self.vals.iteritems():

            # checks if float or list
            try:
                self.vals[k] = float(v)
           
            except:
                try:
                    self.vals[k] =[float(i) for i in  v.split(';')]
                    # print self.vals[k]
                except:
                    pass
        pass

    def plot_all_models(self, update_function, xvalues=None, **kwargs):
        '''
        cycles through all the models and plots the result
        inputs:
            update_function: str
                 the name of the specific function used to update the model_author
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

    def available_models(self, Filter=None, Filter_value=None):
        '''
        Returns a list of all models.
        If this value is provide the list returned will consists of models that 
        have the value in field (indicated by Filter) indicated by "Filter_value"

        inputs: (all optional)
            Filter: (str)
                The model field that is to be checked        
            Filter_value:
                The value to be checked for
        '''

        model_list = self.Models.sections()
        model_list.remove('default')

        # does the filtering
        if Filter is not None:
            for c, model in enumerate(model_list):
                if self.Models.get(model, Filter) not in Filter_value:
                    model_list.remove(model)

        # prints no models available
        if not model_list:
            print 'No models available'

        return model_list

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
