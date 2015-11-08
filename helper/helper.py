
import matplotlib.pylab as plt 

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
            try:
                self.vals[k] = float(v)
            except:
                pass

    def plot_all_models(self, update_function, **kwargs):
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
        for model in self.AvailableModels():
            print update_function, model
        # ax.plot(np.inf,np.inf,'k-',label = 'Auger')
            self.change_model(model)
            result = getattr(self, update_function)(**kwargs)
            ax.plot(result, label=model)
        plt.legend(loc=0)

    def AvailableModels(self):
        '''
        opens the models const file and obtains the available models
        '''
        model_list = self.Models.sections()
        model_list.remove('default')
        return model_list

    def check_doping(self, Na, Nd):
        '''
        checks that Na*Nd is at least bigger than ni^2
        then increases it
        '''

        if Na * Nd < self.ni**2:
            if Na > Nd:
                Nd = self.ni**2 / Na
            else:
                Nd = self.ni**2 / Na
        return Na, Nd