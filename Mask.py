'''
.. :module:: Mask
    :synopsis: Contains the mask class

    
| Developed by: Roland Proud <rp43@st-andrews.ac.uk> 
|               Pelagic Ecology Research Group, University of St Andrews
| Maintained by:
| Modification History
| ?/06/2017   Created class   
| 07/07/2017  Added __list_to_range and layer_to_vector methods
|             Updated build_layer to provide layer output as ping-time-vectors
|             or grid
'''


''' Development

Added ping-time-vector output for layers 
e.g. layer value: ping range: sample time range 

Next:  - method to convert vector output(ptr) to layer grid
       - update build_composite to read ping-time-vectors 
      (therefore providing functionality to merge/add layers
      together of any grid size)

'''
## import packages
import numpy as np
import logging
from operator import itemgetter
import itertools

## get logger
log = logging.getLogger(__name__)
 
class Mask(object):
    '''
    Class that stores layer definitions (method and parameter values) and
    builds single layer masks/composite layer masks.
    
    Each layer is defined as either:
        binary: values of the layer are either 1 or 0
        flag:   values take on a range of int values
        cont:   values are continuous
        
    Binary layers can be combined to form composite masks by either:
        1.) Producing a presence/absence mask: binary layers are added together,
            values that are > 0 are set to 1 resulting in a binary image.
        2.) Merge binary layer values together as strings. Each cell will consist
            of a number (base 2) of length equal to the number of merged masks.
            This binary string is then converted to an integer. Therefore, each 
            unique layer combination (e.g. 0101011) is represented by 
            a unique integer (e.g. 43).
        
    '''
    
    def __init__(self):
        '''
        
        Initialize dictionary object to store layer functions and parameters
        '''
        self.mask = {'binary':{},'flag':{},'cont':{}}   
     
    @staticmethod 
    def __list_to_range(L):
        '''
        '''
        l = []
        for k, g in itertools.groupby(enumerate(L), lambda x: x[1]-x[0] ) :
            v = list(map(itemgetter(1), g))
            if len(v) > 2:
                v = [v[0],v[-1]]
            l += [v]
        return l   

    def layer_to_vector(self,layer_dict):
        '''
        '''
        ## Get list of time vectors and points for each ping      
        layer = layer_dict['layer']
        time  = layer_dict['ping_sample_time']   
        nSamples,nPings  = layer.shape
        pings  = []
        values = []
        times  = []  
        time   = np.round(time,10) 
        for val in np.unique(layer):
            if val == 0:
                continue 
            time_vector = []
            ping        = []   
            for p in range(nPings):
                idx         = np.where(layer[:,p] == val)[0]
                if len(idx) > 0:
                    range_lists = self.__list_to_range(idx)
                    tv          = []
                    for v in range_lists:
                        tv += [[time[x,p] for x in v]]
                    ping        += [p]
                    time_vector += [tv]
            
            ## check if there are any pings with the same time vectors and merge them
            ping = np.array(ping) 
            tvs  = list(np.unique(time_vector))
            ## if there is only one tv/one ping then change structure
            try:
                if type(tvs[0]) != list: 
                    tvs = [[tvs]]
                elif type(tvs[0][0]) != list: 
                    tvs = [tvs]
            except IndexError:
                log.warning('Layer is empty')
                tvs = []
            ping_list        = []
            time_vector_uniq = []   
            for tv in tvs:
                idx               = [int(ping[k]) for k,v in enumerate(time_vector) if v == tv]
                ping_list        += [self.__list_to_range(idx)]
                time_vector_uniq += [tv]
            values += [val]
            pings  += [ping_list]
            times  += [time_vector_uniq]
        
        return {'pings':pings,'times':times,'values':values,\
        'pulse_length':layer_dict['pulse_length'],\
        'ping_time':layer_dict['ping_time']}
    
    
    def add_layer(self, f,params):
        '''
            
        stores layer methods and parameters
        '''
        
        ## get layer type
        layer_type = f.__name__.split('_')[0]
        
        ## check layer type
        if layer_type not in self.mask.keys():
            log.error('Function %s not named correctly, bad type',f.__name__)
            raise TypeError
        
        ## get function unique-name:
        name = f.__name__.split('_')[1]
        
        ## check is already exists, if not, add
        if name not in self.mask[layer_type]:
            log.info('Adding new layer: %s',name)
            self.mask[layer_type][name] = {'params':[params],'func':f}
            return [layer_type, name, 0]
        
        ## update function to most recent   
        self.mask[layer_type][name]['func'] = f     
                
        ## check if param already exist
        values = self.mask[layer_type][name]['params']
        if params in values:
            idx = np.where(np.array(values) == params)[0][0]
            log.warning('layer with these parameter values already exists - layer not added')
            return [layer_type, name, idx]
        else:
            log.info('Appending new parameter values to layer: %s',name)
            values.append(params)
            self.mask[layer_type][name]['params'] = values
            return [layer_type, name, len(values) - 1]
            
    def list_layers(self):
        '''
        print list of all layer functions/parameters stored with ID numbers
        '''
        for ltype in self.mask.keys():
            print(ltype + ' layers:')
            for fname in self.mask[ltype].keys():
                print('\t' + fname + ':')
                for k,param in enumerate(self.mask[ltype][fname]['params']):
                    print('\t\t' + 'ID: ' + str(k) + ', Parameter values:',param)
    
    
    def build_layer(self, obj,layer,output = 'grid'):
        '''
        :param obj: echosounder data and parameters
        :type  obj: raw_reader object
        
        :param layer: layer type, name and ID
        :type  layer: [str, str, int]
        
        :param output: 'grid' for layer, 'ptv' for 
                        ping-time-vectors
        :type  output: str
        
        build layer using a single layer function

        '''
        ## get parameters and function
        params     = self.mask[layer[0]][layer[1]]['params'][layer[2]]
        func       = self.mask[layer[0]][layer[1]]['func']
        ## get layer dictionary
        layer_dict = func(obj,params)   
        if output == 'grid':
            return layer_dict['layer']
        elif output == 'ptv':
            return self.layer_to_vector(layer_dict)
        else:
            log.error('unknown layer output type: %s',output)
            raise ValueError
    
    def build_composite_layer(self,obj,layers,output = 'pa'):
        '''
        :param obj: echosounder data and parameters
        :type  obj: raw_reader object
        
        :param layers: A list of layers, each layer consisting of
                       layer type, name and ID
        :type  layers: [[str, str, int],...,]
        
        :param output: specify whether output should be 
                       presence/absence:'pa' or binary notation
                       stored as integers:'bitwise'
        :type  output: str

        build mask by combining multiple binary layers to either 
        produce a presence/absence mask or binary (base2) representation 

        '''
        ## check layers are all binary
        for l in layers:
            if l[0] != 'binary':
                log.error('One or more layers are not binary, \
                composite layer build requires all binary layers')
                raise TypeError
        
        ## For pa output, add all together and return presence/absence (1 or 0) 
        if output == 'pa':
            comp_layer = self.build_layer(obj,layers[0])
            for l in layers[1:]:
                comp_layer = comp_layer + self.build_layer(obj,l)
            comp_layer[comp_layer > 0] = 1
            return comp_layer
        ## For bitwise, merge all binary layers as strings and convert to int 
        elif output == 'bitwise':
            comp_layer = self.build_layer(obj,layers[0]).astype('S1')
            for l in layers[1:]:
                comp_layer = np.core.defchararray.add(comp_layer, self.build_layer(obj,l).astype('S1'))
            ## return integer value
            return np.reshape(np.array([int(x,2) for x in comp_layer.flatten()]),comp_layer.shape)
        

