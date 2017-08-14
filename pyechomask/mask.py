'''
.. :module:: Mask
    :synopsis: Contains the mask class

    
| Developed by: Roland Proud <rp43@st-andrews.ac.uk> 
|               Pelagic Ecology Research Group, University of St Andrews
| Maintained by:
| Modification History
| ?/06/2017   Created class   

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
    Class that stores mask definitions (method and parameter values) and
    builds single masks/composite masks.
    
    Each mask is defined as either:
        binary: values of the mask are either 1 or 0
        flag:   values take on a range of int values
        cont:   values are continuous
        
    Binary masks can be combined to form composite masks by either:
        1.) Producing a presence/absence mask: binary masks are added together,
            values that are > 0 are set to 1 resulting in a binary mask.
        2.) Merge binary mask values together as strings. Each cell will consist
            of a number (base 2) of length equal to the number of merged masks.
            This binary string is then converted to an integer. Therefore, each 
            unique mask combination (e.g. 0101011) is represented by 
            a unique integer (e.g. 43).
        
    '''
    
    def __init__(self):
        '''
        
        Initialize dictionary object to store mask methods and parameters
        '''
        self.mask = {'binary':{},'flag':{},'cont':{}}   
     
    
    def add_mask(self, f,params):
        '''
            
        stores mask methods and parameters
        '''
        
        ## get mask type
        mask_type = f.__name__.split('_')[0]
        
        ## check mask type
        if mask_type not in self.mask.keys():
            log.error('Function %s not named correctly, bad type',f.__name__)
            raise TypeError
        
        ## get function unique-name:
        name = f.__name__.split('_')[1]
        
        ## check is already exists, if not, add
        if name not in self.mask[mask_type]:
            log.info('Adding new mask: %s',name)
            self.mask[mask_type][name] = {'params':[params],'func':f}
            return [mask_type, name, 0]
        
        ## update function to most recent   
        self.mask[mask_type][name]['func'] = f     
                
        ## check if param already exist
        values = self.mask[mask_type][name]['params']
        if params in values:
            idx = np.where(np.array(values) == params)[0][0]
            log.warning('mask with these parameter values already exists - \
                        mask not added')
            return [mask_type, name, idx]
        else:
            log.info('Appending new parameter values to mask: %s',name)
            values.append(params)
            self.mask[mask_type][name]['params'] = values
            return [mask_type, name, len(values) - 1]
            
    def list_masks(self):
        '''
        print list of all mask methods/parameters stored with ID numbers
        '''
        for mtype in self.mask.keys():
            print(mtype + ' masks:')
            for fname in self.mask[mtype].keys():
                print('\t' + fname + ':')
                for k,param in enumerate(self.mask[mtype][fname]['params']):
                    print('\t\t' + 'ID: ' + str(k) + ', Parameter values:',param)
    
    
    def build_mask(self, Sv_dict, mask_id):
        '''
        :param Sv_dict: contains Sv grid, depth, pulse length
        :type  Sv_dict: dictionary
        
        :param mask_id: mask type, name and ID
        :type  mask_id: [str, str, int]
              
        build mask using a single mask function

        '''
        ## get parameters and function
        params     = self.mask[mask_id[0]][mask_id[1]]['params'][mask_id[2]]
        func       = self.mask[mask_id[0]][mask_id[1]]['func']
        ## get mask dictionary
        Sv_dict   = func(Sv_dict,params)   

        return Sv_dict
    
    
    def build_composite_mask(self,Sv_dict,mask_ids,output_mask_id,output = 'pa'):
        '''
        :param Sv_dict: contains Sv grid, depth, pulse length
        :type  Sv_dict: dictionary
        
        :param mask_ids: A list of mask_ids, each mask_id consisting of
                       mask type, name and ID
        :type  mask_ids: [[str, str, int],...,]
        
        :param output_mask_id: mask type, name and ID relating to desired output grid
        :type  output_mask_id: [str, str, int]
        
        :param output: specify whether output should be 
                       presence/absence:'pa' or binary notation
                       stored as integers:'bitwise'
        :type  output: str

        build mask by combining multiple binary masks to either 
        produce a presence/absence mask or binary (base2) representation 

        '''
        ## check masks are all binary
        for m in mask_ids:
            if m[0] != 'binary':
                log.error('One or more masks are not binary, \
                composite mask build requires all binary masks')
                raise TypeError
        
        ## get grid for mask output
        output_dict     = self.build_mask(Sv_dict,output_mask_id)
        out_row,out_col = output_dict['Sv'].shape
        output_mask     = np.zeros((out_row,out_col))
        pulse_length    = output_dict['pulse_length']
        start_depth     = output_dict['start_depth']
        sample_int      = output_dict['sample_int']
        
        ## if output bitwise change to string grid
        if output == 'bitwise':
            output_mask = np.chararray((out_row,out_col))
        
        ## check pulse length and sample resolution
        ## ADD - allow resolution to be lowered
        ## ADD = allow different ping sets
        for m in mask_ids:
            mask_dict = self.build_mask(Sv_dict,m)
            col_idx   = []
            count     = -1
            ## get pings with the same obs parameters and sample int
            for tPLm,tPLo, sdm,sdo, sim,sio in zip(mask_dict['pulse_length'],pulse_length,
                    mask_dict['start_depth'],start_depth, mask_dict['sample_int'],sample_int):
                count += 1
                if tPLm == tPLo and sdm == sdo and sim == sio:
                    col_idx.append(count)
            ## warn if mask not added
            if len(col_idx) == 0:
                log.warning('mask with ID %s does not have matching observation parameters',m)
                log.info('skipping mask %s',m)
                continue
            ## get shape
            row,col                              = mask_dict['mask'].shape
            mask_new                             = np.zeros((out_row,out_col)) 
            mask_new[0:min(row,out_row),col_idx] = mask_dict['mask'][0:min(row,out_row),col_idx]
            ## For pa output, add all together
            if output == 'pa':                            
                output_mask = mask_new + output_mask
            ## For bitwise, merge all binary masks as strings
            elif output == 'bitwise':
                output_mask = np.core.defchararray.add(output_mask, mask_new.astype('S1'))
        
        ## For pa output, add all together and return presence/absence (1 or 0) 
        if output == 'pa':
            output_mask[output_mask > 0] = 1
            output_dict['mask'] = output_mask
            return output_dict
        ## For bitwise, convert to int 
        elif output == 'bitwise':
            ## return integer value
            output_mask = np.reshape(np.array([int(x,2) for x in \
                                       output_mask.flatten()]),output_mask.shape)
            output_dict['mask'] = output_mask
            return output_dict

        

