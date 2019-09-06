# -*- coding: utf-8 -*-
'''
Created on Mon Oct 19 17:30:27 2015
Provides some help for nipype workflows 
(to be used within connect statements mostly)
@author: Tsuchida
'''



def getElementsFromTxtList(infile_list):
    '''
    For getting the element from a text file for each txt in the list, returning a list of elements.
    If infile_list is in fact not a list but one file (i.e. string), an out_list with one element is returned
    '''
    import warnings
    out_list = []
    
    if isinstance(infile_list, str):
        infile_l = [infile_list]
    else:
        infile_l = infile_list
        
    for infile in infile_l:
        with open(infile) as f:
            outtxt_list = f.read().splitlines()
        
        if len(outtxt_list) > 1:
            warnings.warn("In getElementsFromTxtList, infile contains more than one value, only the first one is returned")
            
        out_list.append(outtxt_list[0])
        
    return out_list


def getElementFromList(inlist,idx,slc=None):
    '''
    For selecting a particular element or slice from a list 
    within a nipype connect statement.
    If the slice is longer than the list, this returns the list
    Also if the list has just one element, return the list.
    '''
    if not isinstance(inlist,list):
        return inlist
    if not slc:
        outlist=inlist[idx]
    else:
        if slc == -1:
            outlist=inlist[idx:]
        else:
            outlist=inlist[idx:slc]
    return outlist


def getValFromDict(dict_key, d):
    '''
    Give a dict key for a dict (d), return corresponding dict val.
    '''
    return d[dict_key]

   
def prependString(string1, string2):
    '''
    Prepend string2 to string1
    '''
    return(string2+string1)


def createTuple2(item1, item2):
    '''
    makes a Tuple from supplied args
    '''
    return (item1, item2)
    
  
    
        
        
        