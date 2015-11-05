#!/usr/bin/env python

# Author: Ni Huang <nihuang at genetics dot wustl dot edu>

import collections
import functools
from numpy import ndarray

class memoized(object):
    '''Decorator. Caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned
    (not reevaluated).
    '''
    def __init__(self, func):
        self.func = func
        self.cache = {}
    def __call__(self, *args):
        if not isinstance(args, collections.Hashable):
            # uncacheable. a list, for instance.
            # better to not cache than blow up.
            return self.func(*args)
        #------------------------------------
        # special treatment for numpy.ndarray
        if type(args[0]) == ndarray:
            key = tuple(args[0])
            if key in self.cache:
                return self.cache[key]
            else:
                value = self.func(*args)
                self.cache[key] = value
                return value
        #------------------------------------
        if args in self.cache:
            return self.cache[args]
        else:
            value = self.func(*args)
            self.cache[args] = value
            return value
    def __repr__(self):
        '''Return the function's docstring.'''
        return self.func.__doc__
    def __get__(self, obj, objtype):
        '''Support instance methods.'''
        return functools.partial(self.__call__, obj)
