# selamibeyhan

from setuptools import setup, Extension

setup(
    name = 'ChebyshevFeatureSelection',
    ext_modules = [Extension('ChebyshevFeatureSelection', ['aisa_feature.c'])])
