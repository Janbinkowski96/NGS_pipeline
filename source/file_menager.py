#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 08:56:35 2020

@author: janbinkowski
"""
import os

class FileMenager:
    def __init__(self, main_directory: str):
        self.main_directory = main_directory
    
    @staticmethod
    def join_paths(path_outer: str, path_inner: str) -> str:
        return os.path.join(path_outer, path_inner)
    
    @staticmethod
    def create_directory(path=str):
        try:
            os.mkdir(path)
            return 0
        except FileExistsError:
            return 1
        
    