#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 10:40:08 2020

@author: janbinkowski
"""
import os
import subprocess as sub

class Shell:
    def __init__(self, main_directory: str):
        self.main_dir = main_directory
        self.cwd = None
        
    @staticmethod
    def prepare_command_(command: str) -> list:
        command = command.split()        
        
        return command
    
    def set_path(self, directory) -> None:
        self.cwd = os.path.join(self.main_dir, directory)
    
    def reset_cwd(self) -> None:
        self.cwd = self.main_dir
        
    def run(self, command: list) -> str:
        process = sub.run(command, cwd=self.cwd, stdout=sub.PIPE)
        return process.stdout