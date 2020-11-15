#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 12:11:26 2020

@author: janbinkowski
"""
import os
import sys
import subprocess
from datetime import datetime

import click

def gunzip_(path: str):
    querry = ["gzip", "-d", "-k", path]
    subprocess.run(querry, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def check_directory_(path: str) -> bool:
    dir_list = os.listdir(path)
    if not dir_list:
        return False
    else:
        return True
    
def update_read_path_(path: str) -> str:
    return path[:-2]
    
def print_(msg: str, bold=False) -> None:
    click.echo(click.style(msg), bold=bold)
    
def clear_() -> None:
    click.clear()
    
def break_() -> None:
    sys.exit(1)

def create_results_dir_(project_name: dir) -> tuple:
    path = os.path.join("Results/", f"{project_name}")
    
    final_path = os.path.join(path, "Final")
    raw_path = os.path.join(path, "Raw")
    
    try:
        os.makedirs(raw_path)
        os.makedirs(final_path)
    except FileExistsError:
        pass
    
    return raw_path, final_path
    
