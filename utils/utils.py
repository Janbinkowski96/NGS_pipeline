#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 12:11:26 2020

@author: janbinkowski
"""
import sys
import click

    
def print_(msg: str, bold=False, blink=False, is_command=False) -> None:
    if is_command:
        msg = " ".join(map(lambda char: str(char), msg))
        msg = "Running: " + msg
    click.echo(click.style(f"{msg}\n", bold=bold, blink=blink))
    
def clear_() -> None:
    click.clear()
    
def break_(msg: str) -> None:
    print_(msg, blink=True)
    sys.exit(1)

def save_raw_file(data: str, path: str, file_name: str) -> None:
    with open(f"{path}/{file_name}.txt", "wb") as handle:
        handle.write(data)

def help_(ctx, param, value):
    if not value or not param:
        return
    click.echo(ctx.get_help())
    ctx.exit()


def check_flags_(ctx, read_1, read_2, reference, reference_genome, regions, algorithm, output, project_name):
    if not read_1 or not read_1 or not reference or not reference_genome  or not regions or not algorithm or not project_name:
        arguments = locals().items()
        missing_args = [flag for flag, value in arguments if value is None]
        missing_args = " | ".join(missing_args)
        click.echo(click.style(f"Missing option/s: {missing_args}, use --help / -h to print options.", bold=True))
        ctx.exit()