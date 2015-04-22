#!/usr/bin/env python

"""
Preprocesses foo.c.in to foo.c. Reads STDIN and writes STDOUT.
"""

import sys
from jinja2 import Template, Environment

env = Environment(block_start_string='/*{',
                  block_end_string='}*/',
                  variable_start_string='{{',
                  variable_end_string='}}')

extra_vars = dict(len=len)
sys.stdout.write(env.from_string(sys.stdin.read()).render(**extra_vars))
