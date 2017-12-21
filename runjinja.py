#!/usr/bin/env python

"""
Preprocesses foo.c.in to foo.c. Reads STDIN and writes STDOUT.
"""

import sys
import hashlib
from jinja2 import Template, Environment

env = Environment(block_start_string='/*{',
                  block_end_string='}*/',
                  variable_start_string='{{',
                  variable_end_string='}}')

extra_vars = dict(len=len)
input = sys.stdin.read()
sys.stdout.write('/* DO NOT EDIT. md5sum of source: %s */' % hashlib.md5(input.encode()).hexdigest())
sys.stdout.write(env.from_string(input).render(**extra_vars))
