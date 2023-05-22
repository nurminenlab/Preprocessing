#!/usr/bin/env python3
# coding=utf-8
import os
for basename in os.listdir('.'):
    if not basename.endswith('.svg'):
        continue
    print(basename)
    with open(basename) as fd:
        s = fd.read()
    with open(basename, 'w') as fd:
        fd.write(s.replace('stroke-miterlimit:100000;', ''))
