#!/usr/bin/env python
# file: grep.py
import re, sys

map(sys.stdout.write,(l for l in sys.stdin if re.search(sys.argv[1],l)))
