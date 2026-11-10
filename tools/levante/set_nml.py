#!/usr/bin/env python3
"""set_nml.py FILE group.key=value ... : set existing keys of a fully expanded par file in place.

Each argument sets one key of one group (case-insensitive group and key names); the value is written as given (quote
strings yourself: 'ice_nor.taub_scale_file="scale.nc"'), a trailing comment on the line is kept. A key that is not in
the group stops with an error: the script never adds keys (add them to the par file, or to the template
input/kryos_defaults.nml for a new key).
"""
import re, sys
f = sys.argv[1]; lines = open(f).read().split("\n")
for kv in sys.argv[2:]:
    gk, val = kv.split("=", 1); grp, key = gk.split(".", 1)
    ing = False; done = False
    for i, l in enumerate(lines):
        s = l.strip()
        if s.lower() == "&" + grp.lower(): ing = True; continue
        if ing and s == "/": break
        if ing and re.match(r"\s*" + re.escape(key) + r"\s*=", l, re.I):
            lines[i] = re.sub(r"=\s*[^!]*?(\s*!.*)?$", "= " + val + r"\1" if re.search("!", l) else "= " + val, l, count=1)
            done = True; break
    if not done: sys.exit(f"set_nml: no key {key} in group &{grp} of {f}")
open(f, "w").write("\n".join(lines))
