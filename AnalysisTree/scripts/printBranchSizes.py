#! /usr/bin/env python

import sys
import ROOT as r

if len(sys.argv) < 2:
    print "  usage: {0} <input root file> <tree name (\"Events\")>".format(sys.argv[0])
    exit(0)

fin = r.TFile(sys.argv[1])
t = fin.Get("Events" if len(sys.argv)<3 else sys.argv[2])

bs = t.GetListOfBranches()

sizes = []
for b in bs:
    sizes.append((b.GetName(), b.GetTotalSize(), b.GetZipBytes(), (1.0*b.GetTotalSize()/b.GetZipBytes() if b.GetZipBytes()>0 else 0)))

max_len = len(sorted(sizes, key=lambda x:len(x[0]))[-1][0])
total_size = sum(s[2] for s in sizes)

sizes = sorted(sizes, key=lambda x:x[2], reverse=True)

fmt_string =  "{{0:{0}s}} {{1:>12s}} {{2:>12s}} {{3:>12s}} {{4:>9s}}".format(max_len)
print fmt_string.format("Branch","Raw Size", "Zip Size", "Compression", "% Total")
fmt_string =  "{{0:{0}s}} {{1:12d}} {{2:12d}} {{3:12.3f}} {{4:8.2f}}%".format(max_len)
for s in sizes:
    print fmt_string.format(s[0], s[1], s[2], s[3], 100.0*s[2]/total_size)


