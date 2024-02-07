Compilation:

make xmaplelcm

Or

make xmaplelcm_debug (if we want assertions to be checked)


RUN: ./xmaplelcm -compute-dip -dip-2clauses -dip-pair-min=5 -produce-proof -dip-type=1 file.cnf

HELP: ./xmaplelcm --help-verb
display extended information about command-line options
