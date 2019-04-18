def fix_TR_fs(in_file, TR_sec):
    # do not use mri_convert interface as it has a bug (already fixed in niyppe master)
    import os
    TR_ms = TR_sec * 1000.
    out_file = os.path.abspath(os.path.basename(in_file))
    cmd = 'mri_convert -tr {TR} {in_file} {out_file}'.format(TR=TR_ms, in_file=in_file, out_file=out_file)
    with open('command.txt', 'w') as fi:
        fi.write(cmd)
    os.system(cmd)
    return out_file