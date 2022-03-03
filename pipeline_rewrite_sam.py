import re

r = re.compile(r"(\w+)\*(\d+)-(\d)\t")

def replaceSam(sam, out):
    f = open(sam)
    f = map(lambda i: re.sub(r, r"\1*\2\t", i), f)
    a = set()
    f = filter(lambda i: i not in a and not a.add(i), f)
    open(out, "w").writelines(f)


replaceSam("/home/linnil1/kir/PING/PING_test/linnil1_syn_wide/linnil1_syn_wide.00.read..sam",
           "/home/linnil1/kir/PING/PING_test/linnil1_syn_wide/linnil1_syn_wide.00.read.sam")

