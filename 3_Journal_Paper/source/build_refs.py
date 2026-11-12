import re, sys
from refs_data import REFS
src=open(sys.argv[1]).read()
order=[]
for m in re.finditer(r"\\cite\{([^}]*)\}",src):
    for k in m.group(1).split(","):
        k=k.strip()
        if k not in order: order.append(k)
miss=[k for k in order if k not in REFS]; assert not miss, miss
out=["\\section*{References}","{\\small\\begin{enumerate}[leftmargin=2.4em,labelsep=0.4em,itemsep=2pt,parsep=0pt,label={\\arabic*.}]"]
for i,k in enumerate(order,1):
    t=REFS[k].replace("&","\\&") if "\\&" not in REFS[k] else REFS[k]
    t=re.sub(r"(https?://\S+)",lambda m:"\\url{"+m.group(1)+"}",t)
    out.append(f"\\item {t}")
out.append("\\end{enumerate}}")
open("refs.tex","w").write("\n".join(out)+"\n")
# numbered \cite -> [n]
cites={k:i for i,k in enumerate(order,1)}
def rep(m):
    nums=sorted(cites[k.strip()] for k in m.group(1).split(","))
    # collapse runs
    parts=[];i=0
    while i<len(nums):
        j=i
        while j+1<len(nums) and nums[j+1]==nums[j]+1: j+=1
        parts.append(str(nums[i]) if j==i else (f"{nums[i]},{nums[j]}" if j==i+1 else f"{nums[i]}--{nums[j]}"))
        i=j+1
    return "["+",".join(parts)+"]"
open("manuscript_cited.tex","w").write(re.sub(r"\\cite\{([^}]*)\}",rep,src))
print(len(order),"references")
