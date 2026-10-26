import csv, json
from base import *
A = build(True); B = build(False)
def ranges(ks):
    ks = sorted(ks); out = []; i = 0
    while i < len(ks):
        j = i
        while j + 1 < len(ks) and ks[j + 1] == ks[j] + 1: j += 1
        out.append(str(ks[i]) if i == j else f"{ks[i]}--{ks[j]}"); i = j + 1
    return ", ".join(out)
rows_csv = []; summary = {}; newlist = []
tex = []
for p in PRIMES:
    T = A[p]
    for t in range(2, TMAX + 1):
        lines = []; ne = na = nn = nno = 0
        for v in range(t + 2, VMAX + 1):
            ks = [k for k in range(t + 1, v // 2 + 1) if adm(p, t, k, v) and v - k >= t]
            if not ks: continue
            ex = [k for k in ks if T.has(t, k, v) and (t, k, v) in B[p].E]
            nw = [k for k in ks if T.has(t, k, v) and (t, k, v) not in B[p].E]
            no = [k for k in ks if (t, k, v) in T.NO]
            op = [k for k in ks if k not in ex and k not in nw and k not in no]
            na += len(ks); ne += len(ex) + len(nw); nn += len(nw); nno += len(no)
            for k in ks:
                st = "exists" if k in ex or k in nw else ("no" if k in no else "open")
                rows_csv.append((p, t, k, v, st, T.E.get((t, k, v), T.NO.get((t, k, v), ""))))
            for k in nw: newlist.append((p, t, k, v, T.E[(t, k, v)]))
            def runs(ks, bold):
                out = []
                for r in ranges(ks).split(", "):
                    if r: out.append((int(r.split("--")[0]), "\\textbf{" + r + "}" if bold else r))
                return out
            exs = ", ".join(x for _, x in sorted(runs(ex, False) + runs(nw, True)))
            lines.append(f"{v} & {exs} & {ranges(op)} \\\\")
        summary[(p, t)] = (na, ne, nn, nno)
        if na == 0: continue
        if ne == na - nno:
            tex.append(f"\\subsection*{{$\\mathrm{{LS}}[{p}]({t},k,v)$}}\nAll {na} admissible parameter sets with $v\\le{VMAX}$ and $k\\le v/2$ are realised" + (f", except the {nno} that are known not to exist" if nno else "") + ".\n")
            if nno:
                tex[-1] += "Non-existent: " + ", ".join(f"$({t},{k},{v})$" for (pp, tt, k, v, st, w) in rows_csv if pp == p and tt == t and st == "no") + ".\n"
            continue
        if ne == 0 and nno == 0:
            tex.append(f"\\subsection*{{$\\mathrm{{LS}}[{p}]({t},k,v)$}}\nNone of the {na} admissible parameter sets with $v\\le{VMAX}$ and $k\\le v/2$ is known to be realised.\n"); continue
        tex.append(f"\\subsection*{{$\\mathrm{{LS}}[{p}]({t},k,v)$}}\n{ne} of {na} admissible parameter sets exist" + (f", {nno} do not exist" if nno else "") + f", {na-ne-nno} are open.\n" + (("Non-existent: " + ", ".join(f"$({t},{k},{v})$" for (pp, tt, k, v, st, w) in rows_csv if pp == p and tt == t and st == "no") + ".\n") if nno else "") + "\n"
                   "{\\footnotesize\\begin{longtable}{r>{\\raggedright\\arraybackslash}p{7.0cm}>{\\raggedright\\arraybackslash}p{7.0cm}}\n\\toprule $v$ & exists, $k=$ & open, $k=$\\\\\\midrule\\endhead\n" + "\n".join(lines) + "\n\\bottomrule\\end{longtable}}\n")
open("tables_body.tex", "w").write("\n".join(tex))
with open("large_sets_table.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["N", "t", "k", "v", "status", "justification"]); w.writerows(rows_csv)
with open("summary.tex", "w") as f:
    for p in PRIMES:
        f.write(f"{p} & " + " & ".join(f"{summary[(p,t)][1]}/{summary[(p,t)][0]}" for t in range(2, TMAX + 1)) + " \\\\\n")
json.dump(newlist, open("newlist.json", "w"))
print(len(rows_csv), sum(1 for r in rows_csv if r[4] == "exists"), sum(1 for r in rows_csv if r[4] == "open"), sum(1 for r in rows_csv if r[4] == "no"))
for n in newlist: print(n)
print([r for r in rows_csv if r[4]=="no"])
