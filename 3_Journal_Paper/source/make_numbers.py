# Builds numbers.tex (macros) and table fragments from the result tables written by Smoking_Preterm_Full_Analysis.R
import pandas as pd, numpy as np, re, os
T="results/tables/"; OUT="gen/"; os.makedirs(OUT,exist_ok=True)
import os as _os
def R(f):
    p=T+f
    if not _os.path.exists(p): p="../run/results/tables/"+f; print("fallback to earlier run for",f)
    return pd.read_csv(p)
def esc(s):
    s=str(s)
    for a,b in (("\\","\\textbackslash{}"),("&","\\&"),("%","\\%"),("_","\\_"),("#","\\#")): s=s.replace(a,b)
    s=s.replace("<","$<$").replace(">","$>$")
    return s
f2=lambda x:f"{x:.2f}"; f1=lambda x:f"{x:.1f}"; f3=lambda x:f"{x:.3f}"; f0=lambda x:f"{x:,.0f}"
sg=lambda x:f"{x:+.2f}".replace("+","+")
M={}
def m(k,v): M[k]=str(v)
# ---- cohort, EDA
flow=R("cohort_flow.csv"); eda=R("eda_numbers.csv").set_index("quantity").value
m("NExtract",f0(flow.n[0])); m("NResident",f0(flow.n[1])); m("NTotal",f0(flow.n[2])); m("PctKept",f2(flow.pct_of_extract[2])); m("NExcluded",f0(flow.n[1]-flow.n[2]))
m("NSmokers",f0(eda["Smokers"])); m("PrevSmoke",f2(eda["Smoking prevalence (%)"])); m("PretermAll",f2(eda["Preterm, all (%)"]))
m("PretermSmk",f1(eda["Preterm, smokers (%)"])); m("PretermNon",f1(eda["Preterm, non-smokers (%)"]))
m("PretermSmkTwo",f2(eda["Preterm, smokers (%)"])); m("PretermNonTwo",f2(eda["Preterm, non-smokers (%)"]))
m("GAsmk",f2(eda["Mean gestational age, smokers (weeks)"])); m("GAnon",f2(eda["Mean gestational age, non-smokers (weeks)"]))
cc=R("complete_case_comparison.csv").set_index("kept"); m("DroppedPreterm",f1(cc.loc[False,"pct_preterm"])); m("KeptPreterm",f1(cc.loc[True,"pct_preterm"]))
m("DroppedSmoke",f2(cc.loc[False,"pct_smokers"])); m("KeptSmoke",f2(cc.loc[True,"pct_smokers"])); m("NDropped",f0(cc.loc[False,"births"]))
# ---- main
mn=R("results_main.csv").set_index("estimator")
a=mn.loc["AIPTW (doubly robust)"]; cr=mn.loc["Crude (unadjusted)"]; att=mn.loc["AIPTW, effect among smokers (ATT)"]; ato=mn.loc["Overlap weights (ATO)"]; ip=mn.loc["IPTW (stabilized)"]; gc=mn.loc["Outcome regression (g-computation)"]
m("CrudeRD",f2(cr.RD)); m("CrudeLo",f2(cr.lo)); m("CrudeHi",f2(cr.hi)); m("RDmain",f2(a.RD)); m("RDlo",f2(a.lo)); m("RDhi",f2(a.hi)); m("SEmain",f3(a.SE))
m("RRmain",f2(a.RR)); m("RRcrude",f2(cr.RR)); m("RDg",f2(gc.RD)); m("RDiptw",f2(ip.RD)); m("RDiptwLo",f2(ip.lo)); m("RDiptwHi",f2(ip.hi))
m("ATT",f2(att.RD)); m("ATTlo",f2(att.lo)); m("ATThi",f2(att.hi)); m("ATO",f2(ato.RD)); m("ATOlo",f2(ato.lo)); m("ATOhi",f2(ato.hi))
m("ReducePct",f0(100*(cr.RD-a.RD)/cr.RD)); m("PerK",f0(10*att.RD))
ev=R("evalues.csv"); m("RRlo",f2(ev.risk_ratio[1])); m("EV",f2(ev.E_value[0])); m("EVlo",f2(ev.E_value[1]))
r0=a.RD/(a.RR-1); m("RRhi",f2((r0+a.hi)/r0))
# ORs
risk1=a.RD  # not used
# ---- spec
rb=R("robustness_summary.csv").set_index("analysis")
sp={"Pre":"Pre-pregnancy set only (no care timing, no WIC)","Ins":"Main plus insurance type","Int":"Main plus two-way interactions","Five":"Clamp at 5 / 95 percentile"}
for k,v in sp.items():
    m("RD"+k,f2(rb.loc[v,"RD"])); m("RD"+k+"Lo",f2(rb.loc[v,"lo"])); m("RD"+k+"Hi",f2(rb.loc[v,"hi"]))
adj=rb.RD.iloc[1:]; m("RangeMin",f2(adj.min())); m("RangeMax",f2(adj.max())); m("PreShift",f2(rb.loc[sp["Pre"],"RD"]-a.RD))
ps=R("sens_ps_handling.csv").set_index("analysis")
nc=ps.loc["No clamp"]; m("RDnoclamp",f2(nc.RD)); tr=[i for i in ps.index if i.startswith("Delete")][0]; m("RDtrim",f2(ps.loc[tr,"RD"])); m("TrimKeep",re.search(r"keeps ([\d.]+)%",tr).group(1))
m("MaxWtClamp",f0(ps.loc["Clamp at 1 / 99 percentile","max_weight_stabilized"])); m("MaxWtNone",f0(nc.max_weight_stabilized))
m("ESSclamp",f0(ps.loc["Clamp at 1 / 99 percentile","effective_n_smokers"])); m("ESSnone",f0(nc.effective_n_smokers))
# ---- overlap / balance
b=R("balance_table.csv"); k=len(b)
m("NLevels",k); m("BalBefore",int((b.before.abs()>.1).sum())); m("BalAfter",int((b.after_iptw.abs()>.1).sum())); m("BalOv",int((b.after_overlap.abs()>.1).sum()))
m("SMDbefore",f2(b.before.abs().max())); m("SMDafter",f2(b.after_iptw.abs().max())); m("SMDov",f2(b.after_overlap.abs().max()))
cs=R("cell_scores.csv"); m("LowNon",f0(100*cs.n0[cs.e<0.01].sum()/cs.n0.sum())); m("LowSmk",f1(100*cs.n1[cs.e<0.01].sum()/cs.n1.sum()))
who=R("who_has_no_lookalike.csv"); ed=who[who.variable=="educ"].set_index("level")
m("LowBach",f0(ed.loc["Bachelor's","pct_of_low_group"])); m("LowMast",f0(ed.loc["Master's","pct_of_low_group"])); m("LowDoc",f0(ed.loc["Doctorate/Prof.","pct_of_low_group"]))
m("LowBachSum",f0(ed.loc[["Bachelor's","Master's","Doctorate/Prof."],"pct_of_low_group"].sum()))
# ---- variance
v=R("variance_check.csv").set_index("quantity").value
m("SEif",f3(v["Standard error, influence function"])); m("SEboot",f3(v["Standard error, full-sample bootstrap"])); m("SEratio",f2(v["Ratio: bootstrap SE / influence-function SE"]))
m("BootLo",f2(v["95% CI, bootstrap percentile (lower)"])); m("BootHi",f2(v["95% CI, bootstrap percentile (upper)"])); m("BootB",f0(v["Bootstrap resamples used"])); m("BootMean",f2(v["Bootstrap mean of the estimate"]))
# ---- crossfit
cf=R("crossfit_summary.csv"); c1=cf[cf.clamp=="1%/99%"].set_index("learner")
for key,nm in (("Log","Logistic_main"),("Lasso","Lasso_interactions"),("Forest","Random_forest"),("Avg","Average_of_3")):
    m("CF"+key,f2(c1.loc[nm,"RD"])); m("CF"+key+"Lo",f2(c1.loc[nm,"lo"])); m("CF"+key+"Hi",f2(c1.loc[nm,"hi"]))
    m("CFw"+key,f2(c1.loc[nm,"RD_weighting_only"])); m("CFo"+key,f2(c1.loc[nm,"RD_outcome_only"]))
cfr=R("crossfit_by_repeat.csv"); m("CFreps",cfr.rep.nunique())
m("CFrange",f2(cf.range_over_repeats.max()))
# ---- history
h=R("history_intensity_results.csv"); hc=R("history_crude.csv").set_index("hist")
A=h[h.contrast.str.startswith("A. ")].iloc[0]; A2=h[h.contrast.str.startswith("A2.")].iloc[0]; Q=h[h.contrast.str.startswith("B2.")].iloc[0]
m("ContQuit",f2(A.RD)); m("ContQuitLo",f2(A.lo)); m("ContQuitHi",f2(A.hi)); m("ContQuitPre",f2(A2.RD))
m("QuitNever",f2(Q.ATT)); m("QuitNeverLo",f2(Q.ATT_lo)); m("QuitNeverHi",f2(Q.ATT_hi)); m("QuitNeverATE",f2(Q.RD))
m("HistNever",f1(hc.loc["Never","preterm_pct"])); m("HistQuit",f1(hc.loc["Quit","preterm_pct"])); m("HistCont",f1(hc.loc["Continued","preterm_pct"])); m("HistStart",f1(hc.loc["Started in pregnancy","preterm_pct"]))
m("NQuit",f0(hc.loc["Quit","births"])); m("NCont",f0(hc.loc["Continued","births"])); m("NStart",f0(hc.loc["Started in pregnancy","births"])); m("NNever",f0(hc.loc["Never","births"]))
Bd=h[h.contrast.str.startswith("B. Continued")].reset_index(drop=True)
m("DoseMin",f2(Bd.ATT.min())); m("DoseMax",f2(Bd.ATT.max()))
# ---- benchmark and ranking
cb=R("confounder_contributions.csv").set_index("variable")
m("EducShare",f0(cb.loc["educ","share_of_smoking_model_pct"])); m("PPShare",f0(cb.loc["ppterm","share_of_preterm_model_pct"]))
m("EducChange",f2(cb.loc["educ","change_in_rd"])); m("PPChange",f2(cb.loc["ppterm","change_in_rd"])); m("PrecareChange",f2(cb.loc["precare","change_in_rd"])); m("WicChange",f2(cb.loc["wic","change_in_rd"]))
m("PrecareShareOut",f0(cb.loc["precare","share_of_preterm_model_pct"]))
ph=R("public_health_translation.csv").set_index("quantity").value
m("ExcessPerK",f0(ph.iloc[0])); m("ExcessBirths",f0(round(ph.iloc[2],-2))); m("OneIn",f0(ph.iloc[3])); m("ExcessBirthsExact",f0(ph.iloc[2]))
m("ExcessLo",f0(round(att.lo/100*eda["Smokers"],-2))); m("ExcessHi",f0(round(att.hi/100*eda["Smokers"],-2)))
# table-1 facts
t1=R("table1_by_smoking.csv"); pp=t1[(t1.variable=="ppterm")&(t1.level=="Yes")].iloc[0]; pn=t1[(t1.variable=="ppterm")&(t1.level=="No")].iloc[0]
m("PPsmk",f1(pp.pct_of_smokers)); m("PPnon",f1(pp.pct_of_nonsmokers)); m("PPrisk",f1(pp.preterm_pct_nonsmoker)); m("NoPPrisk",f1(pn.preterm_pct_nonsmoker))
L=R("levels_by_confounder.csv"); e8=L[L.variable=="educ"].set_index("level"); m("SmkNoDip",f1(e8.loc["9-12th, no diploma","pct_smokers"])); m("SmkDoc",f1(e8.loc["Doctorate/Prof.","pct_smokers"]))
cv=R("cramers_v.csv").set_index("var"); m("VeducSmk",f2(cv.loc["educ","smoker"])); m("VppPre",f2(cv.loc["ppterm","preterm"])); m("VeducWic",f2(cv.loc["wic","educ"]))
open(OUT+"numbers.tex","w").write("\n".join(f"\\newcommand{{\\{k}}}{{{v}}}" for k,v in M.items())+"\n"); print(len(M),"macros")

# ------------------------------------------------------------------ table fragments
def table(name,cols,hdr,rows,align=None,small=True,long=False):
    align=align or ("l"+"r"*(cols-1))
    body="\n".join(" & ".join(r)+" \\\\" for r in rows)
    if long:
        s=f"\\begin{{longtable}}{{{align}}}\n\\toprule\n{hdr} \\\\\n\\midrule\n\\endhead\n{body}\n\\bottomrule\n\\end{{longtable}}\n"
    else:
        s=f"\\begin{{tabular}}{{{align}}}\n\\toprule\n{hdr} \\\\\n\\midrule\n{body}\n\\bottomrule\n\\end{{tabular}}\n"
    open(OUT+name,"w").write(("{\\small\n" if small else "")+s+("}\n" if small else ""))
ci=lambda r,lo="lo",hi="hi": f"{r[lo]:.2f} to {r[hi]:.2f}"
# flow
table("tab_flow.tex",3,"Step & Births & \\% of extract",[[esc(r.step),f0(r.n),f2(r.pct_of_extract)] for r in flow.itertuples()],"lrr")
# dictionary
dd=R("dictionary_analysis_variables.csv").fillna("")
rows=[[f"\\texttt{{{esc(r.variable)}}}",esc(r.type),"" if r.n_levels=="" else str(int(r.n_levels)),esc(r.levels_or_range),f2(r.pct_missing)] for r in dd.itertuples()]
table("tab_dict.tex",5,"Variable & Stored as & Levels & Levels or range & \\% missing",rows,"llrp{7.2cm}r",long=True)
# levels (appendix)
rows=[[f"\\texttt{{{esc(r.variable)}}}",esc(r.level),f0(r.n),f0(r.n_smokers),f2(r.pct_smokers),f2(r.pct_preterm)] for r in L.itertuples()]
table("tab_levels.tex",6,"Variable & Level & Births & Smokers & \\% who smoked & \\% preterm",rows,"llrrrr",long=True)
# main
order=["Crude (unadjusted)","Outcome regression (g-computation)","IPTW (stabilized)","AIPTW (doubly robust)","AIPTW, effect among smokers (ATT)","Overlap weights (ATO)"]
rows=[]
for e in order:
    r=mn.loc[e]; rows.append([esc(e),f2(r.RD),"" if pd.isna(r.lo) else f"{r.lo:.2f} to {r.hi:.2f}","" if pd.isna(r.RR) else f2(r.RR)])
table("tab_main.tex",4,"Estimator & Risk difference (points) & 95\\% CI & Risk ratio",rows,"lrrr")
# robustness
rows=[[esc(i),f2(r.RD),f"{r.lo:.2f} to {r.hi:.2f}"] for i,r in rb.iterrows()]
table("tab_spec.tex",3,"Analysis & Risk difference (points) & 95\\% CI",rows,"lrr")
# ps handling
rows=[[esc(i),f2(r.RD),f"{r.lo:.2f} to {r.hi:.2f}","" if pd.isna(r.max_weight_stabilized) else f0(r.max_weight_stabilized),"" if pd.isna(r.effective_n_smokers) else f0(r.effective_n_smokers)] for i,r in ps.iterrows()]
table("tab_ps.tex",5,"Handling of the propensity score & Risk difference & 95\\% CI & Largest weight & Effective n (smokers)",rows,"p{6.2cm}rrrr")
# balance (largest 12 before)
bt=b.reindex(b.before.abs().sort_values(ascending=False).index).head(12)
rows=[[esc(r.covariate),f2(abs(r.before)),f2(abs(r.after_iptw)),f2(abs(r.after_overlap))] for r in bt.itertuples()]
table("tab_balance.tex",4,"Covariate level (12 most imbalanced) & Before & After weighting & After overlap weights",rows,"lrrr")
# balance full (appendix)
bt=b.reindex(b.before.abs().sort_values(ascending=False).index)
rows=[[esc(r.covariate),f3(r.before),f3(r.after_iptw),f3(r.after_overlap)] for r in bt.itertuples()]
table("tab_balance_all.tex",4,"Covariate level & SMD before & SMD after weighting & SMD after overlap weights",rows,"lrrr",long=True)
# weights
wt=R("weight_summary.csv")
rows=[[esc(r.group),f0(r.births),f2(r.max_weight),f2(r.mean_weight),f0(r.effective_n)] for r in wt.itertuples()]
table("tab_weights.tex",5,"Group & Births & Largest weight & Mean weight & Effective sample size",rows,"lrrrr")
# who has no lookalike
w2=who[who.variable=="educ"]
rows=[[esc(r.level),f1(r.pct_of_low_group),f1(r.pct_of_rest)] for r in w2.itertuples()]
table("tab_who.tex",3,"Education & \\% of non-smokers with score $<$ 0.01 & \\% of all other non-smokers",rows,"lrr")
# variance
vv=R("variance_check.csv")
rows=[[esc(r.quantity),f"{r.value:.3f}" if r.value<100 else f0(r.value)] for r in vv.itertuples()]
table("tab_var.tex",2,"Quantity & Value",rows,"lr")
# crossfit
rows=[]
for nm,lab in (("Logistic_main","Logistic, main effects"),("Lasso_interactions","LASSO with all two-way interactions"),("Random_forest","Random forest"),("Average_of_3","Average of the three")):
    r=c1.loc[nm]; rows.append([lab,f2(r.RD),f"{r.lo:.2f} to {r.hi:.2f}",f2(r.RD_outcome_only),f2(r.RD_weighting_only),f0(r.max_weight),f0(r.ESS)])
table("tab_cf.tex",7,"Learner & Risk difference & 95\\% CI & Outcome part & Weighting part & Largest weight & Effective n",rows,"lrrrrrr")
rows=[]
for nm,lab in (("Logistic_main","Logistic"),("Lasso_interactions","LASSO"),("Random_forest","Forest"),("Average_of_3","Average")):
    r=c1.loc[nm]; rows.append([lab,f2(r.pred_risk_smokers),f2(r.obs_risk_smokers),f2(r.pred_risk_nonsmokers),f2(r.obs_risk_nonsmokers),f2(r.pred_smoking_prob),f2(r.obs_smoking_prob)])
table("tab_cfcal.tex",7,"Learner & Smokers, predicted & observed & Non-smokers, predicted & observed & Smoking, predicted & observed",rows,"lrrrrrr")
# history
lab={"A. Continued vs quit (also adjusted for cigarettes/day before)":"Continued vs quit","A2. Continued vs quit, pre-pregnancy confounders only":"Continued vs quit, pre-pregnancy confounders only","B2. Quit vs never-smokers":"Quit vs never-smokers"}
rows=[]
for r in h.itertuples():
    nm=lab.get(r.contrast) or r.contrast.replace("B. Continued, ","Continued, ").replace(" cigarettes/day before pregnancy, vs never"," cigarettes per day, vs never")
    rows.append([esc(nm),f0(r.n_group),f2(r.crude_RD),f2(r.RD),f"{r.lo:.2f} to {r.hi:.2f}",f2(r.ATT),f"{r.ATT_lo:.2f} to {r.ATT_hi:.2f}"])
table("tab_hist.tex",7,"Comparison & Group size & Crude & ATE & 95\\% CI & ATT & 95\\% CI",rows,"lrrrrrr")
hc2=R("history_crude.csv"); rows=[[esc(r.hist),f0(r.births),f2(r.preterm_pct)] for r in hc2.itertuples()]
table("tab_histcrude.tex",3,"Smoking history & Births & \\% preterm",rows,"lrr")
# ---- trimester of smoking (full natality file, Section 10)
tc=R("trimester_crude.csv").set_index("pattern"); tr=R("trimester_results.csv").set_index("contrast")
m("TriN",f0(tc.births.sum())); m("TriNever",f2(tc.loc["Never smoked","preterm_pct"])); m("TriFirst",f2(tc.loc["First trimester only","preterm_pct"]))
m("TriFirstSecond",f2(tc.loc["First and second trimesters only","preterm_pct"])); m("TriAll",f2(tc.loc["All three trimesters","preterm_pct"])); m("TriLate",f2(tc.loc["Started after the first trimester","preterm_pct"]))
m("TriLateN",f0(tc.loc["Started after the first trimester","births"]))
def tri(k,nm):
    r=tr.loc[k]; m(nm+"ATT",f2(r.ATT)); m(nm+"Lo",f2(r.ATT_lo)); m(nm+"Hi",f2(r.ATT_hi))
tri("First trimester only vs never","TriA"); tri("First and second trimesters only vs never","TriB"); tri("All three trimesters vs never","TriC"); tri("Started after the first trimester vs never","TriD")
dk=[k for k in tr.index if "cigarettes/day in 3rd" in k]
for i,k in enumerate(dk): m("TriDose"+"ABCD"[i],f2(tr.loc[k,"ATT"])); m("TriDoseLo"+"ABCD"[i],f2(tr.loc[k,"ATT_lo"])); m("TriDoseHi"+"ABCD"[i],f2(tr.loc[k,"ATT_hi"]))
m("TriDoseWideATE",f2(tr.loc[dk[3],"RD"])); m("TriDoseWideLo",f2(tr.loc[dk[3],"lo"])); m("TriDoseWideHi",f2(tr.loc[dk[3],"hi"]))
def tl(k):
    return (k.replace(", vs never","").replace(" vs never","").replace("Smoked all trimesters, ","All three trimesters, ").replace(" cigarettes/day in 3rd trimester,"," cigarettes per day in third trimester"))
rows=[[esc(tl(r.Index)),f0(r.n_group),f2(r.crude_RD),f2(r.RD),f"{r.lo:.2f} to {r.hi:.2f}",f2(r.ATT),f"{r.ATT_lo:.2f} to {r.ATT_hi:.2f}"] for r in tr.itertuples()]
table("tab_tri.tex",7,"Smoking pattern, each against never-smokers & Group size & Crude & ATE & 95\\% CI & ATT & 95\\% CI",rows,"lrrrrrr")
# evalues and benchmark
rows=[[esc(r.estimate),f2(r.risk_ratio),f2(r.E_value)] for r in ev.itertuples()]
table("tab_ev.tex",3,"Estimate & Risk ratio & E-value",rows,"lrr")
bm=R("evalue_benchmark.csv")
rows=[[esc(r.covariate),f1(r.pct_of_cohort),f2(r.RR_with_smoking),f2(r.RR_with_preterm),"yes" if r.exceeds_E_value_on_both else "no",sg(r.change_in_rd_when_left_out)] for r in bm.itertuples()]
table("tab_bench.tex",6,"Measured covariate & \\% of cohort & Risk ratio with smoking & Risk ratio with preterm birth (non-smokers) & Above E-value on both & Change in estimate if left out",rows,"lrrrcr")
# contributions
nm={"educ":"Education","race":"Race and ethnicity","age_grp":"Maternal age","precare":"Prenatal-care timing","bmi_grp":"Pre-pregnancy BMI","wic":"WIC","parity":"Parity","ppterm":"Previous preterm birth","phype":"Pre-pregnancy hypertension","pdiab":"Pre-pregnancy diabetes"}
cbo=cb.reset_index().sort_values("share_of_smoking_model_pct",ascending=False)
rows=[[nm[r.variable],f1(r.share_of_smoking_model_pct),f1(r.share_of_preterm_model_pct),f2(r.rd_without_variable),sg(r.change_in_rd)] for r in cbo.itertuples()]
table("tab_contrib.tex",5,"Confounder & Share of smoking model (\\%) & Share of preterm model (\\%) & Estimate without it & Change",rows,"lrrrr")
# public health
rows=[[esc(r.quantity),f1(r.value) if r.value<1000 else f0(r.value)] for r in ph.reset_index().itertuples()]
table("tab_ph.tex",2,"Quantity & Value",rows,"lr")
# missingness
ms=R("missingness.csv"); rows=[[f"\\texttt{{{esc(r.variable)}}}",f0(r.n_missing),f2(r.pct_missing)] for r in ms.itertuples()]
table("tab_miss.tex",3,"Confounder & Births missing & \\% missing",rows,"lrr")
rows=[["Kept (complete data)",f0(cc.loc[True,"births"]),f2(cc.loc[True,"pct_smokers"]),f2(cc.loc[True,"pct_preterm"])],["Dropped (any confounder missing)",f0(cc.loc[False,"births"]),f2(cc.loc[False,"pct_smokers"]),f2(cc.loc[False,"pct_preterm"])]]
table("tab_cc.tex",4,"Group & Births & \\% smoked & \\% preterm",rows,"lrrr")
print("tables written")
_t=open(OUT+"tab_tri.tex").read(); open(OUT+"tab_tri.tex","w").write(re.sub(r"(?<=[ (])-(?=\d)","$-$",_t))

for k in list(M):
    if k.startswith("Tri") and M[k].startswith("-"): M[k]="$-$"+M[k][1:]
open(OUT+"numbers.tex","w").write("\n".join(f"\\newcommand{{\\{k}}}{{{v}}}" for k,v in M.items())+"\n"); print(len(M),"macros (final)")
