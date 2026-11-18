"""Compare ifacekit rotamer categories and chi angles with phenix.rotalyze output (*.rot next to the *.cif files). Development only.
Run: phenix.rotalyze model.cif outliers_only=False > model.rot   for each model, then python compare_rotalyze.py [grids.npz]"""
import sys,re,glob,collections
from ifacekit import rotamer
npz=sys.argv[1] if len(sys.argv)>1 else None
rotamer.use_reference(npz)
tot=collections.Counter(); chi_bad=0; nchk=0
conf=collections.Counter(); bytype=collections.Counter()
for rot in sorted(glob.glob('*.rot')):
    f=rot[:-4]+'.cif'
    ref={}
    for l in open(rot):
        p=l.rstrip('\n').split(':')
        if len(p)<9: continue
        head=p[0].split()
        if len(head)<3 or not re.match(r'^-?\d+$',head[1]): continue
        name=head[2]
        if len(name)==4: continue        # alternate conformer: skip
        try: sc=float(p[2])
        except ValueError: continue
        chis=[float(x) for x in p[3:7] if x]
        if name in ('ASP','PHE','TYR') and len(chis)>1: chis[1]%=180
        if name=='GLU' and len(chis)>2: chis[2]%=180
        ref[(head[0],head[1],name)]=(chis,p[7].lower(),sc)
    _,rows=rotamer.analyze(f)
    mine={tuple(r['residue'].split(':')):r for r in rows}
    n=0;agree=0
    for k,(chis,cat,sc) in ref.items():
        mk=(k[0],k[1],k[2])
        if mk not in mine: continue
        r=mine[mk]; n+=1; nchk+=1
        if len(chis)==len(r['chi']):
            dev=max(min(abs(a-b),360-abs(a-b)) for a,b in zip(chis,r['chi'])) if chis else 0
            if dev>1.5: chi_bad+=1
        conf[(cat,r['category'])]+=1; bytype[(r['table'],cat,r['category'])]+=1
    print(rot[:-4], 'ref',len(ref),'mine',len(mine),'compared',n)
print(dict(conf)); a=sum(v for (x,y),v in conf.items() if x==y); print('agree',a,'of',nchk,f'{100*a/max(nchk,1):.2f}%','chi mismatches',chi_bad)


tabs=sorted({k[0] for k in bytype})
print(f"{'type':7s} {'n':>5s} {'agree%':>7s}  ref->mine disagreements")
for t in tabs:
    n=sum(v for k,v in bytype.items() if k[0]==t); a=sum(v for k,v in bytype.items() if k[0]==t and k[1]==k[2])
    dis={f"{k[1][:3]}>{k[2][:3]}":v for k,v in bytype.items() if k[0]==t and k[1]!=k[2]}
    print(f"{t:7s} {n:5d} {100*a/n:7.1f}  {dis}")
