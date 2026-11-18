"""Compare ifacekit Ramachandran categories with phenix.ramalyze output (*.ram next to the *.cif files). Development only."""
import sys,re,glob,collections
from ifacekit import rama
for ram in sorted(glob.glob('*.ram')):
    f=ram[:-4]+'.cif'
    ref={}
    for l in open(ram):
        m=re.match(r'\s*(\S+)\s+(-?\d+\S*)\s+([A-Z]{3}):([\d.]+):(-?[\d.]+):(-?[\d.]+):(\w+):(.*)',l)
        if m: ref[(m.group(1),m.group(2),m.group(3))]=(float(m.group(5)),float(m.group(6)),m.group(7),m.group(8))
    _,rows=rama.analyze(f)
    mine={}
    for r in rows:
        c,n,nm=r['residue'].split(':'); mine[(c,n,nm)]=r
    common=set(ref)&set(mine); diff=[]; cat=collections.Counter()
    for k in common:
        a=ref[k];b=mine[k]
        cat[(a[2].lower(),b['category'])]+=1
        if a[2].lower()!=b['category'] or a[3].lower().replace('isoleucine or valine','ile/val')[:4]!=rama.LABEL[b['class']].lower()[:4]: diff.append((k,a,b['category'],b['class'],b['phi'],b['psi'],b['density']))
    print(f, 'ref',len(ref),'mine',len(mine),'common',len(common)); print('  ',dict(cat))
    for d in diff[:6]: print('   DIFF',d)
