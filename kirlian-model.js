/* 2D, quasi-static dielectric-breakdown-inspired growth.
 * Dirichlet: original glyphs and grown channels = 1, outer wall = 0.
 * Four-neighbour Laplace stencil; SOR; eight-direction, thresholded growth.
 * Field re-solved every 8 occupied cells. No gas/space-charge/time model.
 */
const KirlianModel=(()=>{
function solve({nx,ny,mask,density=1,length=55,voltage=1,eta=1.6,threshold=.8,seed=52909}){
 const n=nx*ny,phi=new Float64Array(n),fixed=new Uint8Array(n),original=new Uint8Array(n),occupied=new Int32Array(n);occupied.fill(-1);phi.fill(1);
 const dirs=[-1,1,-nx,nx],growthDirs=[-1,1,-nx,nx,-nx-1,-nx+1,nx-1,nx+1],nodes=[],frontier=new Set(),diagonals=new Set();let sweeps=0,solves=0,maxResidual=0;
 for(let y=0;y<ny;y++)for(let x=0;x<nx;x++){const i=y*nx+x;if(x===0||y===0||x===nx-1||y===ny-1){fixed[i]=2;phi[i]=0;}else if(mask[i]){fixed[i]=1;original[i]=1;}}
 // Same-potential closed counters have constant potential and no external field.
 const reachable=new Uint8Array(n),queue=new Int32Array(n);let head=0,tail=0;
 for(let i=0;i<n;i++)if(fixed[i]===2){reachable[i]=1;queue[tail++]=i;}
 while(head<tail){const i=queue[head++],x=i%nx,y=(i/nx)|0;for(const d of dirs){if(d===-1&&x===0||d===1&&x===nx-1||d===-nx&&y===0||d===nx&&y===ny-1)continue;const j=i+d;if(!reachable[j]&&!original[j]){reachable[j]=1;queue[tail++]=j;}}}
 const free=[];for(let i=0;i<n;i++){if(!fixed[i]&&!reachable[i])fixed[i]=3;if(!fixed[i])free.push(i);}
 function relax(maxSweeps=650){let delta=Infinity,k=0;const omega=1.88,tolerance=1e-5;
   while(delta>tolerance&&k<maxSweeps){delta=0;for(let z=0;z<free.length;z++){const i=free[z];if(fixed[i])continue;const correction=omega*((phi[i-1]+phi[i+1]+phi[i-nx]+phi[i+nx])*.25-phi[i]);phi[i]+=correction;delta=Math.max(delta,Math.abs(correction));}k++;}
   let residual=0;for(const i of free)if(!fixed[i])residual=Math.max(residual,Math.abs(phi[i]-(phi[i-1]+phi[i+1]+phi[i-nx]+phi[i+nx])*.25));sweeps+=k;solves++;maxResidual=Math.max(maxResidual,residual);return residual;
 }
 const initialResidual=relax(1500),initialField=phi.slice();
 for(const i of free)if(growthDirs.some(d=>original[i+d]))frontier.add(i);
 let a=seed>>>0;const random=()=>{a=(Math.imul(a,1664525)+1013904223)>>>0;return a/4294967296;};
 const maxDepth=Math.max(2,Math.round(length*(nx-1)/1000)),budget=Math.round(900*density*(.6+.4*length/55));
 function edgeAllowed(i,j){const x=i%nx,y=(i/nx)|0,xx=j%nx,yy=(j/nx)|0;if(x!==xx&&y!==yy){if(original[y*nx+xx]||original[yy*nx+x])return false;const sign=(x-xx)*(y-yy)>0?1:-1;if(diagonals.has(`${Math.min(x,xx)},${Math.min(y,yy)},${-sign}`))return false;}return true;}
 function parentOf(i){let electrode=-1,hasChannel=false;const candidates=[];for(const d of growthDirs){const j=i+d;if(original[j]&&edgeAllowed(i,j))electrode=j;if(occupied[j]>=0){hasChannel=true;const parent=occupied[j],distance=Math.abs(d)===1||Math.abs(d)===nx?1:Math.SQRT2;if(nodes[parent].depth+distance<=maxDepth&&edgeAllowed(i,j))candidates.push({parent,electrode:-1,distance});}}
   if(hasChannel&&electrode>=0)return null;
   if(hasChannel){candidates.sort((a,b)=>a.distance-b.distance||a.parent-b.parent);return candidates[0]||null;}
   return electrode>=0?{parent:-1,electrode,distance:Math.abs(i-electrode)===1||Math.abs(i-electrode)===nx?1:Math.SQRT2}:null;
 }
 if(voltage>0)for(let step=0;step<budget;step++){
   let total=0;const candidates=[],weights=[],parents=[];
   for(const i of frontier){if(fixed[i]){frontier.delete(i);continue;}const par=parentOf(i);if(!par){frontier.delete(i);continue;}const field=Math.max(0,(1-phi[i])*voltage*(nx-1)-threshold),weight=Math.pow(field/par.distance,eta);if(weight>1e-14){candidates.push(i);weights.push(weight);parents.push(par);total+=weight;}}
   if(!total)break;let pick=random()*total,idx=0;for(;idx<weights.length-1;idx++){pick-=weights[idx];if(pick<=0)break;}
   const i=candidates[idx],par=parents[idx],node={cell:i,x:i%nx,y:Math.floor(i/nx),parent:par.parent,rootCell:par.electrode,depth:par.parent<0?par.distance:nodes[par.parent].depth+par.distance};occupied[i]=nodes.length;nodes.push(node);fixed[i]=1;phi[i]=1;frontier.delete(i);
   const parentCell=par.parent<0?par.electrode:nodes[par.parent].cell,px=parentCell%nx,py=Math.floor(parentCell/nx);if(node.x!==px&&node.y!==py)diagonals.add(`${Math.min(node.x,px)},${Math.min(node.y,py)},${(node.x-px)*(node.y-py)>0?1:-1}`);
   for(const d of growthDirs){const j=i+d;if(j>=0&&j<n&&!fixed[j])frontier.add(j);}
   if((step+1)%8===0)relax();
 }
 const finalResidual=relax(1500);
 return {nx,ny,phi,initialField,original,fixed,nodes,stats:{cells:nodes.length,solves,sweeps,initialResidual,finalResidual,maxResidual,voltage,eta,threshold,maxDepth}};
}
// Extract a forest into chains between roots, forks and tips.
function chains(model){const {nodes,nx}=model,children=nodes.map(()=>[]),out=[];nodes.forEach((p,i)=>{if(p.parent>=0)children[p.parent].push(i);});
 function walk(start,first,initial){const pts=initial||[[nodes[start].x,nodes[start].y]];let at=first;while(true){pts.push([nodes[at].x,nodes[at].y]);if(children[at].length!==1)break;at=children[at][0];}out.push({points:pts,depth:nodes[at].depth});}
 nodes.forEach((p,i)=>{if(p.parent<0){const rx=p.rootCell%nx,ry=Math.floor(p.rootCell/nx);walk(i,i,[[(rx+p.x)/2,(ry+p.y)/2]]);}else if(children[p.parent].length!==1)walk(p.parent,i);});return out;
}
function curves(model){
 const raw=chains(model).map(p=>p.points),eps=1e-8;
 const same=(a,b)=>Math.abs(a[0]-b[0])<eps&&Math.abs(a[1]-b[1])<eps;
 function intersects(a,b,c,d){
   const cross=(p,q,r)=>(q[0]-p[0])*(r[1]-p[1])-(q[1]-p[1])*(r[0]-p[0]);
   const ab1=cross(a,b,c),ab2=cross(a,b,d),cd1=cross(c,d,a),cd2=cross(c,d,b);
   if(Math.abs(ab1)<eps&&Math.abs(ab2)<eps){const axis=Math.abs(a[0]-b[0])>Math.abs(a[1]-b[1])?0:1;return Math.min(Math.max(a[axis],b[axis]),Math.max(c[axis],d[axis]))-Math.max(Math.min(a[axis],b[axis]),Math.min(c[axis],d[axis]))>eps;}
   if(same(a,c)||same(a,d)||same(b,c)||same(b,d))return false;
   return ab1*ab2<=eps&&cd1*cd2<=eps;
 }
 function segments(paths){const out=[];paths.forEach((p,id)=>{for(let i=1;i<p.length;i++)out.push({a:p[i-1],b:p[i],id,ord:i});});return out;}
 function keys(s){const out=[];for(let y=Math.floor(Math.min(s.a[1],s.b[1]));y<=Math.floor(Math.max(s.a[1],s.b[1]));y++)for(let x=Math.floor(Math.min(s.a[0],s.b[0]));x<=Math.floor(Math.max(s.a[0],s.b[0]));x++)out.push(x+','+y);return out;}
 function index(segs){const map=new Map();for(const seg of segs)for(const k of keys(seg)){if(!map.has(k))map.set(k,[]);map.get(k).push(seg);}return map;}
 function collisions(segs,map,skipSame){const bad=new Set();for(const s of segs){const seen=new Set();for(const key of keys(s))for(const t of map.get(key)||[]){if(seen.has(t)||s===t||skipSame&&s.id===t.id||s.id===t.id&&Math.abs(s.ord-t.ord)<=1)continue;seen.add(t);if(intersects(s.a,s.b,t.a,t.b)){bad.add(s.id);if(!skipSame)bad.add(t.id);}}}return bad;}
 const smooth=raw.map(pts=>{if(pts.length<3)return pts;const out=[pts[0]];let start=pts[0];for(let i=1;i<pts.length-1;i++){const control=pts[i],end=[(control[0]+pts[i+1][0])/2,(control[1]+pts[i+1][1])/2];for(let k=1;k<=8;k++){const t=k/8,u=1-t;out.push([u*u*start[0]+2*u*t*control[0]+t*t*end[0],u*u*start[1]+2*u*t*control[1]+t*t*end[1]]);}start=end;}out.push(pts[pts.length-1]);return out;});
 const bad=collisions(segments(smooth),index(segments(raw)),true),paths=smooth.map((p,i)=>bad.has(i)?raw[i]:p);
 const more=collisions(segments(paths),index(segments(paths)),false);for(const i of more){paths[i]=raw[i];bad.add(i);}
 // Curves rejected above revert to the certified non-crossing lattice forest.
 const remaining=collisions(segments(paths),index(segments(paths)),false);if(remaining.size)return {paths:raw,rounded:0,fallback:raw.length};
 return {paths,rounded:raw.length-bad.size,fallback:bad.size};
}

return {solve,chains,curves};
})();
if(typeof module!=='undefined'&&module.exports)module.exports=KirlianModel;
