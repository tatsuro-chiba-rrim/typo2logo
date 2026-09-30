const assert=require('node:assert/strict'),{solve,chains}=require('../kirlian-model.js');
function fixture(nx,ny,fn){const mask=Uint8Array.from({length:nx*ny},(_,i)=>+fn(i%nx,Math.floor(i/nx)));return{nx,ny,mask};}
let start=performance.now();
const ring=fixture(96,64,(x,y)=>Math.hypot(x-48,y-32)>=11&&Math.hypot(x-48,y-32)<=15);
const a=solve({...ring,density:.35,length:70,seed:8});
assert(a.stats.maxResidual<1e-4);assert.equal(a.phi[32*96+48],1);assert(!a.nodes.some(n=>Math.hypot(n.x-48,n.y-32)<11),'closed hole stays discharge-free');
for(let y=0;y<a.ny;y++)for(let x=0;x<a.nx;x++){const i=y*a.nx+x;if(!x||!y||x===a.nx-1||y===a.ny-1)assert.equal(a.phi[i],0);if(a.original[i])assert.equal(a.phi[i],1);}
assert.equal(new Set(a.nodes.map(n=>n.cell)).size,a.nodes.length);for(let i=0;i<a.nodes.length;i++){const n=a.nodes[i];if(n.parent>=0){assert(n.parent<i);const p=a.nodes[n.parent];assert(Math.max(Math.abs(n.x-p.x),Math.abs(n.y-p.y))===1);}assert(!a.original[n.cell]);}assert(chains(a).every(p=>p.points.length>=2));
assert.deepEqual(solve({...ring,density:.35,length:70,seed:8}).nodes,a.nodes);assert.equal(solve({...ring,voltage:0}).nodes.length,0);
const two=fixture(128,80,(x,y)=>Math.hypot(x-49,y-40)<=13||Math.hypot(x-79,y-40)<=13),b=solve({...two,density:0});
const E=(m,x,y)=>(1-m.phi[y*m.nx+x])*(m.nx-1);
const outward=E(b,35,40),gap=E(b,63,40);assert(outward>gap*2,'same-potential narrow gap is shielded');
const tip=fixture(128,96,(x,y)=>(x>=37&&x<=90&&y>=48&&y<=65)||(y>=25&&y<48&&Math.abs(x-64)<=Math.floor((y-25)*.27))),t=solve({...tip,density:0});
const pointed=E(t,64,24),flat=E(t,48,47);assert(pointed>flat*1.4,'protrusion concentrates field');
console.log(JSON.stringify({seconds:(performance.now()-start)/1000,ring:a.stats,gap,outward,pointed,flat,tests:'boundary conditions, residual, hole shielding, proximity shielding, tip enhancement, deterministic forest, zero voltage PASS'},null,2));
const wall=solve({...fixture(96,64,(x,y)=>Math.hypot(x-48,y-18)<=10),density:0});
assert(E(wall,48,7)>E(wall,48,29)*1.4,'nearby grounded wall strengthens field');
// Independently inspect all rendered, sampled curve segments, including collinear overlap.
const curves=require('../kirlian-model.js').curves;
for(let seed=1;seed<=5;seed++){
 const m=solve({...ring,density:.18,length:90,seed}),paths=curves(m).paths,segs=[];
 paths.forEach((p,id)=>{for(let k=1;k<p.length;k++)segs.push({a:p[k-1],b:p[k],id,k});});
 const cross=(a,b,c)=>(b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]);const same=(a,b)=>Math.hypot(a[0]-b[0],a[1]-b[1])<1e-8;
 for(let i=0;i<segs.length;i++)for(let j=i+1;j<segs.length;j++){const s=segs[i],t=segs[j];if(s.id===t.id&&Math.abs(s.k-t.k)<=1)continue;const {a,b}=s,{a:c,b:d}=t;const v=[cross(a,b,c),cross(a,b,d),cross(c,d,a),cross(c,d,b)];
 if(v.every(x=>Math.abs(x)<1e-8)){const axis=Math.abs(a[0]-b[0])>Math.abs(a[1]-b[1])?0:1;const overlap=Math.min(Math.max(a[axis],b[axis]),Math.max(c[axis],d[axis]))-Math.max(Math.min(a[axis],b[axis]),Math.min(c[axis],d[axis]));assert(overlap<1e-8,'overlapping branches');}
 else if(!same(a,c)&&!same(a,d)&&!same(b,c)&&!same(b,d))assert(!(v[0]*v[1]<-1e-10&&v[2]*v[3]<-1e-10),'crossing branches');}
}
console.log('Grounded wall proximity and rendered curve intersection checks: PASS');
