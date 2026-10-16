/* Species-level ordinary least squares with an intercept. */
const NSCStatistics = (() => {
 // Two-sided Pearson test: I_(1-r²)((n-2)/2, 1/2).
 function logGamma(z){
  const c=[676.5203681218851,-1259.1392167224028,771.32342877765313,-176.61502916214059,12.507343278686905,-0.13857109526572012,9.984369578019572e-6,1.5056327351493116e-7];
  if(z<.5)return Math.log(Math.PI)-Math.log(Math.sin(Math.PI*z))-logGamma(1-z);
  z-=1;let x=.99999999999980993;for(let i=0;i<c.length;i++)x+=c[i]/(z+i+1);
  const t=z+7.5;return .5*Math.log(2*Math.PI)+(z+.5)*Math.log(t)-t+Math.log(x);
 }
 function betaFraction(a,b,x){
  const tiny=1e-300,safe=v=>Math.abs(v)<tiny?tiny:v;
  let c=1,d=1/safe(1-(a+b)*x/(a+1)),h=d;
  for(let m=1;m<=300;m++){
   let aa=m*(b-m)*x/((a+2*m-1)*(a+2*m));
   d=1/safe(1+aa*d);c=safe(1+aa/c);h*=d*c;
   aa=-(a+m)*(a+b+m)*x/((a+2*m)*(a+2*m+1));
   d=1/safe(1+aa*d);c=safe(1+aa/c);const delta=d*c;h*=delta;
   if(Math.abs(delta-1)<3e-14)return h;
  }
  return NaN;
 }
 function pearsonP(r,n){
  if(n<3||!Number.isFinite(r)||Math.abs(r)>1)return null;
  const x=(1-Math.abs(r))*(1+Math.abs(r)),a=(n-2)/2,b=.5;
  if(x<=0)return 0;if(x>=1)return 1;
  const factor=Math.exp(logGamma(a+b)-logGamma(a)-logGamma(b)+a*Math.log(x)+b*Math.log1p(-x));
  const p=x<(a+1)/(a+b+2)?factor*betaFraction(a,b,x)/a:1-factor*betaFraction(b,a,1-x)/b;
  return Number.isFinite(p)?Math.max(0,Math.min(1,p)):null;
 }
 function linearFit(points) {
  if(points.length<3 || points.some(p=>!Number.isFinite(p.x)||!Number.isFinite(p.y)))return null;
  const n=points.length,mx=points.reduce((s,p)=>s+p.x,0)/n,my=points.reduce((s,p)=>s+p.y,0)/n;
  let xx=0,yy=0,xy=0;
  for(const p of points){const dx=p.x-mx,dy=p.y-my;xx+=dx*dx;yy+=dy*dy;xy+=dx*dy;}
  if(xx===0||yy===0)return null;
  const slope=xy/xx;
  const r=Math.max(-1,Math.min(1,xy/Math.sqrt(xx*yy)));
  return {slope,intercept:my-slope*mx,r,p:pearsonP(r,n)};
 }
 function lengthAdjustedFit(points){
  if(points.length<4||points.some(p=>![p.x,p.y,p.lengthKb].every(Number.isFinite)))return null;
  const n=points.length,mean=k=>points.reduce((s,p)=>s+p[k],0)/n;
  const mx=mean('x'),my=mean('y'),mz=mean('lengthKb');
  let zz=0,xz=0,yz=0,xx=0,yy=0;
  for(const p of points){const x=p.x-mx,y=p.y-my,z=p.lengthKb-mz;zz+=z*z;xz+=x*z;yz+=y*z;xx+=x*x;yy+=y*y;}
  if(zz<=0||xx<=0||yy<=0)return null;
  const residuals=points.map(p=>({x:p.x-mx-xz/zz*(p.lengthKb-mz),y:p.y-my-yz/zz*(p.lengthKb-mz)}));
  const rx=residuals.reduce((s,p)=>s+p.x*p.x,0),ry=residuals.reduce((s,p)=>s+p.y*p.y,0);
  if(rx<=1e-12*xx||ry<=1e-12*yy)return null;
  const fit=linearFit(residuals);if(!fit)return null;
  const lengthSlope=(yz-fit.slope*xz)/zz;
  // Partial-correlation t test has n - 3 degrees of freedom.
  return {...fit,p:pearsonP(fit.r,n-1),intercept:my-fit.slope*mx,lengthSlope,
   modelIntercept:my-fit.slope*mx-lengthSlope*mz,meanLengthKb:mz,df:n-3};
 }
 function logLengthAdjustedFit(points){
  if(points.length<5||points.some(p=>!Number.isFinite(p.lengthKb)||p.lengthKb<=0))return null;
  return lengthAdjustedFit(points.map(p=>({...p,lengthKb:Math.log10(p.lengthKb)})));
 }
 function lengthMassAdjustedFit(points){
  const n=points.length;
  if(n<5||points.some(p=>![p.x,p.y,p.lengthKb,p.logMass].every(Number.isFinite)))return null;
  const mean=k=>points.reduce((a,p)=>a+p[k],0)/n;
  const mx=mean('x'),my=mean('y'),mz=mean('lengthKb'),mw=mean('logMass');
  const z=points.map(p=>p.lengthKb-mz),w=points.map(p=>p.logMass-mw);
  const dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0),zz=dot(z,z),ww=dot(w,w),zw=dot(z,w);
  if(zz<=0||ww<=0)return null;
  const u=w.map((v,i)=>v-zw/zz*z[i]),uu=dot(u,u);
  if(uu<=1e-12*ww)return null;
  const residual=k=>{const v=points.map(p=>p[k]-mean(k)),a=dot(v,z)/zz,b=dot(v,u)/uu;return v.map((x,i)=>x-a*z[i]-b*u[i]);};
  const x=residual('x'),y=residual('y');
  if(dot(x,x)<=1e-12*points.reduce((s,p)=>s+(p.x-mx)**2,0)||dot(y,y)<=1e-12*points.reduce((s,p)=>s+(p.y-my)**2,0))return null;
  const fit=linearFit(x.map((v,i)=>({x:v,y:y[i]})));if(!fit)return null;
  const v=points.map(p=>p.y-my-fit.slope*(p.x-mx)),massSlope=dot(v,u)/uu;
  const lengthSlope=(dot(v,z)-massSlope*zw)/zz;
  return {...fit,p:pearsonP(fit.r,n-2),df:n-4,intercept:my-fit.slope*mx,lengthSlope,massSlope,
   modelIntercept:my-fit.slope*mx-lengthSlope*mz-massSlope*mw};
 }
 return {linearFit,pearsonP,lengthAdjustedFit,logLengthAdjustedFit,lengthMassAdjustedFit};
})();
if(typeof module!=='undefined')module.exports=NSCStatistics;
