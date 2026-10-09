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
 return {linearFit,pearsonP};
})();
if(typeof module!=='undefined')module.exports=NSCStatistics;
