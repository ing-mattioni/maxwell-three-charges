/* Reproducible local PDF build with KaTeX/font and reference checks.
 * npm install --no-save playwright-core@1.57.0
 * node paper/build.cjs [output.pdf]
 * Optional: PLAYWRIGHT_MODULE and EDGE_EXECUTABLE select existing installations.
 */
const fs = require('node:fs');
const path = require('node:path');
const {pathToFileURL} = require('node:url');
const {chromium} = require(process.env.PLAYWRIGHT_MODULE || 'playwright-core');

(async () => {
 const output = path.resolve(process.argv[2] || path.join(__dirname,'paper.pdf'));
 const executablePath=process.env.EDGE_EXECUTABLE || (process.platform==='win32'
   ? 'C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe' : undefined);
 const browser=await chromium.launch({headless:true, executablePath});
 try {
  const page=await browser.newPage({viewport:{width:605,height:1123}});
  const failures=[];
  page.on('pageerror',error=>failures.push(String(error)));
  page.on('requestfailed',request=>failures.push(request.url()+': '+request.failure()?.errorText));
  // Every asset required by the current manuscript is local. Reject accidental
  // external dependencies rather than silently printing missing mathematics.
  await page.route(/^https?:/,route=>route.abort());
  await page.goto(pathToFileURL(path.join(__dirname,'paper.html')).href,{waitUntil:'load'});
  await page.evaluate(()=>document.fonts.ready);
  await page.emulateMedia({media:'print'});
  const checks=await page.evaluate(()=>{
   document.body.style.width='160mm';
   const equations=Array.from(document.querySelectorAll('.eq')).map(el=>{
    const bases=Array.from(el.querySelectorAll('.katex-html > .base'));
    const bounds=bases.map(base=>base.getBoundingClientRect());
    const box=el.getBoundingClientRect();
    const label=el.querySelector('.eqn')?.getBoundingClientRect();
    return {id:el.id,overflow:bounds.some(r=>r.left<box.left-.5 || r.right>box.right+.5),
     labelCollision:!!label&&bounds.some(r=>r.right>label.left-.5&&r.left<label.right&&r.bottom>label.top&&r.top<label.bottom)};
   });
   return {katexErrors:document.body.getAttribute('data-katex-errors'),
    unresolved:document.body.getAttribute('data-unresolved'),mathCount:document.querySelectorAll('.katex').length,
    brokenImages:Array.from(document.images).filter(x=>!x.complete||!x.naturalWidth).map(x=>x.src),
    equations,fonts:Array.from(document.fonts).filter(f=>f.status==='error').map(f=>f.family)};
  });
  checks.runtimeErrors=failures;
  const report=output.replace(/\.pdf$/i,'')+'.layout.json';
  fs.writeFileSync(report,JSON.stringify(checks,null,2));
  if(checks.katexErrors!=='0'||checks.unresolved!=='0'||checks.mathCount<100||
     checks.brokenImages.length||checks.fonts.length||failures.length||
     checks.equations.some(e=>e.overflow||e.labelCollision)) {
   throw new Error('PDF validation failed; inspect '+report);
  }
  // Finish rendering before replacing the existing manuscript.
  const pdf=await page.pdf({preferCSSPageSize:true,printBackground:true});
  fs.writeFileSync(output,pdf);
  console.log('PDF built and checked:',output);
 } finally {await browser.close();}
})().catch(error=>{console.error(error);process.exitCode=1});
