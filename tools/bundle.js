/* Node stand-in for the gitignored bundle.ps1: inline every LOCAL css/js into
   one self-contained file, strip ?v= as it goes, leave Google Fonts external. */
const fs = require('fs');
let html = fs.readFileSync('index.html', 'utf8');
const inlined = [];
html = html.replace(/[ \t]*<link\b[^>]*rel="stylesheet"[^>]*>/g, tag => {
  const m = /href="([^"]+)"/.exec(tag);
  if (!m || /^https?:/.test(m[1])) return tag;
  const file = m[1].split('?')[0];
  inlined.push(file);
  return '<style data-inlined-from="' + file + '">\n' + fs.readFileSync(file, 'utf8').replace(/\s*$/, '') + '\n</style>';
});
/* A deferred script runs after parsing and after every parser-inserted script,
   so inlining it WHERE IT SITS would move it earlier and change the order the
   app boots in. Deferred ones are held back and emitted before </body>. */
const deferred = [];
html = html.replace(/[ \t]*<script\b([^>]*)\bsrc="([^"]+)"([^>]*)><\/script>/g, (tag, a, src, b) => {
  if (/^https?:/.test(src)) return tag;
  const file = src.split('?')[0];
  inlined.push(file);
  const block = '<script data-inlined-from="' + file + '">\n'
    + fs.readFileSync(file, 'utf8').replace(/\s*$/, '') + '\n</script>';
  if (/\bdefer\b/.test(a + b)){ deferred.push(block); return ''; }
  return block;
});
if (deferred.length){
  /* The LAST </body>: app.js carries one inside a template literal that builds
     an export document, and matching that one injects the app into the export. */
  const at = html.lastIndexOf('</body>');
  if (at === -1) throw new Error('no </body> to place the deferred scripts before');
  html = html.slice(0, at) + deferred.join('\n') + '\n' + html.slice(at);
}
fs.writeFileSync('Psychometric_Assistant.html', html);
console.log('inlined:', inlined.join(', '));
