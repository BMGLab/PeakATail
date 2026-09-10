// MathJax 3 config used by the docs site (loaded before tex-mml-chtml.js
// in mkdocs.yml::extra_javascript).  Enables $...$ inline + $$...$$ block
// math syntax, and re-typesets on Material's instant-navigation page swaps
// so equations don't vanish when you switch pages without a full reload.
window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true,
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex",
  },
};

document$.subscribe(() => {
  MathJax.startup.output.clearCache();
  MathJax.typesetClear();
  MathJax.texReset();
  MathJax.typesetPromise();
});
