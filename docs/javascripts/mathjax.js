window.MathJax = {
    tex: {
      inlineMath: [["$","$"], ["\\(","\\)"]],
      displayMath: [["$$","$$"], ["\\[", "\\]"]],
      processEscapes: true,
      processEnvironments: true
    },
    options: {
      ignoreHtmlClass: ".*|",
      processHtmlClass: "arithmatex"
    }
};
  
document$.subscribe(() => {
    MathJax.typesetPromise()
})

MathJax.Hub.Config({
    "tex2jax": { 
        inlineMath: [["$","$"], ["\\(","\\)"]],
        displayMath: [["$$","$$"], ["\\[", "\\]"]],
    }
});
MathJax.Hub.Config({
    config: ["MMLorHTML.js"],
    jax: ["input/TeX", "output/HTML-CSS", "output/NativeMML"],
    extensions: ["MathMenu.js", "MathZoom.js"]
});