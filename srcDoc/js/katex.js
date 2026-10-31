// Render the math that pymdownx.arithmatex (generic mode) leaves in the page:
// inline math as <span class="arithmatex">\(...\)</span>, blocks as
// <div class="arithmatex">\[...\]</div>.
document$.subscribe(({ body }) => {
  body.querySelectorAll(".arithmatex").forEach((el) => {
    const tex = el.textContent.trim();
    const display = tex.startsWith("\\[");
    katex.render(tex.slice(2, -2), el, {
      displayMode: display,
      throwOnError: false,
    });
  });
});
